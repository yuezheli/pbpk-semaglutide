# date: 4/18/2026 
# author: Yuezhe Li 
# purpose of this code: to simulate for different dosing

using ProjectRoot
using Pkg; Pkg.activate(@projectroot);
using ModelingToolkit
using ModelingToolkit: getdefault
using DifferentialEquations
using DataFrames
using Plots
using Plots.PlotMeasures
using SymbolicIndexingInterface
using CSV
using Sundials 

## observed PK
# semaglutide, sc, 0.5 mg
# https://pubmed.ncbi.nlm.nih.gov/30788808/, Figure 2
plasma_sc = DataFrame(
    time_d = [0.25, 0.51, 0.66, 1.0, 1.28, 1.5, 1.78, 2.0, 2.27, 2.52, 2.99, 4.01, 5.01, 6.01, 7.0],
    conc_nM = [8.88, 10.73, 11.4, 12.1, 12.2, 11.99, 11.88, 11.88, 11.77, 11.67, 11.27, 10.87, 9.87, 9.12, 8.08]
);

# semaglutide, oral, 10 mg x10
plasma_oral = DataFrame(
    time_d = [1.06, 3.06, 5.05, 7.05, 8.06, 9.12, 10.05, 11.07, 13.05, 16.09, 20.02, 23.07, 30.05],
    conc_nM = [1.59, 4.7, 8.5, 10.79, 11.68, 11.61, 13.58, 11.55, 9.2, 6.73, 4.31, 3.36, 1.52]
);

## Setup
# load constants 
include(@projectroot("script/constants.jl"));

# load model 
include(@projectroot("model/twopore_albumin_peptide.jl"));
@named pbpk = albumin_peptide_pbpk_mtk();
pbpk_sys = structural_simplify(pbpk);

# Set the time span 
tspan = (0.0, 800.0);

# update parameters from mouse to human 
include(@projectroot("script/helper-volume-mouse2human.jl"));
p_homo = updatevolume_human(pbpk_sys, kon_Alb = 1E6, koff_Alb = 3.1); 
prob_homo = ODEProblem(pbpk_sys, p_homo, tspan);

## Dosing  
include(@projectroot("script/helper-infusion.jl"));
# sc
u0_sc = Dict([pbpk_sys.A_SC => 1E-3/MW_semaglutide]); # mol
p_abs_sc = Dict([pbpk_sys.bioavailability => 0.84, pbpk_sys.k_a => 0.0253]);  # https://pmc.ncbi.nlm.nih.gov/articles/PMC6437231/
prob_sc = remake(prob_homo, u0=u0_sc, p=p_abs_sc);
affect_sc!(integrator) = integrator[pbpk_sys.A_SC] += 1E-3/MW_semaglutide;
cb_sc = PeriodicCallback(affect_sc!, hr_per_day * 7.0); 

# oral
u0_oral = Dict([pbpk_sys.A_SC => 10E-3/MW_semaglutide]); # mol 
p_abs_oral = Dict([pbpk_sys.bioavailability => 0.01, pbpk_sys.k_a => 0.0286]);
prob_oral = remake(prob_homo, u0=u0_oral, p = p_abs_oral);
oral_dose_times = hr_per_day:hr_per_day:(hr_per_day * 9);
affect_oral!(integrator) = integrator[pbpk_sys.A_SC] += 10E-3/MW_semaglutide;
cb_oral = PresetTimeCallback(oral_dose_times, affect_oral!);

## Simulation
sol_sc = solve(prob_sc, Rodas5(), reltol=1e-6, abstol=1e-9, saveat = 1.0, callback = cb_sc);
sol_oral = solve(prob_oral, Rodas5(), reltol=1e-6, abstol=1e-9, saveat = 1.0, callback = cb_oral);

## Visualization
plt_sc = Plots.plot(xlabel="Time (hr)", ylabel="Plasma semaglutude concentration (nM)", 
         yaxis = :log10, yticks = [1E-2, 1E-1, 1, 10, 1E2, 1E3], ylims = [1E-2, 1E3], 
         xticks = [0, 24, 48, 168, 336, 504, 672], xlims = [0, 672]);
Plots.plot!(sol_sc.t, (sol_sc[pbpk.C_Plasma_Albud] + sol_sc[pbpk.C_Plasma_Exo]) * nmol_per_mol, label = "sims");
Plots.plot!(plasma_sc.time_d * hr_per_day, plasma_sc.conc_nM, seriestype = :scatter, label = "Overgaard et al., 2019"); 
display(plt_sc)

savefig(plt_sc, @projectroot("deliv/figure/plasma-pk-semaglutide-sc.png"));

plt_oral = Plots.plot(xlabel="Time (hr)", ylabel="Plasma semaglutude concentration (nM)", 
         yaxis = :log10, yticks = [1E-2, 1E-1, 1, 10, 1E2, 1E3], ylims = [1E-2, 1E3], 
         xticks = [0, 24, 48, 72, 168, 336, 504, 672], xlims = [0, 672]);
Plots.plot!(sol_oral.t, (sol_oral[pbpk.C_Plasma_Albud] + sol_oral[pbpk.C_Plasma_Exo]) * nmol_per_mol, label = "sims");
Plots.plot!(plasma_oral.time_d * hr_per_day, plasma_oral.conc_nM, seriestype = :scatter, label = "Overgaard et al., 2021"); 
display(plt_oral)

savefig(plt_oral, @projectroot("deliv/figure/plasma-pk-semaglutide-oral.png"));

## Additional comparison 
# set up dosing of 1.5 mg oral everyday, vs sc 0.25 mg 
u0_sc_point25 = Dict([pbpk_sys.A_SC => 0.25E-3/MW_semaglutide]); # mol
affect_sc_point25!(integrator) = integrator[pbpk_sys.A_SC] += 0.25E-3/MW_semaglutide;
cb_sc_point25 = PeriodicCallback(affect_sc_point25!, hr_per_day * 7.0); 

u0_oral_1point5 = Dict([pbpk_sys.A_SC => 1.5E-3/MW_semaglutide]); # mol 
affect_oral_1point5!(integrator) = integrator[pbpk_sys.A_SC] += 1.5E-3/MW_semaglutide;
cb_oral_1point5 = PeriodicCallback(affect_oral_1point5!, hr_per_day); 

sol_sc_point25 = solve(remake(prob_sc, u0 = u0_sc_point25), Rodas5(), reltol=1e-6, abstol=1e-9, saveat = 1.0, callback = cb_sc_point25);
sol_oral_1point5 = solve(remake(prob_oral, u0 = u0_oral_1point5), Rodas5(), reltol=1e-6, abstol=1e-9, saveat = 1.0, callback = cb_oral_1point5);

plt_oral_sc = Plots.plot(xlabel="Time (hr)", ylabel="Semaglutude concentration (nM)", 
         yaxis = :log10, yticks = [1E-2, 1E-1, 1, 10, 1E2, 1E3], ylims = [1E-2, 1E3], 
         xticks = [0, 24, 48, 72, 168, 336, 504, 672], xlims = [0, 800]);
Plots.plot!(sol_sc_point25.t, (sol_sc_point25[pbpk.C_Plasma_Albud] + sol_sc_point25[pbpk.C_Plasma_Exo]) * nmol_per_mol, label = "sc, 0.25 mg, plasma", color = :red);
Plots.plot!(sol_oral_1point5.t, (sol_oral_1point5[pbpk.C_Plasma_Albud] + sol_oral_1point5[pbpk.C_Plasma_Exo]) * nmol_per_mol, label = "oral, 1.5 mg, plasma", color = :blue);
Plots.plot!(sol_sc_point25.t, (sol_sc_point25[pbpk.C_IS_Albud_SI] + sol_sc_point25[pbpk.C_IS_Exo_SI]) * nmol_per_mol, label = "sc, 0.25 mg, small interstine", color = :red, linestyle = :dash);
Plots.plot!(sol_oral_1point5.t, (sol_oral_1point5[pbpk.C_IS_Albud_SI] + sol_oral_1point5[pbpk.C_IS_Exo_SI]) * nmol_per_mol, label = "oral, 1.5 mg, small interstine", color = :blue, linestyle = :dash);
Plots.plot!(sol_sc_point25.t, (sol_sc_point25[pbpk.C_IS_Albud_LI] + sol_sc_point25[pbpk.C_IS_Exo_LI]) * nmol_per_mol, label = "sc, 0.25 mg, large interstine", color = :red, linestyle = :dot);
Plots.plot!(sol_oral_1point5.t, (sol_oral_1point5[pbpk.C_IS_Albud_LI] + sol_oral_1point5[pbpk.C_IS_Exo_LI]) * nmol_per_mol, label = "oral, 1.5 mg, large interstine", color = :blue, linestyle = :dot);
display(plt_oral_sc)

savefig(plt_oral_sc, @projectroot("deliv/figure/pk-semaglutide-sc-oral-gi.png"));
