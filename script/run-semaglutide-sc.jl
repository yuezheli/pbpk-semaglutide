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

## observed PK, semaglutide, sc, 0.5 mg
# https://pubmed.ncbi.nlm.nih.gov/30788808/, Figure 2
plasma_sc = DataFrame(
    time_d = [0.25, 0.51, 0.66, 1.0, 1.28, 1.5, 1.78, 2.0, 2.27, 2.52, 2.99, 4.01, 5.01, 6.01, 7.0],
    conc_nM = [8.88, 10.73, 11.4, 12.1, 12.2, 11.99, 11.88, 11.88, 11.77, 11.67, 11.27, 10.87, 9.87, 9.12, 8.08]
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
affect_sc!(integrator) = integrator[pbpk_sys.A_SC] += 0.5E-3/MW_semaglutide;
cb_sc = PeriodicCallback(affect_sc!, 24.0 * 7.0); 

# oral
u0_oral = Dict([pbpk_sys.A_SC => 1.5E-3/MW_semaglutide]); # mol 
p_abs_oral = Dict([pbpk_sys.bioavailability => 0.01, pbpk_sys.k_a => 0.0286]);
prob_oral = remake(prob_homo, u0=u0_oral, p = p_abs_oral);
affect_oral!(integrator) = integrator[pbpk_sys.A_SC] += 1.5E-3/MW_semaglutide;
cb_oral = PeriodicCallback(affect_oral!, 24.0);

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
