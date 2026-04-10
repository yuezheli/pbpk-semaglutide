# date:2/23/2026 
# author: Yuezhe Li 
# purpose of this: to test run model 

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

# observed PK, semaglutide, IV, 0.5mg
# https://pubmed.ncbi.nlm.nih.gov/30788808/, Figure 2
plasma_data = DataFrame(
    time_d = [0.03, 0.33, 0.42, 0.61, 0.67, 0.85, 1.01, 1.27, 1.5, 1.76, 2.02, 2.27, 2.51, 3.02, 4.01, 5.01, 6.01, 7.0, 9.98, 13.94, 17.95, 20.94],
    conc_nM = [18.66, 13.33, 12.2, 11.16, 10.11, 9.34, 9.43, 8.55, 7.9, 7.15, 7.22, 6.67, 6.35, 6.23, 5.75, 5.06, 4.54, 4.07, 2.79, 1.77, 1.04, 1.12]
);

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

u0_exg = Dict([pbpk_sys.C_Plasma_Exo => 0.5E-3/4113/3.126]);
prob_iv = ODEProblem(pbpk_sys, merge(p_homo, u0_exg), tspan);

sol = solve(prob_iv, Rodas5(), reltol=1e-6, abstol=1e-9, saveat = 1.);

# visualization 
plt_plasma = plot(xlabel="Time (hr)", ylabel="Plasma semaglutude concentration (nM)", yaxis = :log10, yticks = [1E-1, 1, 10, 1E2, 1E3], ylims = [1E-5, 1E3]);
plot!(sol.t, (sol[pbpk.C_Plasma_Albud] + sol[pbpk.C_Plasma_Exo]) * nmol_per_mol, label = "sims");
plot!(plasma_data.time_d * hr_per_day, plasma_data.conc_nM, seriestype = :scatter, label = "Overgaard et al., 2019"); 
display(plt_plasma)

savefig(plt_plasma, @projectroot("deliv/figure/plasma-pk-semaglutide-iv.png"));
