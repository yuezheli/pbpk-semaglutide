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
using GLMakie 

## observed PK, semaglutide, IV, 0.5mg
# https://pubmed.ncbi.nlm.nih.gov/30788808/, Figure 2
plasma_data = DataFrame(
    time_d = [0.03, 0.33, 0.42, 0.61, 0.67, 0.85, 1.01, 1.27, 1.5, 1.76, 2.02, 2.27, 2.51, 3.02, 4.01, 5.01, 6.01, 7.0, 9.98, 13.94, 17.95, 20.94],
    conc_nM = [18.66, 13.33, 12.2, 11.16, 10.11, 9.34, 9.43, 8.55, 7.9, 7.15, 7.22, 6.67, 6.35, 6.23, 5.75, 5.06, 4.54, 4.07, 2.79, 1.77, 1.04, 1.12]
);

## simulation
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

u0_exg = Dict([pbpk_sys.C_Plasma_Exo => 0.5E-3/MW_semaglutide/3.126]);
prob_iv = ODEProblem(pbpk_sys, merge(p_homo, u0_exg), tspan);

sol = solve(prob_iv, Rodas5(), reltol=1e-6, abstol=1e-9, saveat = 1.);

## visualization 
plt_plasma = plot(xlabel="Time (hr)", ylabel="Plasma semaglutude concentration (nM)", yaxis = :log10, yticks = [1E-1, 1, 10, 1E2, 1E3], ylims = [1E-5, 1E3]);
plot!(sol.t, (sol[pbpk.C_Plasma_Albud] + sol[pbpk.C_Plasma_Exo]) * nmol_per_mol, label = "sims");
plot!(plasma_data.time_d * hr_per_day, plasma_data.conc_nM, seriestype = :scatter, label = "Overgaard et al., 2019"); 
display(plt_plasma)

savefig(plt_plasma, @projectroot("deliv/figure/plasma-pk-semaglutide-iv.png"));


## GUI

function run_pbpk_iv_app()
    # Create the main window
    fig = Figure(size = (1200, 800), fontsize = 18)
    
    # --- UI Layout ---
    # Create sliders for Dose, Kon, and Koff
    sg = SliderGrid(fig[1, 1:2],
        (label = "Dose (mg)", range = 0.1:0.1:10.0, startvalue = 0.5, format = "{:.1f} mg"),
        (label = "Albumin:Semaglutide Koff (1/hr)", range = 0.1:0.1:10.0, startvalue = 3.1, format = "{:.1f}")
    )
    
    dose_slider = sg.sliders[1].value
    koff_slider = sg.sliders[2].value

    # --- Data Observables ---
    # These act as dynamic variables. When updated, the plot updates automatically.
    time_obs   = Observable(Float64[])
    plasma_obs = Observable(Float64[])
    si_obs     = Observable(Float64[])
    li_obs     = Observable(Float64[])

    # --- Plot Setup ---
    ax = Axis(fig[2, 1:2], 
        title = "Semaglutide PBPK Simulation",
        xlabel = "Time (hr)", 
        ylabel = "Concentration (nM)", 
        yscale = log10
    )
    
    # Bind observables to plot lines
    lines!(ax, time_obs, plasma_obs, label = "Plasma (Total)", linewidth = 3, color = :blue)
    lines!(ax, time_obs, si_obs, label = "Small Intestine (IS)", linewidth = 3, color = :orange)
    lines!(ax, time_obs, li_obs, label = "Large Intestine (IS)", linewidth = 3, color = :green)
    
    axislegend(ax, position = :rt)
    GLMakie.ylims!(ax, 1e-2, 1e3)

    # --- The Update Engine ---
    # This block runs every time any slider is moved
    onany(dose_slider, koff_slider) do dose, koff
        
        new_prob = deepcopy(prob_iv)

        new_prob[pbpk.C_Plasma_Exo] = (dose * 1E-3) / MW_semaglutide / 3.126
        
        new_prob.ps[pbpk_sys.koff_Alb] = koff
        
        sol = solve(new_prob, Rodas5(), reltol=1e-6, abstol=1e-9, saveat=1.0) 
        
        # 4. Extract data and push to Observables
        time_obs[] = sol.t
        
        # Calculate Total Plasma (Free + Bound) [cite: 19]
        plasma_obs[] = (sol[pbpk.C_Plasma_Albud] .+ sol[pbpk.C_Plasma_Exo]) .* 1e9
        
        # Calculate Total Interstitial Small Intestine (Free + Bound)
        si_obs[] = (sol[pbpk.C_IS_Albud_SI] .+ sol[pbpk.C_IS_Exo_SI]) .* 1e9
        
        # Calculate Total Interstitial Large Intestine (Free + Bound)
        li_obs[] = (sol[pbpk.C_IS_Albud_LI] .+ sol[pbpk.C_IS_Exo_LI]) .* 1e9
    end

    # Trigger the initial solve to populate the graph on startup
    notify(dose_slider)

    # Launch the App
    display(fig)
end

# Run the app
run_pbpk_iv_app()

