# date: 4/19/2026 
# author: Yuezhe Li 
# purpose of this code: to build a web GUI to compare sc injection and oral 

using ProjectRoot
using Pkg; Pkg.activate(@projectroot);
using ModelingToolkit
using ModelingToolkit: getdefault
using DifferentialEquations
using SymbolicIndexingInterface
using Sundials 
using WGLMakie
using Bonito

## Setup 
include(@projectroot("script/constants.jl"))
include(@projectroot("model/twopore_albumin_peptide.jl"))
include(@projectroot("script/helper-volume-mouse2human.jl"))
@named pbpk = albumin_peptide_pbpk_mtk();
pbpk_sys = structural_simplify(pbpk);

## App
app = App() do session::Session
    # --- 1. Define UI Widgets ---
    # Sliders for Dosing and Time
    sc_dose_slider   = Slider(0.05:0.05:1.0, value = 0.25)
    oral_dose_slider = Slider(0.5:0.5:10.0, value = 1.5)
    time_slider      = Slider(100.0:100.0:1000.0, value = 800.0)
    
    # Sliders for PK Parameters
    bio_sc_slider   = Slider(0.1:0.01:1.0, value = 0.84)
    ka_sc_slider    = Slider(0.01:0.0001:0.1, value = 0.0253)
    bio_oral_slider = Slider(0.005:0.001:0.05, value = 0.01)
    ka_oral_slider  = Slider(0.01:0.0001:0.1, value = 0.0286)
    kon_slider      = Slider(1E5:1E5:5E6, value = 1E6)
    koff_slider     = Slider(0.5:0.1:10.0, value = 3.1)
    
    # Dropdown for Organs
    organ_dropdown = Dropdown(["Plasma", "SI", "LI", "Liver", "Lung", "Kidney", "Heart", "Brain"])

    # --- 2. Observables for Plotting ---
    # Line plot observables (FIXED: Separate time observables for SC and Oral)
    t_data_sc    = Observable(Float64[])
    t_data_oral  = Observable(Float64[])
    sc_conc      = Observable(Float64[])
    oral_conc    = Observable(Float64[])
    plot_title   = Observable("Organ Concentration: Plasma")
    
    # Scatter plot observables (for experimental data)
    obs_sc_scatter_t   = Observable(Float64[])
    obs_sc_scatter_y   = Observable(Float64[])
    obs_oral_scatter_t = Observable(Float64[])
    obs_oral_scatter_y = Observable(Float64[])

    # --- 3. Reactive Simulation Logic ---
    onany(sc_dose_slider.value, oral_dose_slider.value, time_slider.value, 
          bio_sc_slider.value, ka_sc_slider.value, bio_oral_slider.value, ka_oral_slider.value, 
          kon_slider.value, koff_slider.value, organ_dropdown.value) do sc_dose, oral_dose, t_end, bio_sc, ka_sc, bio_or, ka_or, kon, koff, organ

        # Update human parameters
        p_homo = updatevolume_human(pbpk_sys, kon_Alb = kon, koff_Alb = koff)
        prob_homo = ODEProblem(pbpk_sys, p_homo, (0.0, t_end))

        # Setup SC Problem and Callbacks
        u0_sc = Dict([pbpk_sys.A_SC => (sc_dose * 1E-3) / MW_semaglutide])
        p_abs_sc = Dict([pbpk_sys.bioavailability => bio_sc, pbpk_sys.k_a => ka_sc])
        prob_sc = remake(prob_homo, u0=u0_sc, p=p_abs_sc)
        
        # Use anonymous functions to prevent closure bugs
        affect_sc! = integrator -> integrator[pbpk_sys.A_SC] += (sc_dose * 1E-3) / MW_semaglutide
        cb_sc = PeriodicCallback(affect_sc!, hr_per_day * 7.0)

        # Setup Oral Problem and Callbacks
        u0_oral = Dict([pbpk_sys.A_SC => (oral_dose * 1E-3) / MW_semaglutide])
        p_abs_oral = Dict([pbpk_sys.bioavailability => bio_or, pbpk_sys.k_a => ka_or])
        prob_oral = remake(prob_homo, u0=u0_oral, p=p_abs_oral)
        
        affect_oral! = integrator -> integrator[pbpk_sys.A_SC] += (oral_dose * 1E-3) / MW_semaglutide
        cb_oral = PeriodicCallback(affect_oral!, hr_per_day)

        # Solve
        sol_sc = solve(prob_sc, Rodas5(), reltol=1e-6, abstol=1e-9, saveat=1.0, callback=cb_sc)
        sol_oral = solve(prob_oral, Rodas5(), reltol=1e-6, abstol=1e-9, saveat=1.0, callback=cb_oral)

        # --- Dynamic Organ Selection & Scatter Plot Logic ---
        if organ == "Plasma"
            raw_y_sc = (sol_sc[pbpk.C_Plasma_Albud] .+ sol_sc[pbpk.C_Plasma_Exo]) .* nmol_per_mol
            raw_y_oral = (sol_oral[pbpk.C_Plasma_Albud] .+ sol_oral[pbpk.C_Plasma_Exo]) .* nmol_per_mol
            
            # Show scatter data if defined in the global scope
            if @isdefined(plasma_sc) && @isdefined(plasma_oral)
                obs_sc_scatter_t[] = plasma_sc.time_d .* hr_per_day
                obs_sc_scatter_y[] = plasma_sc.conc_nM
                obs_oral_scatter_t[] = plasma_oral.time_d .* hr_per_day
                obs_oral_scatter_y[] = plasma_oral.conc_nM
            end
        else
            var_albud = getproperty(pbpk, Symbol("C_IS_Albud_$organ"))
            var_exo   = getproperty(pbpk, Symbol("C_IS_Exo_$organ"))
            raw_y_sc = (sol_sc[var_albud] .+ sol_sc[var_exo]) .* nmol_per_mol
            raw_y_oral = (sol_oral[var_albud] .+ sol_oral[var_exo]) .* nmol_per_mol
            
            # Hide scatter data for other organs
            obs_sc_scatter_t[] = Float64[]
            obs_sc_scatter_y[] = Float64[]
            obs_oral_scatter_t[] = Float64[]
            obs_oral_scatter_y[] = Float64[]
        end

        # Clamp values to 1e-9 to prevent log10(0.0) from returning -Inf and crashing WGLMakie
        y_sc = max.(1e-9, raw_y_sc)
        y_oral = max.(1e-9, raw_y_oral)

        # Update Observables (triggers plot update)
        t_data_sc[] = sol_sc.t
        t_data_oral[] = sol_oral.t
        sc_conc[] = y_sc
        oral_conc[] = y_oral
        plot_title[] = "Organ Concentration: $organ"
    end

    # Trigger the simulation once on startup to populate the plot
    notify(organ_dropdown.value)

    # --- 4. Plotting with WGLMakie ---
    # Using 'size' instead of 'resolution'
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1], 
              xlabel = "Time (hr)", 
              ylabel = "Semaglutide concentration (nM)", 
              yscale = log10,
              title  = plot_title)
    
    ylims!(ax, 1E-2, 1E3)

    # Plot simulation lines (FIXED: Matched independent time arrays)
    lines!(ax, t_data_sc, sc_conc, label="SC Dosing (Sim)", color=:red)
    lines!(ax, t_data_oral, oral_conc, label="Oral Dosing (Sim)", color=:blue)
    
    axislegend(ax, position=:rt)

    # --- 5. Front-End Layout ---
    return DOM.div(
        DOM.h2("Semaglutide PBPK Interactive Simulation"),
        DOM.div(
            DOM.div("Organ Target: ", organ_dropdown),
            DOM.div("Sim Time (hr): ", time_slider),
            DOM.h4("Subcutaneous Parameters"),
            DOM.div("Dose (mg): ", sc_dose_slider),
            DOM.div("Bioavailability: ", bio_sc_slider),
            DOM.div("Absorption Rate (ka): ", ka_sc_slider),
            DOM.h4("Oral Parameters"),
            DOM.div("Dose (mg): ", oral_dose_slider),
            DOM.div("Bioavailability: ", bio_oral_slider),
            DOM.div("Absorption Rate (ka): ", ka_oral_slider),
            DOM.h4("Albumin Binding Kinetics"),
            DOM.div("k_on: ", kon_slider),
            DOM.div("k_off: ", koff_slider),
            style="width: 300px; float: left; padding-right: 20px;"
        ),
        DOM.div(fig, style="float: left;")
    )
end


## Server setup
if @isdefined(server)
    try
        close(server)
        println("Closed existing server.")
    catch
        # Ignore if it was already closed
    end
end

println("Starting new server at http://127.0.0.1:8080")
server = Bonito.Server(app, "127.0.0.1", 8080)