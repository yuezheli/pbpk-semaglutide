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
    # Time slider in units of Days
    time_slider      = Slider(1.0:1.0:100.0, value = 35.0) 
    
    # Sliders for PK Parameters
    bio_sc_slider   = Slider(0.1:0.01:1.0, value = 0.84)
    ka_sc_slider    = Slider(0.01:0.0001:0.1, value = 0.0253)
    bio_oral_slider = Slider(0.005:0.001:0.05, value = 0.01)
    ka_oral_slider  = Slider(0.01:0.0001:0.1, value = 0.0286)
    
    # KD slider replaces kon and koff
    KD_slider       = Slider(100.0:100.0:10000.0, value = 3100.0)
    
    # Dropdown for Organs 
    organ_options = ["Plasma", "small intestine (SI)", "large intestine (LI)", "Liver", "Lung", "Kidney", "Heart", "Brain"]
    organ_dropdown = Dropdown(organ_options)

    # --- 2. Observables for Plotting ---
    # Line plot observables
    t_data_sc    = Observable(Float64[])
    t_data_oral  = Observable(Float64[])
    sc_conc      = Observable(Float64[])
    oral_conc    = Observable(Float64[])
    plot_title   = Observable("Organ Concentration: Plasma")

    # --- 3. Reactive Simulation Logic ---
    onany(sc_dose_slider.value, oral_dose_slider.value, time_slider.value, 
          bio_sc_slider.value, ka_sc_slider.value, bio_oral_slider.value, ka_oral_slider.value, 
          KD_slider.value, organ_dropdown.value) do sc_dose, oral_dose, t_end_days, bio_sc, ka_sc, bio_or, ka_or, KD_nM, organ_ui_name

        # Convert user input Time (days) to Time (hours) for the ODE solver
        t_end_hr = t_end_days * hr_per_day

        # Compute koff dynamically while keeping kon fixed
        kon_fixed = 1E6
        koff_computed = KD_nM * 1E-9 * kon_fixed

        # Update human parameters with the computed koff
        p_homo = updatevolume_human(pbpk_sys, kon_Alb = kon_fixed, koff_Alb = koff_computed)
        prob_homo = ODEProblem(pbpk_sys, p_homo, (0.0, t_end_hr))

        # Setup SC Problem and Callbacks
        u0_sc = Dict([pbpk_sys.A_SC => (sc_dose * 1E-3) / MW_semaglutide])
        p_abs_sc = Dict([pbpk_sys.bioavailability => bio_sc, pbpk_sys.k_a => ka_sc])
        prob_sc = remake(prob_homo, u0=u0_sc, p=p_abs_sc)
        
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

        # Map UI Dropdown name to Symbolic MTK parameter name
        organ_symbol_map = Dict(
            "Plasma" => "Plasma",
            "small intestine (SI)" => "SI",
            "large intestine (LI)" => "LI",
            "Liver" => "Liver",
            "Lung" => "Lung",
            "Kidney" => "Kidney",
            "Heart" => "Heart",
            "Brain" => "Brain"
        )
        organ_key = organ_symbol_map[organ_ui_name]

        # --- Dynamic Organ Selection ---
        if organ_key == "Plasma"
            raw_y_sc = (sol_sc[pbpk.C_Plasma_Albud] .+ sol_sc[pbpk.C_Plasma_Exo]) .* nmol_per_mol
            raw_y_oral = (sol_oral[pbpk.C_Plasma_Albud] .+ sol_oral[pbpk.C_Plasma_Exo]) .* nmol_per_mol
        else
            var_albud = getproperty(pbpk, Symbol("C_IS_Albud_$organ_key"))
            var_exo   = getproperty(pbpk, Symbol("C_IS_Exo_$organ_key"))
            raw_y_sc = (sol_sc[var_albud] .+ sol_sc[var_exo]) .* nmol_per_mol
            raw_y_oral = (sol_oral[var_albud] .+ sol_oral[var_exo]) .* nmol_per_mol
        end

        # Clamp values to 1e-9 to prevent log10(0.0) from returning -Inf and crashing WGLMakie
        y_sc = max.(1e-9, raw_y_sc)
        y_oral = max.(1e-9, raw_y_oral)

        # Update Observables (triggers plot update)
        # Convert the hours array output from the solver back into days for plotting
        t_data_sc[] = sol_sc.t ./ hr_per_day
        t_data_oral[] = sol_oral.t ./ hr_per_day
        
        sc_conc[] = y_sc
        oral_conc[] = y_oral
        plot_title[] = "Organ Concentration: $organ_ui_name"
    end

    # Trigger the simulation once on startup to populate the plot
    notify(organ_dropdown.value)

    # --- 4. Plotting with WGLMakie ---
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1], 
              xlabel = "Time (days)",   
              ylabel = "Semaglutide concentration (nM)", 
              yscale = log10,
              xticks = 0:7:100,
              title  = plot_title)
    
    ylims!(ax, 1E-2, 1E3)

    # Plot simulation lines
    lines!(ax, t_data_sc, sc_conc, label="SC Dosing", color=:red)
    lines!(ax, t_data_oral, oral_conc, label="Oral Dosing", color=:blue)
    
    axislegend(ax, position=:rt)

    # --- 5. Front-End Layout ---
    # Using Flexbox CSS to align the label, slider, and the observable value nicely inline
    return DOM.div(
        DOM.h2("Semaglutide PBPK Interactive Simulation"),
        DOM.div(
            DOM.div(style="margin-bottom: 15px;", "Organ Target: ", organ_dropdown),
            DOM.div(style="display: flex; align-items: center; gap: 10px; margin-bottom: 5px;", "Sim Time (days): ", time_slider, time_slider.value),  
            DOM.h4("Subcutaneous Parameters"),
            DOM.div(style="display: flex; align-items: center; gap: 10px; margin-bottom: 5px;", "Dose (mg): ", sc_dose_slider, sc_dose_slider.value),
            DOM.div(style="display: flex; align-items: center; gap: 10px; margin-bottom: 5px;", "Bioavailability: ", bio_sc_slider, bio_sc_slider.value),
            DOM.div(style="display: flex; align-items: center; gap: 10px; margin-bottom: 5px;", "Absorption Rate (ka): ", ka_sc_slider, ka_sc_slider.value),
            DOM.h4("Oral Parameters"),
            DOM.div(style="display: flex; align-items: center; gap: 10px; margin-bottom: 5px;", "Dose (mg): ", oral_dose_slider, oral_dose_slider.value),
            DOM.div(style="display: flex; align-items: center; gap: 10px; margin-bottom: 5px;", "Bioavailability: ", bio_oral_slider, bio_oral_slider.value),
            DOM.div(style="display: flex; align-items: center; gap: 10px; margin-bottom: 5px;", "Absorption Rate (ka): ", ka_oral_slider, ka_oral_slider.value),
            DOM.h4("Albumin Binding Kinetics"),
            DOM.div(style="display: flex; align-items: center; gap: 10px; margin-bottom: 5px;", "Binding Affinity (K_D in nM): ", KD_slider, KD_slider.value),
            style="width: 400px; float: left; padding-right: 20px;" # Increased width to fit the numbers
        ),
        DOM.div(fig, style="float: left;")
    )
end

# Check if 'server' exists in the current session and close it safely
if @isdefined(server)
    try
        close(server)
        println("Closed existing server.")
    catch
        # Ignore if it was already closed
    end
end

# Start the local server
println("Starting new server at http://127.0.0.1:8080")
server = Bonito.Server(app, "127.0.0.1", 8080)