# date: 9/24/2025
# author: Yuezhe Li 
# purpose of this script: IV dosing function 

# define infusion function 
function InfusionCallback(Dose_in_mg, mdl; 
    infusion_d = [0, 21], infusion_hr = 1, V_Plasma = 3.126, MW_Drug = MW_semaglutide, hr_per_day = hr_per_day)
    
    total_dose_M = Dose_in_mg*1E-3/V_Plasma/MW_Drug; # [M]
    infusion_rates = total_dose_M / infusion_hr
    infusion_end_times = infusion_d*hr_per_day .+ infusion_hr;

    cbs = [];
    p_index = parameter_index(mdl, mdl.infusion)
    for i in 1:lastindex(infusion_d)
        function infusion_on!(integrator)
            t = integrator.t
            p_tmp = integrator.p
            setindex!(p_tmp, infusion_rates, p_index)
            integrator.p = p_tmp
        end
        function infusion_off!(integrator)
            p_tmp = integrator.p
            setindex!(p_tmp, 0.0, p_index)
            integrator.p = p_tmp
        end
        cb_on = PresetTimeCallback(infusion_d[i]*hr_per_day, infusion_on!, save_positions = (false,false))
        cb_off = PresetTimeCallback(infusion_end_times[i], infusion_off!,save_positions = (false,false))
        push!(cbs, cb_on)
        push!(cbs, cb_off)
    end
    
    cbset = CallbackSet(cbs...);

    return cbset
end