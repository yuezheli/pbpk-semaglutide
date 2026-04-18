# date: 2/23/2026 
# author: Yuezhe Li 
# purpose of this code: to code up albumin-binding exogenous PK model (in preparation of incorporating albumin:PDC binding in the future)
# default parameters were for mouse 
# https://pubmed.ncbi.nlm.nih.gov/38691205/

function albumin_peptide_pbpk_mtk(; name)
    @independent_variables t # [unit = u"hr"]
    D = Differential(t)

    pars = @parameters begin      
        # --- Drug Properties (Mouse Albumin) ---
        MW = 67.0E3              # [unit = u"Dalton"]
        MW_albud = 4113          # [unit = u"Dalton"]
        BW = 0.028               # [unit = u"kg"]
        
        # --- Albumin-Specific Disposition Factors ---
        Spino_alb = 0.579        # [unit = u""] Pinocytosis modification factor for albumin
        Frea_alb = 0.97          # [unit = u""] Fraction of filtered albumin reabsorbed (97%)
        Ksyn_alb = 2.8E-4        # [unit = u"M/hr"] Endogenous synthesis rate
        
        # --- Endosomal Processing (FcRn) ---
        k_deg = 15.3             # [unit = u"1/hr"] Lysosomal degradation rate
        CL_up_baseline = 1.22    # [unit = u"L/hr/L"] Pinocytosis rate per endosomal volume
        FR = 0.715               # [unit = u""] Recycling fraction to vascular space
        kon_FcRn = 5.8E7         # [unit = u"1/M/hr"] FcRn association rate
        koff_FcRn = 43.9         # [unit = u"1/hr"] FcRn dissociation rate
        FcRn_total = 1.2534E-4   # [unit = u"M"] Total endosomal FcRn concentration
        C_0_endo_Alb = 5.2E-4    # [unit = u"M"] endogenous albumin steady state concentration
        Epsilon = 1E-16          # a small volume to avoid numeric issue
        
        # --- Two-Pore & Physiological Constants ---
        GFR = 0.0167             # [unit = u"L/hr"] Glomerular filtration rate
        sigma_IS = 0.2           # [unit = u""] Lymphatic reflection coefficient
        x_j = 0.38               # [unit = u""] Isogravimetric flow constant
        X_p = 13197.0            # [unit = u"nm^3"] Two-pore proportionality constant
        alpha_L = 0.042          # [unit = u""] Large pore hydraulic conductance
        alpha_S = 0.958          # [unit = u""] Small pore hydraulic conductance
        r_L = 22.85              # [unit = u"nm"] Large pore radius
        r_S = 4.44               # [unit = u"nm"] Small pore radius

        # --- Albumin:exogenous binding ---
        kon_Alb = 3.6E9        # [unit = u"1/M/hr"] 
        koff_Alb = 3.1         # [unit = u"1/hr"] 

        # --- Drug absorption ---
        bioavailability = 0.8       # drug bioavailability
        k_a = 0.1                   # drug absorption rate [unit = u"1/hr"] 

        # --- Drug infusion ---
        infusion = 0.       # drug infusion rate 
        
        # --- Volumes [L] ---
        V_Plasma = 0.000944
        V_LN = 0.000113
        V_V_Heart = 0.00000585
        V_V_Lung = 0.0000295
        V_V_Liver = 0.000164
        V_V_Muscle = 0.000249
        V_V_Skin = 0.000188
        V_V_Adipose = 0.0000218
        V_V_Bone = 0.0000621
        V_V_Brain = 0.0000107
        V_V_Kidney = 0.0000289
        V_V_SI = 0.0000116
        V_V_LI = 0.000005
        V_V_Pancreas = 0.00000534
        V_V_Thymus = 0.0000005
        V_V_Spleen = 0.0000154
        V_V_Other = 0.0000195
        
        V_IS_Heart = 0.0000217
        V_IS_Lung = 0.0000384
        V_IS_Muscle = 0.00147
        V_IS_Skin = 0.00166
        V_IS_Adipose = 0.000337
        V_IS_Bone = 0.000525
        V_IS_Brain = 0.0000873
        V_IS_Kidney = 0.0000788
        V_IS_Liver = 0.000385
        V_IS_SI = 0.000127
        V_IS_LI = 0.0000545
        V_IS_Pancreas = 0.0000169
        V_IS_Thymus = 0.00000153
        V_IS_Spleen = 0.0000254
        V_IS_Other = 0.0000797
        
        V_E_Heart = 7.6E-7
        V_E_Lung = 1.02E-6
        V_E_Muscle = 5.66E-5
        V_E_Skin = 2.51E-5
        V_E_Adipose = 9.91E-6
        V_E_Bone = 1.41E-5
        V_E_Brain = 2.43E-6
        V_E_Kidney = 2.63E-6
        V_E_Liver = 9.63E-6
        V_E_SI = 3.64E-6
        V_E_LI = 1.57E-6
        V_E_Pancreas = 4.85E-7
        V_E_Thymus = 5.0E-8
        V_E_Spleen = 6.35E-7
        V_E_Other = 2.33E-6
        
        # --- Plasma Flows [L/hr] ---
        Q_Heart = 0.0365
        Q_Muscle = 0.0861
        Q_Skin = 0.0278
        Q_Adipose = 0.0134
        Q_Bone = 0.0152
        Q_Brain = 0.0118
        Q_Kidney = 0.0685
        Q_Liver = 0.0103
        Q_SI = 0.0581
        Q_LI = 0.0173
        Q_Pancreas = 0.00624
        Q_Thymus = 0.00119
        Q_Spleen = 0.00818
        Q_Other = 0.0109
        # Q_Lung = 0.35171
    end

    # Total Venous Return (Lung Arterial Flow)
    global Q_Lung = Q_Heart + Q_Muscle + Q_Skin + Q_Adipose + Q_Bone + Q_Brain + Q_Kidney + Q_Liver + Q_SI + Q_LI + Q_Pancreas + Q_Thymus + Q_Spleen + Q_Other
    global L_Lymph_drain = (Q_Lung + Q_Heart + Q_Muscle + Q_Skin + Q_Adipose + Q_Bone + Q_Brain + Q_Kidney + Q_Thymus + Q_Other + Q_Liver)*0.002;  # [L/hr]

    # enumerate organ
    Lung = 1
    Liver = 2
    Heart = 3
    Muscle = 4
    Skin = 5
    Adipose = 6
    Bone = 7
    Brain = 8
    Kidney = 9
    SI = 10
    LI = 11
    Pancreas = 12
    Thymus = 13
    Spleen = 14
    Other = 15

    # Variables
    # 2. Variable Enumeration Strategy 
    # (Fixed the LoadError by defining variables outside the macro block)
    orgs = [:Lung, :Heart, :Kidney, :Muscle, :Skin, :Liver, :Brain, :Adipose, 
            :Thymus, :Bone, :SI, :LI, :Spleen, :Pancreas, :Other]

    Qs = [Q_Lung, Q_Heart, Q_Kidney, Q_Muscle, Q_Skin, Q_Liver, Q_Brain, Q_Adipose, Q_Thymus, Q_Bone, Q_SI, Q_LI, Q_Spleen, Q_Pancreas, Q_Other]
    VVs = [V_V_Lung, V_V_Heart, V_V_Kidney, V_V_Muscle, V_V_Skin, V_V_Liver, V_V_Brain, V_V_Adipose, V_V_Thymus, V_V_Bone, V_V_SI, V_V_LI, V_V_Spleen, V_V_Pancreas, V_V_Other]
    VISs = [V_IS_Lung, V_IS_Heart, V_IS_Kidney, V_IS_Muscle, V_IS_Skin, V_IS_Liver, V_IS_Brain, V_IS_Adipose, V_IS_Thymus, V_IS_Bone, V_IS_SI, V_IS_LI, V_IS_Spleen, V_IS_Pancreas, V_IS_Other]
    VEs = [V_E_Lung, V_E_Heart, V_E_Kidney, V_E_Muscle, V_E_Skin, V_E_Liver, V_E_Brain, V_E_Adipose, V_E_Thymus, V_E_Bone, V_E_SI, V_E_LI, V_E_Spleen, V_E_Pancreas, V_E_Other]
    
    @variables A_SC_Albud(t) = 0.0                             # Amount of drug in the depot [mol]
    @variables C_Plasma(t)=C_0_endo_Alb C_LN(t)=Epsilon        # Endogenous albumin
    @variables C_Plasma_Exo(t)=0.0 C_LN_Exo(t)=0.0             # Exogenous albumin
    @variables C_Plasma_Albud(t)=Epsilon C_LN_Albud(t)=Epsilon # Exogenous albumin-binding protein
    
    # endogenous albumin
    CV_vars  = [(@variables $(Symbol("C_V_", o))(t) = Epsilon)[1] for o in orgs] 
    CIS_vars = [(@variables $(Symbol("C_IS_", o))(t) = Epsilon)[1] for o in orgs] 
    UB_vars  = [(@variables $(Symbol("C_E_UB_", o))(t) = Epsilon)[1] for o in orgs] 
    B_vars   = [(@variables $(Symbol("C_E_B_", o))(t) = Epsilon)[1] for o in orgs] 
    # exogenous albumin:peptide 
    CV_Exo_vars  = [(@variables $(Symbol("C_V_Exo_", o))(t) = 0.0)[1] for o in orgs] 
    CIS_Exo_vars = [(@variables $(Symbol("C_IS_Exo_", o))(t) = 0.0)[1] for o in orgs] 
    UB_Exo_vars  = [(@variables $(Symbol("C_E_UB_Exo_", o))(t) = 0.0)[1] for o in orgs] 
    B_Exo_vars   = [(@variables $(Symbol("C_E_B_Exo_", o))(t) = 0.0)[1] for o in orgs]
    # exogenous albumin-binding proteins 
    CV_Albud_vars  = [(@variables $(Symbol("C_V_Albud_", o))(t) = Epsilon)[1] for o in orgs] 
    CIS_Albud_vars = [(@variables $(Symbol("C_IS_Albud_", o))(t) = Epsilon)[1] for o in orgs] 
    CE_Albud_vars = [(@variables $(Symbol("C_E_Albud_", o))(t) = Epsilon)[1] for o in orgs] 
    # FcRn
    F_vars   = [(@variables $(Symbol("FcRn_", o))(t) = FcRn_total)[1] for o in orgs]
    vars = [C_Plasma; C_Plasma_Exo; A_SC_Albud; C_LN; C_LN_Exo; CV_vars; CV_Exo_vars; CIS_vars; CIS_Exo_vars; UB_vars; UB_Exo_vars; B_vars; B_Exo_vars; F_vars; C_Plasma_Albud; C_LN_Albud; CV_Albud_vars; CE_Albud_vars; CIS_Albud_vars]

    # Derived Molecular Constants (albumin)
    a_e = 0.5614 * (MW/1000)^(1/3) + 0.09611 * (MW/1000)^(2/3)
    theta = (1 / (1 + (a_e/2.95)^7.11))^3.8  # Sepp et al. (Eq 21)
    sigma_L = 3.5E-5 * MW^0.717
    sigma_S = 1 - 0.8489 * exp(-4.0E-5 * MW)

    # Derived Molecular Constants (exogenous albumin-binding protein)
    a_e_albud = 0.5614 * (MW_albud/1000)^(1/3) + 0.09611 * (MW_albud/1000)^(2/3)
    theta_albud = (1 / (1 + (a_e_albud/2.95)^7.11))^3.8  # Sepp et al. (Eq 21)
    sigma_L_albud = 3.5E-5 * MW_albud^0.717
    sigma_S_albud = 1 - 0.8489 * exp(-4.0E-5 * MW_albud)

    eqs = Equation[]
    
    for (i, org) in enumerate(orgs)
        Q_vals = [Q_Lung, Q_Heart, Q_Kidney, Q_Muscle, Q_Skin, Q_Liver, Q_Brain, Q_Adipose, 
                  Q_Thymus, Q_Bone, Q_SI, Q_LI, Q_Spleen, Q_Pancreas, Q_Other]
        Q_org = Q_vals[i]

        J_org = Q_org * 0.002
        J_L = x_j * J_org + alpha_L * J_org
        J_S = -x_j * J_org + alpha_S * J_org
        
        # Two-pore permeability products (PS) and clearances (CLtp) calculation (albumin)
        PS_L = X_p * (1/a_e) * (0.3429*exp(-1.2175E-4*MW + 0.6571*exp(-4.21E-6*MW))) * alpha_L / (r_L^2) * J_org
        PS_S = X_p * (1/a_e) * (0.2353*exp(-8.295E-5*MW + 0.7648*exp(-5.3095E-4*MW))) * alpha_S / (r_S^2) * J_org

        # Two-pore permeability products (PS) and clearances (CLtp) calculation (exogenous albumin-binding protein)
        PS_L_albud = X_p * (1/a_e_albud) * (0.3429*exp(-1.2175E-4*MW_albud + 0.6571*exp(-4.21E-6*MW_albud))) * alpha_L / (r_L^2) * J_org
        PS_S_albud = X_p * (1/a_e_albud) * (0.2353*exp(-8.295E-5*MW_albud + 0.7648*exp(-5.3095E-4*MW_albud))) * alpha_S / (r_S^2) * J_org
        
        # State, endogenous albumin 
        CV   = CV_vars[i]   
        CIS  = CIS_vars[i]   
        CEUB = UB_vars[i]    
        CEB  = B_vars[i]     
        # State, exogenous albumin 
        CV_X = CV_Exo_vars[i]
        CIS_X = CIS_Exo_vars[i]
        CEUB_X = UB_Exo_vars[i]
        CEB_X = B_Exo_vars[i]
        # FcRn
        FRN  = F_vars[i]   
        # State, exogenous albumin-binding peptide 
        CV_Albud = CV_Albud_vars[i]
        CE_Albud = CE_Albud_vars[i]
        CIS_Albud = CIS_Albud_vars[i]

        # volume
        VV = VVs[i]
        VIS = VISs[i]
        VE = VEs[i]

        CLup = CL_up_baseline * VE

        Pe_L = J_L*(1-sigma_L) / (PS_L + Epsilon)
        Pe_S = J_S*(1-sigma_S) / (PS_S + Epsilon)
        CL_tp_L = PS_L * ( 1 - CIS/CV ) * Pe_L/(exp(Pe_L)-1) .+ J_L*(1-sigma_L);
        CL_tp_S = PS_S * ( 1 - CIS/CV ) * Pe_S/(exp(Pe_S)-1) .+ J_S*(1-sigma_S);
        # incorporate conditions when drug conc = 0
        CL_tp_L = ifelse.(isnan.(CL_tp_L), 0.0, CL_tp_L)
        CL_tp_S = ifelse.(isnan.(CL_tp_S), 0.0, CL_tp_S)
        CLtp = CL_tp_L + CL_tp_S

        Pe_L_albud = J_L*(1-sigma_L_albud) / (PS_L_albud + Epsilon)
        Pe_S_albud = J_S*(1-sigma_S_albud) / (PS_S_albud + Epsilon)
        CLtp_albud = (
            PS_L_albud * ( 1 - CIS_Albud/CV_Albud ) * Pe_L_albud/(exp(Pe_L_albud)-1 + Epsilon) .+ J_L*(1-sigma_L_albud) + 
            PS_S_albud * ( 1 - CIS_Albud/CV_Albud ) * Pe_S_albud/(exp(Pe_S_albud)-1 + Epsilon) .+ J_S*(1-sigma_S_albud)
        )

        # --- Vascular & Interstitial Equations ---
        if org == :Lung
            push!(eqs, D(CV) ~ ((Q_org + J_org)*C_Plasma - Q_org*CV - (CLtp + Spino_alb*CLup)*CV + 2*CLup*FR*CEB)/VV - kon_Alb*CV_Albud*CV + koff_Alb*CV_X )
            push!(eqs, D(CV_X) ~ ((Q_org + J_org)*C_Plasma_Exo - Q_org*CV_X - (CLtp + Spino_alb*CLup)*CV_X + 2*CLup*FR*CEB_X)/VV + kon_Alb*CV_Albud*CV - koff_Alb*CV_X )
            push!(eqs, D(CV_Albud) ~ ((Q_org + J_org)*C_Plasma_Albud - Q_org*CV_Albud - (CLtp_albud + CLup)*CV_Albud)/VV - kon_Alb*CV_Albud*CV + koff_Alb*CV_X )
        elseif org == :Kidney
            push!(eqs, D(CV) ~ (Q_org*CV_vars[Lung] - (Q_org - J_org)*CV - (CLtp + Spino_alb*CLup)*CV + 2*CLup*FR*CEB - GFR*theta*(1-Frea_alb)*CV)/VV - kon_Alb*CV_Albud*CV + koff_Alb*CV_X)
            push!(eqs, D(CV_X) ~ (Q_org*CV_Exo_vars[Lung] - (Q_org - J_org)*CV_X - (CLtp + Spino_alb*CLup)*CV_X + 2*CLup*FR*CEB_X - GFR*theta*(1-Frea_alb)*CV_X)/VV + kon_Alb*CV_Albud*CV - koff_Alb*CV_X)
            push!(eqs, D(CV_Albud) ~ (Q_org*CV_Albud_vars[Lung] - (Q_org - J_org)*CV_Albud - (CLtp_albud + CLup)*CV_Albud - GFR*theta_albud*CV_Albud)/VV - kon_Alb*CV_Albud*CV + koff_Alb*CV_X)
        elseif org == :Liver
            # Splanchnic return into liver
            influx_endo = (
                Q_org*CV_vars[Lung]
                + (Q_vals[Pancreas]-Q_vals[Pancreas]*0.002)*CV_vars[Pancreas]
                + (Q_vals[SI]-Q_vals[SI]*0.002)*CV_vars[SI]
                + (Q_vals[LI]-Q_vals[LI]*0.002)*CV_vars[LI]
                + (Q_vals[Spleen]-Q_vals[Spleen]*0.002)*CV_vars[Spleen] 
            )
            influx_exo = (
                Q_org*CV_Exo_vars[Lung]
                + (Q_vals[Pancreas]-Q_vals[Pancreas]*0.002)*CV_Exo_vars[Pancreas]
                + (Q_vals[SI]-Q_vals[SI]*0.002)*CV_Exo_vars[SI]
                + (Q_vals[LI]-Q_vals[LI]*0.002)*CV_Exo_vars[LI]
                + (Q_vals[Spleen]-Q_vals[Spleen]*0.002)*CV_Exo_vars[Spleen] 
            )
            influx_albud = (
                Q_org*CV_Albud_vars[Lung]
                + (Q_vals[Pancreas]-Q_vals[Pancreas]*0.002)*CV_Albud_vars[Pancreas]
                + (Q_vals[SI]-Q_vals[SI]*0.002)*CV_Albud_vars[SI]
                + (Q_vals[LI]-Q_vals[LI]*0.002)*CV_Albud_vars[LI]
                + (Q_vals[Spleen]-Q_vals[Spleen]*0.002)*CV_Albud_vars[Spleen] 
            )
            efflux_endo = (
                (Q_org-J_org)
                + (Q_vals[Pancreas]-Q_vals[Pancreas]*0.002)
                + (Q_vals[Spleen]-Q_vals[Spleen]*0.002)
                + (Q_vals[SI]-Q_vals[SI]*0.002)
                + (Q_vals[LI]-Q_vals[LI]*0.002)
            )*CV
            efflux_exo = (
                (Q_org-J_org)
                + (Q_vals[Pancreas]-Q_vals[Pancreas]*0.002)
                + (Q_vals[Spleen]-Q_vals[Spleen]*0.002)
                + (Q_vals[SI]-Q_vals[SI]*0.002)
                + (Q_vals[LI]-Q_vals[LI]*0.002)
            )*CV_X
            efflux_albud = (
                (Q_org-J_org)
                + (Q_vals[Pancreas]-Q_vals[Pancreas]*0.002)
                + (Q_vals[Spleen]-Q_vals[Spleen]*0.002)
                + (Q_vals[SI]-Q_vals[SI]*0.002)
                + (Q_vals[LI]-Q_vals[LI]*0.002)
            )*CV_Albud
            push!(eqs, D(CV) ~ ( influx_endo - efflux_endo - (CLtp + Spino_alb*CLup)*CV + 2*CLup*FR*CEB)/VV - kon_Alb*CV_Albud*CV + koff_Alb*CV_X )
            push!(eqs, D(CV_X) ~ ( influx_exo - efflux_exo - (CLtp + Spino_alb*CLup)*CV_X + 2*CLup*FR*CEB_X)/VV + kon_Alb*CV_Albud*CV - koff_Alb*CV_X )
            push!(eqs, D(CV_Albud) ~ ( influx_albud - efflux_albud - (CLtp + CLup)*CV_Albud)/VV - kon_Alb*CV_Albud*CV + koff_Alb*CV_X )
        else
            # normal organs
            push!(eqs, D(CV) ~ (Q_org*CV_vars[Lung] - (Q_org - J_org)*CV - (CLtp + Spino_alb*CLup)*CV + 2*CLup*FR*CEB)/VV - kon_Alb*CV_Albud*CV + koff_Alb*CV_X )
            push!(eqs, D(CV_X) ~ (Q_org*CV_Exo_vars[Lung] - (Q_org - J_org)*CV_X - (CLtp + Spino_alb*CLup)*CV_X + 2*CLup*FR*CEB_X)/VV + kon_Alb*CV_Albud*CV - koff_Alb*CV_X )
            push!(eqs, D(CV_Albud) ~ (Q_org*CV_Albud_vars[Lung] - (Q_org - J_org)*CV_Albud - (CLtp + CLup)*CV_Albud)/VV - kon_Alb*CV_Albud*CV + koff_Alb*CV_X )
        end
        
        push!(eqs, D(CIS) ~ (CLtp*CV - (1-sigma_IS)*J_org*CIS - Spino_alb*CLup*CIS + 2*CLup*(1-FR)*CEB)/VIS - kon_Alb*CIS_Albud*CIS + koff_Alb*CIS_X )
        push!(eqs, D(CIS_X) ~ (CLtp*CV_X - (1-sigma_IS)*J_org*CIS_X - Spino_alb*CLup*CIS_X + 2*CLup*(1-FR)*CEB_X)/VIS + kon_Alb*CIS_Albud*CIS - koff_Alb*CIS_X )
        push!(eqs, D(CIS_Albud) ~ (CLtp*CV_Albud - (1-sigma_IS)*J_org*CIS_Albud - CLup*CIS_Albud)/VIS - kon_Alb*CIS_Albud*CIS + koff_Alb*CIS_X )

        # --- Endosomal Space Equations ---
        push!(eqs, D(CEUB) ~ Spino_alb*CL_up_baseline*(CV + CIS) - kon_FcRn*CEUB*FRN + koff_FcRn*CEB - k_deg*CEUB)
        push!(eqs, D(CEUB_X) ~ Spino_alb*CL_up_baseline*(CV_X + CIS_X) - kon_FcRn*CEUB_X*FRN + koff_FcRn*CEB_X - k_deg*CEUB_X)
        push!(eqs, D(CEB) ~ kon_FcRn*CEUB*FRN - koff_FcRn*CEB - 2*CL_up_baseline*CEB)
        push!(eqs, D(CEB_X) ~ kon_FcRn*CEUB_X*FRN - koff_FcRn*CEB_X - 2*CL_up_baseline*CEB_X)
        push!(eqs, D(FRN) ~ 2*CL_up_baseline*CEB + 2*CL_up_baseline*CEB_X + koff_FcRn*CEB + koff_FcRn*CEB_X - kon_FcRn*CEUB*FRN - kon_FcRn*CEUB_X*FRN)

        push!(eqs, D(CE_Albud) ~ CL_up_baseline*(CV_Albud + CIS_Albud) - k_deg*CE_Albud)
    end

    # --- Depot of drug ---
    push!(eqs, D(A_SC_Albud) ~ (-k_a * A_SC_Albud))

    # --- Central Plasma and Lymph Node Equations ---
    venous_return_endo = sum(1:length(orgs)) do j
        if orgs[j] == :Lung
            return 0.0 # Skip lung in the summation as per your logic 
        else
            Q_vals = [Q_Lung, Q_Heart, Q_Kidney, Q_Muscle, Q_Skin, Q_Liver, Q_Brain, Q_Adipose, 
                      Q_Thymus, Q_Bone, Q_SI, Q_LI, Q_Spleen, Q_Pancreas, Q_Other]
            return (Q_vals[j] - Q_vals[j]*0.002) * CV_vars[j]
        end
    end
    push!(eqs, D(C_Plasma) ~ (venous_return_endo + L_Lymph_drain*C_LN - (Qs[Lung] + Qs[Lung]*0.002)*C_Plasma + Ksyn_alb*V_Plasma)/V_Plasma
                                - kon_Alb*C_Plasma*C_Plasma_Albud + koff_Alb*C_Plasma_Exo
                                )

    venous_return_exo = sum(1:length(orgs)) do j
        if orgs[j] == :Lung
            return 0.0 # Skip lung in the summation as per your logic 
        else
            Q_vals = [Q_Lung, Q_Heart, Q_Kidney, Q_Muscle, Q_Skin, Q_Liver, Q_Brain, Q_Adipose, 
                      Q_Thymus, Q_Bone, Q_SI, Q_LI, Q_Spleen, Q_Pancreas, Q_Other]
            return (Q_vals[j] - Q_vals[j]*0.002) * CV_Exo_vars[j]
        end
    end
    push!(eqs, D(C_Plasma_Exo) ~ (venous_return_exo + L_Lymph_drain*C_LN_Exo - (Qs[Lung] + Qs[Lung]*0.002)*C_Plasma_Exo)/V_Plasma
                                 + kon_Alb*C_Plasma*C_Plasma_Albud - koff_Alb*C_Plasma_Exo
                                 + k_a * A_SC_Albud * bioavailability / V_Plasma 
                                 + infusion)

    venous_return_albud = sum(1:length(orgs)) do j
        if orgs[j] == :Lung
            return 0.0 # Skip lung in the summation as per your logic 
        else
            Q_vals = [Q_Lung, Q_Heart, Q_Kidney, Q_Muscle, Q_Skin, Q_Liver, Q_Brain, Q_Adipose, 
                      Q_Thymus, Q_Bone, Q_SI, Q_LI, Q_Spleen, Q_Pancreas, Q_Other]
            return (Q_vals[j] - Q_vals[j]*0.002) * CV_Albud_vars[j]
        end
    end
    push!(eqs, D(C_Plasma_Albud) ~ (venous_return_albud + L_Lymph_drain*C_LN_Albud - (Qs[Lung] + Qs[Lung]*0.002)*C_Plasma_Albud)/V_Plasma
                                     - kon_Alb*C_Plasma*C_Plasma_Albud + koff_Alb*C_Plasma_Exo )

    Q_map = Dict(:Lung=>Q_Lung, :Heart=>Q_Heart, :Kidney=>Q_Kidney, :Muscle=>Q_Muscle, :Skin=>Q_Skin, :Liver=>Q_Liver, :Brain=>Q_Brain, :Adipose=>Q_Adipose, :Thymus=>Q_Thymus, :Bone=>Q_Bone, :SI=>Q_SI, :LI=>Q_LI, :Spleen=>Q_Spleen, :Pancreas=>Q_Pancreas, :Other=>Q_Other)

    CIS_map = Dict(zip(orgs, CIS_vars))
    l_drain = sum([(1-sigma_IS)*(Q_map[o]*0.002)*CIS_map[o] for o in orgs])
    push!(eqs, D(C_LN) ~ (l_drain - L_Lymph_drain*C_LN)/V_LN - kon_Alb*C_LN*C_LN_Albud + koff_Alb*C_LN_Exo)

    CIS_map_exo = Dict(zip(orgs, CIS_Exo_vars))
    l_drain_exo = sum([(1-sigma_IS)*(Q_map[o]*0.002)*CIS_map_exo[o] for o in orgs])
    push!(eqs, D(C_LN_Exo) ~ (l_drain_exo - L_Lymph_drain*C_LN_Exo)/V_LN + kon_Alb*C_LN*C_LN_Albud - koff_Alb*C_LN_Exo)

    CIS_map_albud = Dict(zip(orgs, CIS_Albud_vars))
    l_drain_albud = sum([(1-sigma_IS)*(Q_map[o]*0.002)*CIS_map_albud[o] for o in orgs])
    push!(eqs, D(C_LN_Albud) ~ (l_drain_albud - L_Lymph_drain*C_LN_Albud)/V_LN - kon_Alb*C_LN*C_LN_Albud + koff_Alb*C_LN_Exo)

    return ODESystem(eqs, t, vars, pars; name = name)
end