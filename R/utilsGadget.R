# Model Playground Utils

#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`


#' Model outputs for bottom right box
#'
#' @export
pars_dropdown_names <- c(
  pv_CO2 = 45, #46,      #  [mmHg] venous CO2 partial pressure [h]
  pv_O2  = 35, # 40,     #  [mmHg] venous O2 partial pressure [h]
  pi_O2  = 159,          #  [mmHg] air O2 partial pressure at sea level
  pi_CO2 = 0.3,          #  [mmHg] air CO2 partial pressure at sea level
  Vtid = 600,         # [ml]
  resp2CC_ratio = 0.2 # ratio of the respiratory cycle to cardiac cycle duration: slow breathing: 0.08; //hyperventilation: 0.8; 
)

#' Model outputs
#'
#' @export
extraOutputs <- c(
  "Qtv", "Qpv", "Qmv", "Qav", "Qaa", "Qp", "Qs",
  "Qca", "Qc", "Qcv",
  "Ppa", "Ppv", "Paorta", "Psa", "Psv", "Pca", "Pcv",
  "Plv", "Prv",
  "Vtot",
  "A_aorta_0", "A_aorta", "Aortic_strain",
  "delta_strain", "f_Baro", "HR", "P_ra", 
  # "alpha_map", "alpha_rap", "rsna", "beta_rsna",
  # "R_aa", "R_ea", "R_r", "Phi_rb",
  # "P_gh", "P_f", "Phi_gfilt", "Sigma_tgf", "Phi_filsod",
  # "gamma_filsod", "gamma_at", "gamma_rsna", "eta_ptsodreab",
  # "Phi_ptsodreab", "Phi_mdsod", "Psi_al",
  # "eta_dtsodreab", "Phi_dtsodreab", "Phi_dtsod",
  # "lambda_dt", "lambda_anp", "eta_cdsodreab",
  # "Phi_cdsodreab", "Phi_usod", "mu_al", "mu_adh",
  # "Phi_twreab", "Phi_u", "xi_ksod", "xi_map",
  # "xi_at", "N_als", "V_ecf", "C_al", "C_sod",
  # "Phi_sodin", "N_adhs", "Phi_win", "C_adh",
  "Pb_0","Pa_0", "Rc_pB_CO2", "Rs_pa_CO2", "Rs", "Rbrain", "Raorta", "Rtot",
  "Csa", "Csv", "ca_CO2", "ca_O2", "cv_CO2", "cv_O2", "cS_CO2", "cS_O2",
  "cB_CO2", "cB_O2", "Va", "d_Va", "Rav", "Rmv", "Rpv", "Rtv", "Raa"
)


#' @export
plot_options_params <- c(
  "time",               
  "VolAorta",            
  "VolSystemicArtery",                
  "VolSystemicVein",               
  "VolLeftVentricle", 
  "Vtot",
  "CardiacOutput",                 
  "SympNervActivity",        
  "MeanArterialPressure", 
  "Qaa", 
  "Qp", 
  "Qs",
  "Qca", 
  "Qc",  
  "Qcv",
  "Ppa", 
  "Ppv", 
  "Paorta", 
  "Psa", 
  "Psv",
  "Pca", 
  "Pcv",
  "Plv", 
  "Prv", 
  "HR"
)


# initial values for the gadget

#'@export
state_gadget <- c(
  #initial values for cardiovascular model
  HeartCycle = 0,                   # starting heart cycle fraction
  VolPulmonaryArtery =  112.23655,  # [ml] pulmonary arterial blood volume
  VolPulmonaryVein =  570.7688,     # [ml] pulmonary venous blood volume
  VolAorta = 138.2129,              # [ml] aortic blood volume
  VolSystemicArtery =  1125.199,    # [ml] systemic arterial blood volume
  VolSystemicVein = 3033.680,       # [ml] systemic venous blood volume
  VolCerebralArtery = 102.2822,     # 0.0237 * 4600,  # [mL] cerebral arterial blood volume
  VolCerebralVein = 386.9102,       # [mL] cerebral venous blood volume
  VolLeftVentricle =  121.000,      # [mL] Left ventricular volume [59]
  VolRightVentricle =  122.500,     # [mL] Right ventricular volume [8, 29]
  Autoregulation = 0.83,            # autoregulation
  CardiacOutput  = 5000,            # [mL/min]Cardiac output
  MeanAorticStrain = 1.3,           # Mean aortic strain
  Baro_on = 0.91,                   # Active baroreceptors
  SympNervActivity = 0.5,           # Sympathetic nervous system activity
  Renin_activity    = 0.18,
  AT2_activity      = 0.18,
  N_al              = 1,
  N_adh             = 0.96,
  N_sod             = 2163,
  Sigma_tgf_delayed = 1,
  Phi_mdsod_delayed = 3.6,
  MeanArterialPressure = 100,       # mean arterial pressure


  # initial values for respiratory model
  VolAlveolar      = 2200, #150,       # alveolar volume (ml) set up to the minimum value (see FRC calculation) --> physiol. range: 140-160 ml
  ConcSystemicCO2  = 0.5357938,        # [mLSTPD/mL] Systemic tissue CO2 concentration
  ConcSystemicO2   = 0.1357565,        # [mLSTPD/mL] Systemic tissue O2 concentration
  ConcCerebralCO2  = 0.5372,           # [mLSTPD/mL] Cerebral tissue CO2 concentration
  ConcCerebralO2   = 0.112,            # [mLSTPD/mL] Cerebral tissue O2 concentration
  PressDeadspaceCO2 = 5,               # [mmHg] CO2 partial pressure lung deadspace 1
  PressDeadspaceO2  = 159,             # [mmHg] O2 partial pressure lung deadspace 1
  PressAlveolarCO2  = 40,              # [mmHg] Alveolar CO2 partial pressure
  PressAlveolarO2   = 100              # [mmHg] Alveolar O2 partial pressure
)

