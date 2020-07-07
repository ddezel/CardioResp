// Ellwein et. al preprint
// Isabelle Rudolf, Diane de Zélicourt and David Granjon

// The cardiovascular system is represented by a closed circuit with
// 3 arterial compartments, 3 venous compartments, and 
// 2 ventricular compartments -> NO atrium

#include <R.h>
#include <math.h>
#include <stdio.h>

// Load deSolve specific modules
#include "init.h"

// Load model specific modules
#include "elastances.h"
// #include "heart.h"



/*******************************************************************
 * NAME : void derivs (int *neq, double *t, double *y,
 *                    double *ydot, double *yout, int *ip)
 *
 * DESCRIPTION :     Calculates the values of the derivatives
 *
 * INPUTS :
 *       PARAMETERS:
 *           int     *neq             the number of equations
 *           double  *t               the value of the independent variable
 *           double  *y               points to a double precision array of
 *                                    length *neq that contains the current
 *                                    value of the state variables
 *           double  *ydot            points to an array that will contain
 *                                    the calculated derivatives
 *           double  *yout            points to a double precision vector whose
 *                                    first nout values are other output variables
 *                                    (different from the state variables y),
 *                                    and the next values are double precision
 *                                    values as passed by parameter rpar when
 *                                    calling the integrator. The key to the
 *                                    elements of *yout is set in *ip
 *           int     *ip              points to an integer vector whose length
 *                                    is at least 3; the first element (ip[0])
 *                                    contains the number of output values
 *                                    (which should be equal or larger than nout),
 *                                    its second element contains the length of
 *                                    *yout, and the third element contains the
 *                                    length of *ip; next are integer values,
 *                                    as passed by parameter ipar when calling
 *                                    the integrator
 */
void derivs (int *neq, double *t, double *y, double *ydot,
             double *yout, int *ip)
{
  if (ip[0] < 1) error("nout should be at least 1");
  
  
  /*------------------------Define state variables--------------------------*
   |  
   |        This part is to make equations more human readable
   |
   *------------------------------------------------------------------------*/
  
  int ny = 0;
  
  double HeartCycle;
  HeartCycle = y[ny]; ny++;
  
  // V for the cardiovascular side / 8 compartments 
  double VolPulmonaryArtery, VolPulmonaryVein, VolAorta, VolSystemicArtery, 
  VolSystemicVein, VolCerebralArtery, VolCerebralVein, VolLeftVentricle, VolRightVentricle;
  VolPulmonaryArtery    = y[ny]; ny++;
  VolPulmonaryVein      = y[ny]; ny++;
  VolAorta              = y[ny]; ny++;
  VolSystemicArtery     = y[ny]; ny++;
  VolSystemicVein       = y[ny]; ny++;
  VolCerebralArtery     = y[ny]; ny++;
  VolCerebralVein       = y[ny]; ny++;
  VolLeftVentricle      = y[ny]; ny++;
  VolRightVentricle     = y[ny]; ny++;
  
  // Define autoregulation 
  double Autoregulation;
  Autoregulation = y[ny]; ny++;
  
  // Define cardiac output
  double CardiacOutput;
  CardiacOutput     = y[ny]; ny++;
  
  // Define average aortic strain
  double MeanAorticStrain;
  MeanAorticStrain = y[ny]; ny++;
  
  // Define number of active baroreceptors
  double Baro_on;
  Baro_on = y[ny]; ny++;
  
  // Sympathetic nervous sytem activity
  double SympNervActivity;
  SympNervActivity = y[ny]; ny++;
  
  // RAAS system
  double Renin_activity, AT2_activity, N_al;
  Renin_activity = y[ny]; ny++;
  AT2_activity = y[ny]; ny++;
  N_al = y[ny]; ny++;
  
  // Vasopressin
  double N_adh;
  N_adh = y[ny]; ny++;
  
  // Sodium quantity
  double N_sod;
  N_sod = y[ny]; ny++;
  
  // Delayed terms ( only for numerical reasons)
  double Sigma_tgf_delayed, Phi_mdsod_delayed;
  Sigma_tgf_delayed = y[ny]; ny++;
  Phi_mdsod_delayed = y[ny]; ny++;
  
  // Mean arterial pressure
  double MeanArterialPressure;
  MeanArterialPressure = y[ny]; ny++;
  
  // state variables for the respiration
  double VolAlveolar, ConcSystemicCO2, ConcSystemicO2, ConcCerebralCO2, ConcCerebralO2, PressDeadspaceCO2, PressDeadspaceO2, // pD2_CO2, pD2_O2,pD3_CO2, pD3_O2, 
  PressAlveolarO2, PressAlveolarCO2;
  VolAlveolar       = y[ny]; ny++;
  ConcSystemicCO2   = y[ny]; ny++;
  ConcSystemicO2    = y[ny]; ny++; 
  ConcCerebralCO2            = y[ny]; ny++; 
  ConcCerebralO2             = y[ny]; ny++; 
  PressDeadspaceCO2 = y[ny]; ny++;
  PressDeadspaceO2  = y[ny]; ny++;
  PressAlveolarCO2  = y[ny]; ny++; 
  PressAlveolarO2   = y[ny]; ny++;  
  
  
  /*------------------------Simulations--------------------------------*
   |  Below are the conditions used in simulations: blood infusion
   |  hemorrhage, baroreflex stimulations and others
   *-------------------------------------------------------------------*/
  double Q_in = 0;
  double f_stim = 0;
  
  if (p.context == 1) {
    if (*t >= 0 && *t <= 17.5) {
      /* hemorrhage */
      Q_in = -111;
    } else {
      /* hemorrhage */
      Q_in = 0;
    };
  } else if (p.context == 2) {
    if (*t >= 0 && *t <= 5) {
      /* infusion */
      Q_in = 500;
    } else {
      /* infusion */
      Q_in = 0;
    }
  } else if (p.context == 3) {
    if (*t >= 50 && *t <= 100) {
      /* Baroreflex stimulation */
      f_stim = 1;
    } else {
      /* Baroreflex release */
      f_stim = 0;
    }
  } else if (p.context == 4) {
    // Sympathetic nervous system inihibition/activation during anesthesia
    if (*t >= 0 && *t <= 5) {
      SympNervActivity = 0.25;
    } else if (*t >= 5 && *t <= 15) {
      SympNervActivity = 0.05;
    } else if (*t >= 15 && *t <= 17.5) {
      SympNervActivity = 0.5;
    } else {
      SympNervActivity = 0.05;
    }
  }
  
  
  // Define internal parameters -> 
  double SN_baseline = 0.5;
  double H0 = 70;                   //  [bpm] mean heart rate  
  double H1 = 100;                  //  [bpm] Regulation of heart rate by SN system
  double Rmv_open = 0.000136; //0.001 / 60;     //  [mmHg min/mL] Mitral valve resistance 0.001 estimated      
  double Rav_open = 0.000068; //0.001 / 60;     //  [mmHg min/mL] Aortic valve resistance 0.001 estimated     
  double Rtv_open = 0.000136; //0.001 / 60;     //  [mmHg min/mL] Tricuspid valve resistance 0.001 estimated        
  double Rpv_open = 0.000068; //0.001 / 60;     //  [mmHg min/mL] Tricuspid valve resistance 0.001 estimated
  double Rmv_closed = 83.33333; //5000 / 60;    //  [mmHg min/mL] Mitral valve resistance 5000 estimated      
  double Rav_closed = 83.33333; //5000 / 60;    //  [mmHg min/mL] Aortic valve resistance 5000 estimated     
  double Rtv_closed = 83.33333; //5000 / 60;    //  [mmHg min/mL] Tricuspid valve resistance 5000 estimated        
  double Rpv_closed = 83.33333; //5000 / 60;    //  [mmHg min/mL] Tricuspid valve resistance 5000 estimated
 
  double Cpa = 3.5;       //4.3;           //  [mL/mmHg] Pulmonary arterial compliance 
  double Cpv = 78.9;      //64;            //  [mL/mmHg] Pulmonary venous compliance 
  double Csa_0 = 9.466;   //21;   //24;    //  [mL/mmHg] Systemic arteries compliance [baseline arterial compliance was 2.56 ± 0.18 ]
  double Csv_0 = 377.55;  //385;  //413;   //  [mL/mmHg] Systemic venous compliance 
  double Cca =  1;        //2.96;          //  [mL/mmHg] Cerebral arterial compliance 
  double Ccv = 35.48;     //30;            //  [mL/mmHg] Cerebral venous compliance

  double Caorta_0 = 1.5;      //1.75             // [mL/mmHg] Aortic compliance 
  double Raorta_0 = 0.01394; //0.0004166667;     // [mmHg min/mL] Aortic resistance (added /60 to correct Psa)
  //double Rs_0   = 0.0014;                     // [mmHg min/mL] Nominal systemic arterial resistance

  /*------------------------Equations for the Heart-------------------------*
   |  Below are the equations from the modified Ellwein 2016 model
   |  This model is not officially published but is promising and well 
   |  documented.
   *------------------------------------------------------------------------*/
  
  //  previous parameter values:
  // ---------------------------------------------------------
  // resistances
  // ---------------------------------------------------------
   
  // Adding CO2 regulation to Rc and Rs (coupling Respiratory and Cardiovascular Model)
  // cf Ellwein PhD Thesis Chapter 7.2.2 Regulation p. 139 and 65
  double Rc_pB_CO2, Rc_max, Rc_min, Pb_0, Pb_CO2, Pb_CO2_goal, Rc_s, k_c;
  double Pa_0, Pa_goal, Rs_s, Rs_min, Rs_max, k_s;  //Rs_pam;
  double Rs_pa_CO2, pa_CO2_goal;
  
  // calculate ConcCerebralCO2 in advance for resistance calculcations
  double ca_CO2, ca_O2, cv_CO2, cv_O2;
  ca_CO2 = p.KCO2 * PressAlveolarCO2 + p.kCO2;               // ~ 0.493 
  ca_O2 = p.K1 * pow((1- exp(-p.K2 * PressAlveolarO2)), 2);  // ~ 0.197
  
  // values used for Rc_min/max and Rs_min/max are from Ellwein phD thesis p. 145
  // the values commented out are from calculations based on estimations
  
  Rc_min = 0.0578; //0.006; //0.0658;              // 1.829 / 60;    // minimal cerebral resistance [mmHg*min/mL]
  Rc_max = 0.1735; //0.06; //0.1991667;           // 3.251 / 60;    // maximal cerebral resistance [mmHg*min/mL]  
  Rc_s = 0.5 * (Rc_min + Rc_max);   // Cerebral Systemic Resistance [mmHg*min/mL]
  Pb_CO2 = 40;                  // Partial Pressure of CO2 in brain [mmHg]
  Pb_CO2_goal = 40;             // Goal Partial Pressure of CO2 in brain [mmHg]
  pa_CO2_goal = 40;             // Arterial Partial Pressure CO2 Goal [mmHg]
  k_c = 50;                     // Fixed parameter for cerebral part, given by Ellweins phD thesis    
  Pa_goal = 85;                 // Arterial Pressure Goal [mmHg]
  Rs_min = 0.00453; //0.002319047; //=0.01623333/ 7;          // 0.01125;       // Min Systemic Resistance [mmHg*min/mL]
  Rs_max = 0.01361; //0.005621429; //=0.03935/ 7;             // 0.020  ;       // Max Systemic Resistance [mmHg*min/mL] 
  Rs_s = 0.5 * (Rs_min + Rs_max);    // Systemic Resistance [mmHg*min/mL]
  k_s = 5;                      // Fixed parameter for systemic part, given by Ellweins phD thesis 
  
  Pb_0 = 0.5; // Pb_CO2_goal* pow( (Rc_s - Rc_min) / (Rc_max - Rc_s), (1/k_c) ); // replace by its nominal value
  //Pa_0 = Pa_goal * pow( (Rs_s - Rs_min) / (Rs_max - Rs_s), (1/k_s) );
  Pa_0 = 0.5; // pa_CO2_goal * pow( (Rs_max - Rs_s) / (Rs_s - Rs_min), (1/k_s) ); // replace by its nominal value
  
  Rc_pB_CO2 = (Rc_max - Rc_min) * (Pb_0 / (ca_CO2 + Pb_0)) + Rc_min;  
  // Rs_pam = (Rs_max - Rs_min) * (Pa_0 / (MeanArterialPressure + Pa_0)) + Rs_min; //-> not necessary, MeanArterialPressure is already regulated by kidney
  Rs_pa_CO2 = (Rs_max - Rs_min) * (ca_CO2 / (ca_CO2 + Pa_0)) + Rs_min; 
  
  
  double Rca, Rcv, Rbrain, Rtot;
  double Rs, Rp, Raorta, Rc, Raa; 
  double Rs_coeff;
  Rs_coeff = 1 + p.alpha2 * SN_baseline;
  
  Rp    = 0.004;        // DDZ 0.002;     //0.004; // [mmHg min/ml] Pulmonary resistance --> phys. range: 0.00025-0.0016 mmHg min/ml (wiki), 1/7 of Rs
  Raorta = Raorta_0 * (1 + 0.5 * (SympNervActivity-SN_baseline)/SN_baseline); // DDZ: the 1-phi without a multiplicative factor is kind of awkward // * phi_periph) + (1 - phi_periph));
  Raa    = Raorta * Autoregulation;
  //Rs     = Rs_pa_CO2/Rs_coeff * (1 + p.alpha2 * (SympNervActivity)); 
  Rs     = Rs_s/Rs_coeff * (1 + p.alpha2 * (SympNervActivity)); // DDZ --> if the effect of CO2 is included in SympNervActivity, then we should not have it the base Rs calculation
  
  Rca = 0.20 * Rc_pB_CO2;       //1e-005; // [mmHg min/ml] Cerebral arterial resistance (Psa − Pca)/Qca
  Rc  = 0.79 * Rc_pB_CO2;       //0.092;  // [mmHg min/ml] Cerebral resistance (Pca − Pcv)/Qc
  Rcv = 0.01 * Rc_pB_CO2;       //4e-005; // [mmHg min/ml] Cerebral venous resistance (Pcv − Psv)/Qcv
  Rbrain = Rca + Rc + Rcv;                // [mmHg min/ml] total brain resistance
  
  // autoregulation
  double Autoregulation_inf;
  Autoregulation_inf = 0.5 * (1 + tanh((CardiacOutput- p.CO_0) / p.CO_1));
  
  // vascular compliances
  // ---------------------------------------------------------
  double Csa, Csv, Caorta;
  Csa = Csa_0 / (1 + p.alpha1 * SympNervActivity); // * (1 + p.alpha3 * AT2_activity));  
  Csv = Csv_0 / (1 + p.alpha1 * SympNervActivity); // * (1 + p.alpha3 * AT2_activity)); 
  Caorta = Caorta_0; // / (1 + p.alpha1 * SympNervActivity); // * (1 + p.alpha3 * AT2_activity));  
  
  // Pressures
  double Ppa, Ppv, Paorta, Psa, Psv, Pca, Pcv, Plv, Prv;
  Ppa = VolPulmonaryArtery / Cpa;
  Ppv = VolPulmonaryVein / Cpv; 
  Paorta = VolAorta / Caorta;
  Psa = VolSystemicArtery / Csa;
  Psv = VolSystemicVein / Csv; 
  Pca = VolCerebralArtery / Cca;
  Pcv = VolCerebralVein / Ccv;
  
  Plv = TimeVaryingElastance_L(HeartCycle, VolLeftVentricle, SympNervActivity, SN_baseline); 
  Prv = TimeVaryingElastance_R(HeartCycle, VolRightVentricle, SympNervActivity, SN_baseline);
  
  // Update ventricular heart valves resistances (open/closed/regurgitation)
  // k = 2, is the speed of the transition from open to closed
  // k = 10 in the original PhD thesis
  double Rav, Rmv, Rpv, Rtv;
  Rav = min(Rav_open + exp(-2 * (Plv - Paorta)), Rav_closed);
  Rmv = min(Rmv_open + exp(-2 * (Ppv - Plv)), Rmv_closed);
  Rpv = min(Rpv_open + exp(-2 * (Prv - Ppa)), Rpv_closed);
  Rtv = min(Rtv_open + exp(-2 * (Psv - Prv)), Rtv_closed);  
  
  Rav = Rav_open;
  Rmv = Rmv_open;
  Rpv = Rpv_open;
  Rtv = Rtv_open;
  
  // Update fluxes (flow rates between compartments)
  double Qtv, Qpv, Qmv, Qav, Qaa, Qp, Qs, Qca, Qc, Qcv;
  Qtv = max(0, (Psv - Prv) / Rtv);     // [mL/min] flow rate through the tricuspid valve 
  Qpv = max(0, (Prv - Ppa) / Rpv);     // [mL/min] flow rate through the pulmonary valve 
  Qmv = max(0, (Ppv - Plv) / Rmv);     // [mL/min] flow rate through the mitral valve 
  Qav = max(0, (Plv - Paorta) / Rav);  // [mL/min] flow rate through the aortic valve 
  
  Qaa = (Paorta - Psa) / Raa;   // [mL/min] flow rate through the aorta to systemic arteries
  Qp  = (Ppa - Ppv) / Rp;                        // [mL/min] flow rate through the pulmonary vascular bed, 82ml/s oder 4920 ml/min
  Qs  = (Psa - Psv) / Rs;                        // [mL/min] flow rate through the systemic vascular bed
  
  Qca = (Paorta - Pca) / Rca;          // [mL/min] flow rate through the cerebral arteries
  Qc  = (Pca - Pcv) / Rc;           // [mL/min] flow rate through the cerebral vascular bed
  Qcv = (Pcv - Psv) / Rcv;          // [mL/min] flow rate through the cerebral veins
  
  // Total blood volume
  double Vtot;
  Vtot = VolPulmonaryArtery + VolPulmonaryVein + VolAorta + VolSystemicArtery + 
    VolSystemicVein + VolLeftVentricle + VolRightVentricle + VolCerebralArtery + VolCerebralVein;
  
  // Baroreflex
  double A_aorta_0, A_aorta, Aortic_strain, delta_strain, f_Baro;
  A_aorta_0 = pi * pow(p.d_aorta_0, 2) / 4;
  A_aorta = 1000 * VolAorta / p.L_aorta; // convert from mL to mm3 since area is in mm2 and V in mL
  Aortic_strain = (A_aorta - A_aorta_0) / A_aorta_0;
  delta_strain = max(Aortic_strain - MeanAorticStrain, 0);
  
  f_Baro = p.f_0 * Baro_on * delta_strain / (delta_strain + p.delta_strain_0);
  
  // Heart rate [bpm]
  double HR;
  HR = H0 + H1 * (SympNervActivity - SN_baseline)/SN_baseline; //* (SympNervActivity - 0.25);
  
  /*------------------------Equations for the kidney-------------------*
   |  Below are equations of the model of Karaaslan et al., 2005
   |  It contains only the kidney process, since the heart equations
   |  are taken from another model (Beard). 
   *-------------------------------------------------------------------*/
  double g, Renin_activity_inf, C_sod, C_al, alpha_map, alpha_rap, rsna, beta_rsna, 
  R_aa, R_ea, R_r, Phi_rb,
  P_gh, P_f, Phi_gfilt, Sigma_tgf, Phi_filsod, gamma_filsod, gamma_at,
  gamma_rsna, eta_ptsodreab, Phi_ptsodreab, Phi_mdsod, Psi_al, eta_dtsodreab,
  Phi_dtsodreab, Phi_dtsod, lambda_dt, lambda_anp, eta_cdsodreab, Phi_cdsodreab,
  Phi_usod, Phi_sodin, C_anp, C_adh, C_k, mu_al, mu_adh, Phi_twreab, Phi_u,
  xi_ksod, xi_at, N_als, xi_map, V_ecf, C_at, N_adhs, Phi_win, P_ra, CO_scaled,
  epsilon_CO_scaled, epsilon_Pra, epsilon_Pma, log10_Cat, epsilon_At_Rea,
  epsilon_Phimdsod_tgf, epsilon_Philsod, epsilon_At_Ptsodreab, epsilon_rsna,
  log10_Cal, epsilon_Al_Dtsodreab, epsilon_dtsod_Cdsodreab, pow_ADH, log10_Cadh,
  epsilon_Al_wreab, epsilon_ADH_wreab, epsilon_Pma_Nals, epsilon_At_Nals;
  
  // define internal parameters (kidney)
  double N_rsna = 1;           // Normalized renal sympathetic nerve activity
  double R_aass = 31.67;       //
  double P_B = 18;             // [mmHg] Bowman hydrostatic pressure
  double P_go = 28;            // [mmHg] Glomerular osmotic pressure
  double C_gcf = 0.00781;      // Glomerular capillary filtration coefficient
  double n_etapt = 0.8;        // Fractional proximal sodium reabsorption
  double n_epsilondt = 0.5;    // Fractional distal tubule sodium reabsorption
  double n_etacd = 0.93;       // Fractional collecting duct sodium reabsorption
  
  // ECF in L (convert from mL to L as blood volume is in mL)                                                          
  V_ecf = Vtot * 2.76 /1000; 
  
  // Renin
  g = 4 * (100 - p.P1 - 0.737 * p.P2);
  Renin_activity_inf = (1 - tanh((MeanArterialPressure -g * SympNervActivity - p.P1) / p.P2)) / 2;
  
  // we also scale the angiotensine concentration. y[14] is the normalized
  // ATII concentration returned by the cardiac submodel (0.19). C_at is 
  // 20ng/l in Karaaslan model 2005. Need a factor 105.3 to scale.
  C_at = AT2_activity * 105.3;
  
  // Calculation of sodium, aldosterone and ADH concentrations
  C_sod = N_sod / V_ecf; 
  C_al = N_al * 85;
  C_adh = 4 * N_adh;
  
  // the following terms representing sodium intake, normalized ANP concentration, 
  // potassium concentration are taken constants. 
  
  Phi_sodin = 0.126;
  C_k = 5;
  
  // normalize atrial natriuretic peptide concentration (C_anp), 
  //as a function of right atrial pressure (P_ra)
  CO_scaled = CardiacOutput/ 1000; // convert from mL/min to L/min
  epsilon_CO_scaled = exp(CO_scaled * 0.2281);
  P_ra = 0.2787 * epsilon_CO_scaled; // right atrial pressure in Karaaslan
  epsilon_Pra = exp((P_ra - 3.762) / 1);
  C_anp = 7.427 - 6.554 / (1 + epsilon_Pra);
  
  // Block 1: Renal sympathetic nerve activity
  // decreased by P_ma (arterial pressure) as well as
  // P_ra (right atrial pressure)
  // No effect of Angiotensin (assumed)
  epsilon_Pma = exp((MeanArterialPressure - 100) / 15);
  alpha_map = 0.5 + 1.1 / (1 + epsilon_Pma);
  alpha_rap = 1 - 0.08 * P_ra;
  rsna = N_rsna * alpha_map * alpha_rap;      // Renal sympathetic nerve activity
  
  /* rsna stimulation  from Karaaslan et al. 2005 */
  //if (*t >= 0 && *t <= 17680) {
  //  rsna = 1.4;
  //} else {
  //  rsna = N_rsna * alpha_map * alpha_rap;
  //}
  
  /* renal denervation Karaaslan et al. 2005 */
  //rsna = 1;
  
  // Block 2: Renal vascular resistance, that is 
  // sum of afferent (R_aa) and efferent (R_ea) arteriolar resistances.
  // Rea is a function of Angiotensine, Raa a constance times effects of 
  // rsna and tubuloglomerular feedback (tgf). Since, Sigma_tgf is defined after,
  // I replaced it by y[17] which is Sigma_tgf_delayed, defined in the ODE part.
  beta_rsna = 1.5 * (rsna - 1) + 1;
  R_aa = 1.038414 * R_aass * beta_rsna * Sigma_tgf_delayed; // add 0.93 to correct the pressure
  log10_Cat = log10(C_at);
  epsilon_At_Rea = exp(3.108 - 1.785 * log10_Cat);
  R_ea = (0.9432 + 0.1363 / (0.2069 + epsilon_At_Rea)) * 51.66;
  R_r = R_aa + R_ea;
  Rtot = 1 / (1 / Rbrain + 1 / (Raa + Rs) );   // total body resistance
  
  // Block3: Renal blood flow, which is mean arterial pressure (MeanArterialPressure) 
  // over renal vascular resistance (R_r)
  Phi_rb = MeanArterialPressure / R_r;
  
  // renal artery constriction, does not work
  //if (*t >= 0 && *t <= 1440) {
  //  Phi_rb = 0;
  //} else {
  //  Phi_rb = MeanArterialPressure / R_r;
  //}
  
  
  // Block 4: Glomerular filtration rate (Phi_gfilt), product of net 
  // filtration pressure (P_f), times glomerular filtration constant (C_gcf).
  // P_gh is the glomerular hydrostatic pressure, P_B the hydrostatic 
  // pressure in Bowman, P_go the glomerular osmotic pressure. P_gh is the 
  // difference between mean arterial pressure and afferent arteriolar 
  // resistance
  P_gh = MeanArterialPressure - Phi_rb * R_aa;
  P_f = P_gh - (P_B + P_go);
  Phi_gfilt = P_f * C_gcf;
  
  /* Simulations of GFR decrease*/
  //if (*t <=50 && *t < 300) {
  //  Phi_gfilt = P_f * C_gcf;
  //} else {
  //  Phi_gfilt = 0.2 * P_f * C_gcf;
  //}
  
  
  // Block 5: Tubuloglomerular feedback (Sigma_tgf), function of the macula 
  // densa sodium flow (Phi_mdsod). As the equation of Phi_mdsod is defined
  // after, I replace it by a delayed term y[18] that is Phi_mdsod_delayed.
  epsilon_Phimdsod_tgf = exp((Phi_mdsod_delayed - 3.859) / -0.9617);
  Sigma_tgf = 0.3408 + 3.449 / (3.88 + epsilon_Phimdsod_tgf);
  
  // Block 6: Filtered sodium load (Phi_filsod)
  Phi_filsod = Phi_gfilt * C_sod;
  
  // Block 7: Sodium reabs in proximal tubule (Phi_ptsodreab), product of 
  // what is filtered (Phi_filsod) and a fraction coefficient (eta_ptsodreab).
  // eta_ptsodreab is a function of filtered sodium load (gamma_filsod), 
  // angiotensine (gamma_at) and rsna (gamma_rsna). Its normal value is set 
  // to 0.8 (n_etapt)
  epsilon_Philsod = exp(Phi_filsod - 14);
  gamma_filsod = 0.8 + 0.3 / (1 + epsilon_Philsod / 138);
  epsilon_At_Ptsodreab = exp(2.6 - 1.8 * log10_Cat);
  gamma_at = 0.95 + 0.12 / (1 + epsilon_At_Ptsodreab);
  epsilon_rsna = exp(1 - rsna);
  gamma_rsna = 0.5 + 0.7 / (1 + epsilon_rsna / 2.18);
  eta_ptsodreab = n_etapt * gamma_filsod * gamma_at * gamma_rsna * 0.9952413;
  Phi_ptsodreab = Phi_filsod * eta_ptsodreab;
  
  // Block 8: Macula densa sodium flow (Phi_mdsod), difference between what
  // is filtered and reabsorbed in the proximal tubule
  Phi_mdsod = Phi_filsod - Phi_ptsodreab;
  
  // Block 9: Sodium reabs in distal tubule (Phi_dtsodreab) product of 
  // what is filtered in the macula densa (Phi_mdsod) and a 
  // fraction coefficient (eta_dtsodreab). eta_dtsodreab is affected by 
  // aldosterone concentration (Psi_al).
  log10_Cal = log10(C_al);
  epsilon_Al_Dtsodreab = exp((0.48 - 1.2 * log10_Cal) / 0.88);
  Psi_al = 0.17 + 0.94 / (1 + epsilon_Al_Dtsodreab);
  eta_dtsodreab = n_epsilondt * Psi_al * 0.9939141;
  Phi_dtsodreab = Phi_mdsod * eta_dtsodreab;
  
  // Block 10: Distal tubule sodium outflow (Phi_dfsod), difference between what
  // was in macula densa and what is reabsorbed in the distal tubule.
  Phi_dtsod = Phi_mdsod - Phi_dtsodreab;
  
  // Block 11: Sodium reabs in the collecting duct (Phi_cdsodreab), product of 
  // what is filtered in the distal tubule (Phi_dtsod) and a 
  // fraction coefficient (eta_cdsodreab). eta_cdsodreab is a function of its
  // normal value times the effect of distal tubule sodium outflow (lambda_dt)
  // and the effect of the atrial natriuretic peptide (lambda_anp)
  epsilon_dtsod_Cdsodreab = exp((Phi_dtsod - 1.6) / 2);
  lambda_dt = 0.82 + 0.39 / (1 + epsilon_dtsod_Cdsodreab);
  lambda_anp = -0.1 * C_anp + 1.1199;
  eta_cdsodreab = n_etacd * lambda_dt * lambda_anp * 0.9753;
  Phi_cdsodreab = Phi_dtsod * eta_cdsodreab;
  
  // Block 12:  Urine sodium flow (Phi_usod), equals the sodium that is in the
  // collecting duct less what has been reabsorbed
  Phi_usod = Phi_dtsod - Phi_cdsodreab;
  
  // Block 13: Water intake (Phi_win), a function of vasopressin concentration
  // (C_adh). Multiply by 0.77 to correct for steady-state value
  // add max to prevent negative values
  pow_ADH = pow(C_adh, -1.607);
  Phi_win = max(0, (0.008 / (1 + 1.822 * pow_ADH) - 0.0053) * 0.77 * 0.9363);
  
  // Dehydration example
  //if (*t >= 1440 && *t <= 2880) {
  //  Phi_win = 0;
  //} else {
  //  Phi_win = (0.008 / (1 + 1.822 * pow_ADH) - 0.0053) * 0.77 * 0.9363; 
  //}
  
  // block 27: Tubular water reabsorption rate (Phi_twreab), affected by 
  // alsdosterone (mu_al) and vasopressin concentration (mu_adh).
  log10_Cadh = log10(C_adh);
  epsilon_Al_wreab = exp((0.48 - 1.2 * log10_Cal) / 0.88);
  epsilon_ADH_wreab = exp(0.6 - 3.7 * log10_Cadh);
  mu_al = 0.17 + 0.94 / (1 + epsilon_Al_wreab);
  mu_adh = 0.37 + 0.8 / (1 + epsilon_ADH_wreab);
  Phi_twreab = 0.025 - 0.001 / (mu_al * mu_adh) + 0.8 * Phi_gfilt;
  
  // Block 28: urine flow rate (Phi_u). Difference between GFR (Phi_gfilt) and
  // the tubular water reabsorption (Phi_twreab).
  Phi_u = max(0, Phi_gfilt - Phi_twreab);
  
  // Block 34: Aldosterone concentration (N_al). Effects of potassium/sodium ratio
  // (xi_ksod), mean arterial pressure (xi_map) and angiotensin concentration
  // (xi_at).
  
  //# the following equations is corrected for steady-state value (0.80). Moreover,
  // I set a max(0, equation) to ensure that xi_ksod remains positive. It was
  // a problem simulating the hemorrhage of Beard et al., trigering numerical 
  // errors
  xi_ksod = max(0, ((C_k / C_sod) / 0.003525 - 9) / 0.85); 
  if (MeanArterialPressure <= 100) {
    epsilon_Pma_Nals = exp(-0.0425 * MeanArterialPressure);
    xi_map = 69.03 * epsilon_Pma_Nals;
  } else {
    xi_map = 1;
  }
  epsilon_At_Nals = exp((2.82 - 1.5 * log10_Cat) / 0.8);
  xi_at = 0.4 + 2.4 / (1 + epsilon_At_Nals);
  N_als = xi_ksod * xi_map * xi_at;
  
  // Block 26: Vasopressin concentration (C_adh). Secretion rate is mainly 
  // controlled by sodium concentration (C_sod), autonomic multiplier effect
  // (epsilon_aum) and the right atrial pressure (delta_ra).
  // find a solution for delta_ra
  N_adhs = max(0, (C_sod - 141 + 2.42) / 3); 
  // I removed the effects of rsna and Pra and add a term to correct for the initial decrease
  
  
  /*------------------------ODEs-------------------------*
   |
   |               Below are ODEs equations
   |
   *------------------------------------------------------*/
  
  // fraction of the cardiac cycle (normalized time)
  double dHeartCycle;
  dHeartCycle = HR;
  
  // Update volumes in each compartments
  double d_VolPulmonaryArtery, d_VolPulmonaryVein, d_VolAorta, d_VolSystemicArtery, 
  d_VolSystemicVein, d_VolCerebralArtery, d_VolCerebralVein, d_VolLeftVentricle, d_VolRightVentricle;
  d_VolPulmonaryArtery = Qpv - Qp;  
  d_VolPulmonaryVein   = Qp  - Qmv;
  d_VolAorta           = Qav - Qaa;
  d_VolSystemicArtery  = Qaa - Qs - Qca;
  d_VolSystemicVein    = Qs - Qtv + Qcv; //+ Phi_win - Phi_u + Q_in
  d_VolCerebralArtery  = Qca - Qc;
  d_VolCerebralVein    = Qc  - Qcv;
  d_VolLeftVentricle   = Qmv - Qav;
  d_VolRightVentricle  = Qtv - Qpv;
  
  // autoregulation
  double d_Autoregulation;
  d_Autoregulation = 1 / p.tau_autoreg * (Autoregulation_inf - Autoregulation);

  // Cardiac Output 
  double d_CardiacOutput;
  // d_CardiacOutput= 1 / p.tau_CardiacOutput* ((Paorta - Psa) / (Raorta * Autoregulation) - CO);
  d_CardiacOutput= 1 / p.tau_CO* (Qaa - CardiacOutput);  //simplification
  
  // Mean arterial pressure
  double d_MeanArterialPressure;
  d_MeanArterialPressure = 1 / 0.25 * (Paorta - MeanArterialPressure);
  
  // Mean aortic strain
  double d_MeanAorticStrain;
  d_MeanAorticStrain = 1 / p.tau_strain * (Aortic_strain - MeanAorticStrain);
  
  // Active baroreceptors dynamic
  double d_Baro_on;
  d_Baro_on = p.a * (1 - Baro_on) - p.b * f_Baro / p.f_0;
  
  // SN system
  double d_SympNervActivity;
  if (p.context == 4) {
    // in case of anesthesia
    d_SympNervActivity = 0;
  } else {
    // otherwise
    //d_SympNervActivity = p.f_SN * (1 - SympNervActivity) - (f_Baro + f_stim) * SympNervActivity;
    //d_SympNervActivity = p.f_SN * (0.5 - SympNervActivity) + 0.08333333 * ((ca_CO2 * SympNervActivity));   // 0.083 = 1/12
    d_SympNervActivity = p.f_SN * (ca_CO2 - SympNervActivity); // + 0.15 * (ca_CO2 - Pa_0_1) / Pa_0_1* SympNervActivity;
  }
  
  // RAAS sytem
  double d_Renin_activity, d_AT2_activity, d_N_al;
  d_Renin_activity     = (Renin_activity_inf - Renin_activity) / p.tau_renin;        
  d_AT2_activity       = (Renin_activity - AT2_activity) / p.tau_at;                
  d_N_al               = 1 / p.tau_al * (N_als - N_al);                                  
  
  // Vasopressin
  double d_N_adh;
  d_N_adh              = 1 / p.tau_adh * (N_adhs - N_adh);                               
  
  // Sodium homeostasis
  double d_N_sod;
  d_N_sod              = Phi_sodin - Phi_usod;                                     
  
  // needed to use Sigma_tgf and Phi_mdsod before they are defined
  double d_Sigma_tgf_delayed, d_Phi_mdsod_delayed;
  d_Sigma_tgf_delayed  = Sigma_tgf - Sigma_tgf_delayed;                            
  d_Phi_mdsod_delayed  = Phi_mdsod - Phi_mdsod_delayed;                            
  
  
  /*--------------------- Respiratory Model ---------------------*
   * 
   *-------------------------------------------------------------*/
  
  cv_CO2 = ( Qs * ConcSystemicCO2 + Qc * ConcCerebralCO2 )  / (Qs + Qc);   //ConcSystemicCO2 [mlSTPD/ml] ~ 0.562       
  cv_O2  = ( Qs * ConcSystemicO2  + Qc * ConcCerebralO2  )  / (Qs + Qc);   //ConcSystemicO2  [mlSTPD/ml] ~ 0.117    
  
  // ventilation process
  //---------------------------------------------
  // Ve = ventilation volume (mL)  --> d_Ve = ventilation rate [mL/min]  ~ ca 6'000 mL/min (until 9000 mL/min)
  // VolAlveolar= alveolar volume (mL)     --> d_VolAlveolar= alveolar ventilation rate [mL/min] (approximated as ~d_Ve)
  double d_Ve, d_VolAlveolar;
  //double tidal_volume;    // [mL]
  // double resp2CC_ratio;   // ratio of the respiratory cycle to cardiac cycle duration 
  // resp2CC_ratio = 0.2;    // slow breathing: 0.08; //hyperventilation: 0.8; 
  //tidal_volume  = 600;    // hyperventilation: 800; 
  
  // simulate hyperventilation after 1000seconds
  //if (HeartCycle < 2000) {
  //  resp2CC_ratio = 0.2;   // slow breathing: 0.08; //hyperventilation: 0.8; 
  //  tidal_volume  = 600;    // hyperventilation: 800; 
  //}
  //else {
  //  resp2CC_ratio = 0.8;   // slow breathing: 0.08; //hyperventilation: 0.8; 
  //  tidal_volume  = 800;    // hyperventilation: 800; 
  //}
  
  d_Ve = (2 * pi * HR * p.resp2CC_ratio) * 0.5 * p.Vtid * sin(2 * pi * HeartCycle * p.resp2CC_ratio);
  d_VolAlveolar= d_Ve;
  
  // intrathoracic pressure
  // dpit = - 1.8 * d_VolAlveolar/ Va_diff;
  
  // gas concentrations
  double d_ConcSystemicCO2, d_ConcSystemicO2, d_ConcCerebralCO2, d_ConcCerebralO2;
  d_ConcSystemicCO2 = ( p.MS_CO2 + Qs * (ca_CO2 - ConcSystemicCO2) ) / p.VS_CO2;   //ConcSystemicCO2 [mlSTPD/ml] ~ 0.562       
  d_ConcSystemicO2  = (-p.MS_O2  + Qs * (ca_O2  - ConcSystemicO2 )) / p.VS_O2;    //ConcSystemicO2  [mlSTPD/ml] ~ 0.117    
  d_ConcCerebralCO2 = ( p.MB_CO2 + Qc * (ca_CO2 - ConcCerebralCO2) ) / p.VB_CO2;   //ConcCerebralCO2 [mlSTPD/ml] ~ 0.547       
  d_ConcCerebralO2  = (-p.MB_O2  + Qc * (ca_O2  - ConcCerebralO2 )) / p.VB_O2;    //ConcCerebralO2  [mlSTPD/ml] ~ 0.157   
  
  
  // Inspiration/expiration cycle
  double d_PressDeadspaceCO2, d_PressDeadspaceO2, d_PressAlveolarCO2, d_PressAlveolarO2; //d_pD2_CO2, d_pD2_O2, d_pD3_CO2, d_pD3_O2, d_PressAlveolarCO2, d_PressAlveolarO2;
  
  // conversion from gas fractions to partial pressure (for an ambient P=760mmHg satured --> 47mmHg partial P water)
  // gas_fraction_A = partial_pressure_A/(ambient_pressure − partial_pressure_water) 
  //                = partial_pressure_A/(760-47)
  //                = partial_pressure_A//713
  // or partial_pressure_A = gas_fraction_A * 713
  
  // conversion STPD (T=0°C=273.15°K and P=760mmHg) to BTPS (T=37°C=273.15+37°K=310.15°K and Pdry=Psaturated-Pwater=760-47=713mmHg)
  // V = nRT/P
  // so Vbtps = nR(310/713)  and Vstpd = nR(273/760)
  // Vstpd = 713/310*273/760 * Vbtps = 0.8261842 * Vbtps
  // Vbtps = 760/273*310/713 * Vstpd = 1.210653753 * Vstpd
  
  double pulmShunting;      // portion of blood shunting the pulmonary circulation ~ 2%
  double Vstpd2Vbtps;
  double gasFrac2partialP;
  
  pulmShunting     = 0.02;
  Vstpd2Vbtps      = 1.210653753;
  gasFrac2partialP = 713;
  
  if (d_VolAlveolar> 0) {
    // Inspiration with one deadspace
    d_PressDeadspaceCO2 = d_Ve * (p.pi_CO2 - PressDeadspaceCO2) / p.VD;
    d_PressDeadspaceO2  = d_Ve * (p.pi_O2  - PressDeadspaceO2)  / p.VD;
    
    d_PressAlveolarCO2  = ( ((1-pulmShunting) * gasFrac2partialP * Vstpd2Vbtps) * Qp * (cv_CO2 - ca_CO2) + d_Ve * (PressDeadspaceCO2 - PressAlveolarCO2) ) / VolAlveolar; 
    d_PressAlveolarO2   = ( ((1-pulmShunting) * gasFrac2partialP * Vstpd2Vbtps) * Qp * (cv_O2  - ca_O2 ) + d_Ve * (PressDeadspaceO2  - PressAlveolarO2 ) ) / VolAlveolar;
    
    //d_PressDeadspaceCO2 = d_Ve * (p.pi_CO2 - PressDeadspaceCO2) / VD_3;      
    //d_PressDeadspaceO2  = d_Ve * (p.pi_O2  - PressDeadspaceO2)  / VD_3;
    //d_pD2_CO2 = d_Ve * (PressDeadspaceCO2  - pD2_CO2) / VD_3;
    //d_pD2_O2  = d_Ve * (PressDeadspaceO2   - pD2_O2)  / VD_3;
    //d_pD3_CO2 = d_Ve * (pD2_CO2  - pD3_CO2) / VD_3;
    //d_pD3_O2  = d_Ve * (pD2_O2   - pD3_O2)  / VD_3;
    
  } else {
    // Expiration with one deadspace
    d_PressDeadspaceCO2 = d_Ve * (PressDeadspaceCO2 - PressAlveolarCO2) / p.VD;
    d_PressDeadspaceO2  = d_Ve * (PressDeadspaceO2  - PressAlveolarO2)  / p.VD;
    
    //d_PressDeadspaceCO2 = d_Ve * (PressDeadspaceCO2 - pD2_CO2) / VD_3;
    //d_PressDeadspaceO2  = d_Ve * (PressDeadspaceO2  - pD2_O2)  / VD_3;
    //d_pD2_CO2 = d_Ve * (pD2_CO2 - pD3_CO2) / VD_3;
    //d_pD2_O2  = d_Ve * (pD2_O2  - pD3_O2)  / VD_3;
    //d_pD3_CO2 = d_Ve * (pD3_CO2 - PressAlveolarCO2)  / VD_3;
    //d_pD3_O2  = d_Ve * (pD3_O2  - PressAlveolarO2)   / VD_3;
    
    d_PressAlveolarCO2  = ( ((1-pulmShunting) * gasFrac2partialP * Vstpd2Vbtps) * Qp * (cv_CO2 - ca_CO2) ) / VolAlveolar; 
    d_PressAlveolarO2   = ( ((1-pulmShunting) * gasFrac2partialP * Vstpd2Vbtps) * Qp * (cv_O2  - ca_O2 ) ) / VolAlveolar;
  }
  
  
  // ydot is the array for ODE output
  int nydot = 0;
  
  ydot[nydot] = HeartCycle; nydot++;
  ydot[nydot] = d_VolPulmonaryArtery; nydot++;
  ydot[nydot] = d_VolPulmonaryVein;   nydot++;
  ydot[nydot] = d_VolAorta; nydot++;
  ydot[nydot] = d_VolSystemicArtery;  nydot++;
  ydot[nydot] = d_VolSystemicVein;    nydot++;
  ydot[nydot] = d_VolCerebralArtery;  nydot++;
  ydot[nydot] = d_VolCerebralVein;    nydot++;
  ydot[nydot] = d_VolLeftVentricle;   nydot++;
  ydot[nydot] = d_VolRightVentricle;  nydot++;
  ydot[nydot] = d_Autoregulation;     nydot++;
  ydot[nydot] = d_CardiacOutput;      nydot++;
  ydot[nydot] = d_MeanAorticStrain;   nydot++;
  ydot[nydot] = d_Baro_on;            nydot++;
  ydot[nydot] = d_SympNervActivity;   nydot++;
  ydot[nydot] = d_Renin_activity;     nydot++;
  ydot[nydot] = d_AT2_activity;       nydot++;
  ydot[nydot] = d_N_al;               nydot++;
  ydot[nydot] = d_N_adh;              nydot++;
  ydot[nydot] = d_N_sod;              nydot++;
  ydot[nydot] = d_Sigma_tgf_delayed;  nydot++;
  ydot[nydot] = d_Phi_mdsod_delayed;  nydot++;
  ydot[nydot] = d_MeanArterialPressure; nydot++;
  
  // respiratory output
  ydot[nydot] = d_VolAlveolar;        nydot++;
  ydot[nydot] = d_ConcSystemicCO2 ;   nydot++;
  ydot[nydot] = d_ConcSystemicO2 ;    nydot++;
  ydot[nydot] = d_ConcCerebralCO2;    nydot++;
  ydot[nydot] = d_ConcCerebralO2;     nydot++;
  ydot[nydot] = d_PressDeadspaceCO2;  nydot++;
  ydot[nydot] = d_PressDeadspaceO2;   nydot++;
  //ydot[nydot] = d_pD2_CO2;          nydot++;
  //ydot[nydot] = d_pD2_O2;           nydot++;
  //ydot[nydot] = d_pD3_CO2;          nydot++;
  //ydot[nydot] = d_pD3_O2;           nydot++;
  ydot[nydot] = d_PressAlveolarCO2;   nydot++;
  ydot[nydot] = d_PressAlveolarO2;    nydot++;
  
  
  /*-------------------------YOUT ARRAY: OUTPUT VARIABLES-------------------*
   |  This part contains the output variables returned by the model
   |  It is of course possible to include or delete variables.
   *------------------------------------------------------------------------*/
  
  int nyout = 0;
  
  // cardiac output variables
  yout[nyout] = Qtv; nyout++;
  yout[nyout] = Qpv; nyout++;
  yout[nyout] = Qmv; nyout++;
  yout[nyout] = Qav; nyout++;
  yout[nyout] = Qaa; nyout++;
  yout[nyout] = Qp;  nyout++;
  yout[nyout] = Qs;  nyout++;
  yout[nyout] = Qca; nyout++;
  yout[nyout] = Qc ; nyout++;
  yout[nyout] = Qcv; nyout++;
  yout[nyout] = Ppa; nyout++;
  yout[nyout] = Ppv; nyout++;
  yout[nyout] = Paorta; nyout++;
  yout[nyout] = Psa; nyout++;
  yout[nyout] = Psv; nyout++;
  yout[nyout] = Pca; nyout++;
  yout[nyout] = Pcv; nyout++;
  yout[nyout] = Plv; nyout++;
  yout[nyout] = Prv; nyout++;
  yout[nyout] = Vtot;nyout++;
  yout[nyout] = A_aorta_0; nyout++; 
  yout[nyout] = A_aorta;   nyout++; 
  yout[nyout] = Aortic_strain; nyout++;
  yout[nyout] = delta_strain;  nyout++;
  yout[nyout] = f_Baro; nyout++;
  yout[nyout] = HR;     nyout++;
  yout[nyout] = P_ra;   nyout++;
  
  // renal output variables
  // yout[nyout] = alpha_map; nyout++;
  // yout[nyout] = alpha_rap; nyout++;
  // yout[nyout] = rsna; nyout++;
  // yout[nyout] = beta_rsna; nyout++;
  // yout[nyout] = R_aa; nyout++;
  // yout[nyout] = R_ea; nyout++;
  // yout[nyout] = R_r; nyout++;
  // yout[nyout] = Phi_rb; nyout++;
  // yout[nyout] = P_gh; nyout++;
  // yout[nyout] = P_f; nyout++;
  // yout[nyout] = Phi_gfilt; nyout++;
  // yout[nyout] = Sigma_tgf; nyout++;
  // yout[nyout] = Phi_filsod; nyout++;
  // yout[nyout] = gamma_filsod; nyout++;
  // yout[nyout] = gamma_at; nyout++;
  // yout[nyout] = gamma_rsna; nyout++;
  // yout[nyout] = eta_ptsodreab; nyout++;
  // yout[nyout] = Phi_ptsodreab; nyout++;
  // yout[nyout] = Phi_mdsod; nyout++; 
  // yout[nyout] = Psi_al; nyout++;
  // yout[nyout] = eta_dtsodreab; nyout++;
  // yout[nyout] = Phi_dtsodreab; nyout++;
  // yout[nyout] = Phi_dtsod; nyout++;
  // yout[nyout] = lambda_dt; nyout++;
  // yout[nyout] = lambda_anp; nyout++;
  // yout[nyout] = eta_cdsodreab; nyout++; 
  // yout[nyout] = Phi_cdsodreab; nyout++;
  // yout[nyout] = Phi_usod; nyout++; 
  // yout[nyout] = mu_al; nyout++;
  // yout[nyout] = mu_adh; nyout++;
  // yout[nyout] = Phi_twreab; nyout++;
  // yout[nyout] = Phi_u; nyout++; 
  // yout[nyout] = xi_ksod; nyout++;
  // yout[nyout] = xi_map; nyout++;
  // yout[nyout] = xi_at; nyout++;
  // yout[nyout] = N_als; nyout++;
  // yout[nyout] = V_ecf; nyout++;
  // yout[nyout] = C_al; nyout++;
  // yout[nyout] = C_sod; nyout++;
  // yout[nyout] = Phi_sodin; nyout++;
  // yout[nyout] = N_adhs; nyout++; 
  // yout[nyout] = Phi_win; nyout++;
  // yout[nyout] = C_adh; nyout++;
  
  // respiratory and cerebral output variables
  yout[nyout] = Pb_0;   nyout++;
  yout[nyout] = Pa_0;   nyout++;
  yout[nyout] = Rc_pB_CO2; nyout++;
  //yout[nyout] = Rs_pam; nyout++;
  yout[nyout] = Rs_pa_CO2; nyout++;
  yout[nyout] = Rs;     nyout++;
  yout[nyout] = Rbrain; nyout++;
  yout[nyout] = Raorta; nyout++;
  yout[nyout] = Rtot;   nyout++;
  yout[nyout] = Csa;    nyout++;
  yout[nyout] = Csv;    nyout++;
  yout[nyout] = ca_CO2; nyout++;
  yout[nyout] = ca_O2 ; nyout++;
  yout[nyout] = cv_CO2; nyout++;
  yout[nyout] = cv_O2 ; nyout++;
  yout[nyout] = ConcSystemicCO2; nyout++;
  yout[nyout] = ConcSystemicO2 ; nyout++;
  yout[nyout] = ConcCerebralCO2; nyout++;
  yout[nyout] = ConcCerebralO2 ; nyout++;
  yout[nyout] = VolAlveolar;     nyout++;
  yout[nyout] = d_VolAlveolar;   nyout++;
  yout[nyout] = Rav;    nyout++;
  yout[nyout] = Rmv;    nyout++;
  yout[nyout] = Rpv;    nyout++;
  yout[nyout] = Rtv;    nyout++;
  yout[nyout] = Raa;    nyout++;
    
}
