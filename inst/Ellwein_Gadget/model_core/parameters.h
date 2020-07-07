union parameters {
  struct {double context,       
          // cardiovascular parameters
          //H0,              
          //H1,            
          //Rmv_open,        
          //Rav_open,        
          //Rtv_open,        
          //Rpv_open,        
          //Rmv_closed,      
          //Rav_closed,      
          //Rtv_closed,      
          //Rpv_closed,     
          
          
          // Parameters rates, tissue volumes, and gas dissociation constants
          M_CO2,           
          M_O2,           
          MB_CO2,          
          MB_O2,           
          MS_CO2,         
          MS_O2,           
          VT_CO2,          
          VT_O2,           
          VB_CO2,          
          VB_O2,           
          VS_CO2,          
          VS_O2,           
          VA_CO2,          
          VA_O2,   // 15        
          VD,              
          K1,              
          K2,              
          KCO2,            
          kCO2,           
          qp,              
          pv_CO2,          
          pv_O2,           
          pi_O2,           
          pi_CO2,          
          
          // compliances
          //Cpa,             
          //Cpv,             
          //Csa_0,           
          //Csv_0,           
          //Cca,             
          //Ccv,             
          
          tau_CO,          
          //Caorta,          
          //Raorta,          
          d_aorta_0,       
          L_aorta,         
          tau_strain,  // 40    
          f_0,             
          delta_strain_0,  
          a,               
          b,               
          f_SN,
          
          tau_autoreg,
          CO_0,
          CO_1,
          
          // kidney parameters
          //N_rsna,          
          //R_aass,          
          //P_B,             
          //P_go,           
          //C_gcf,           
          //n_etapt,         
          //n_epsilondt,     
          //n_etacd,         

          P1,              
          P2,              
          tau_renin,       
          tau_at,          
          tau_al,          
          tau_adh,         
          tau_map,         
          
          alpha1,          
          alpha2,          
          //alpha3,          
          //alpha4,        
          //Rs_0
          
          Vtid,
          resp2CC_ratio; };             
          //Vdead,           
          //alpha,           
          //beta,            
          //gamma,           
          //delta,           
          //epsilon,        
          //p_CO2tr;
  double value[48];
} p;