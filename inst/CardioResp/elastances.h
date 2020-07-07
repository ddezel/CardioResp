/*******************************************************************
 * NAME : double TimeVaryingElastance_R/L
 *
 * DESCRIPTION: returns Plv/Prv = pressure of left/right ventricle
 *
 * INPUTS: 
 *         Theta = fraction of cardiac cycle [sec] 
 *         Vlv/Vrv = Volume of the left/right ventricle [ml]
 *         SN_activity = Sympathetic nervous system activity, between 0 and 1 (?)
 *         
 * OUTPUTS:
 *        Plv/Prv = pressure left/right ventricle 
 * 
 */

#include "utils.h"

double TimeVaryingElastance_L(double Theta, double Vlv, double SN_activity, double SN_baseline) {
  
  // variables names
  double Elv_t, E1, E2, Eml, EMl, Vld, Plv, TM, TR, tTilde;
  
  Vld = 20;                                   // [ml] zero pressure volume 
  Eml = 0.1; //0.07; //0.0715; //0.06; //0.25;       //Plv_dia / (EDVlv - Vld); about 5mmHg for 120 mL=0.04 // minimum elastance of left ventricle -> 0.01923077 mmHg
  EMl = 4.25; //5;    //4.00; //8;         //Plv_sys / (ESV_l - Vld); about 110 mmHg for 70mL=1.6 // maximum elastance of left ventricle -> 0.4331254 mmHg                                     // heart period, 60sec/62.5bpm = 0.96
  TM = 0.3;                                   // time to maximum (systolic) elastance, adapted from Davids Renal Model
  TR = 0.15;                                  // the remaining time to relaxation, adapted from Davids Renal Model
  Elv_t = 0;
  tTilde = fmod(Theta, 1);                    // fraction of cardiac cycle (seconds) 
  // duration of a heart cycle: 0.6 - 1 sec                                              
  
  // make cosinus calculations before
  E1 = cos(pi * tTilde / TM);
  E2 = cos(pi* (tTilde - TM) / TR);
  
  if (tTilde < TM) {            // heart cycle is not further than 0.3 seconds = in systole
    Elv_t = Eml + ( (0.75 + 0.25*SN_activity/SN_baseline ) * EMl - Eml) * (1 - E1) / 2;
  }
  else 
  {
    if (tTilde < TM + TR) {     // heart cycle is not further than 0.45 seconds = in systole or early diastole
      Elv_t = Eml + ( (0.75 + 0.25*SN_activity/SN_baseline ) * EMl - Eml) * (E2 + 1) / 2;
    }
    else                        // heart cycle is further than 0.45 seconds = late diastole
    {
      Elv_t = Eml;
    }
  }
  
  /* left ventricular pressure */
  Plv = Elv_t * (Vlv - Vld);
  
  
  return (Plv);
}


double TimeVaryingElastance_R(double Theta, double Vrv, double SN_activity, double SN_baseline) {
  
  // variables names
  double Erv_t, E1, E2, Vrd, Emr, EMr, Prv, TM, TR, tTilde;
  
  Vrd = 20;                           // [ml] zero pressure volume 
  Emr = 0.1; //0.073; //0.0715; //0.074; //0.25;         //Prv_dia / (EDVrv - Vrd);     about 5mmHg for 100 mL // minimum elastance of right ventricle
  EMr = 1.125; //1.27; //1.27; //8;            //Prv_sys / (ESV_r - Vrd);  about 35 mmHg for 50mL=0.7    // maximum elastance of right ventricle
  TM = 0.3;                          // maximum (systolic) elastance, adapted from Davids Renal Model
  TR = 0.15;                         // the remaining time to relaxation, adapted from Davids Renal Model
  Erv_t = 0;
  tTilde = fmod(Theta, 1); // fraction of cardiac cycle (seconds) 
  
  // make cosinus calculations before
  E1 = cos(pi * tTilde / TM);
  E2 = cos(pi * (tTilde - TM) / TR);
  
  if (tTilde < TM) {
    Erv_t = Emr + ( (0.75 + 0.25*SN_activity/SN_baseline) * EMr - Emr) * (1 - E1) / 2;
  }
  else 
  {
    if (tTilde < TM + TR) 
    {
      Erv_t = Emr + ( (0.75 + 0.25*SN_activity/SN_baseline) * EMr - Emr) * (E2 + 1) / 2;
    }
    else 
    {
      Erv_t = Emr;
    }
  }
  
  /* right ventricular pressure */
  Prv = Erv_t * max(Vrv-Vrd,0);
  return (Prv);
}