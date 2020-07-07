/*------------------------Equations for the Heart-------------------------*
 |  Below are the equations from the modified Ellwein 2016 model
 |  This model is not officially published but is promising and well 
 |  documented.
 *------------------------------------------------------------------------*/
//#include "parameters.h"
//#include "utils.h"
//#include "elastances.h"
//
//
//// resistances
//struct Resistances_Outputs {
//  double Rp, Rs, Rca, Rc, Rcv, Rbrain, Rtot;
//};
//
//// Function to return multiple values using struct
//struct Resistances_Outputs initialize()
//{
//  struct Resistances_Outputs output = {
//    
//  };
//  
//  return output;
//}
//
//
//void heart_module(double SN_activity, double AT2_activity, double Vpa, 
//                  double Vpv, double Vaorta, double Vsa, double Vsv,
//                  double Vca, double Vcv, double Ccycle_fraction,
//                  double Vlv, double Vrv, double Mean_aortic_strain,
//                  double Baro_on) {
//  
//}


//#include <stdio.h>
//
//
//// create the patient wrapper
//typedef struct
//{
//  // heart
//  struct heart
//  {
//    struct resistances
//    {
//      double Rs;
//      double Ra;
//    } res;
//    
//    struct compliances
//    {
//      double Cs;
//      double Ca;
//    } compli;
//    
//    double HR;
//    
//  } heart;
//  
//  // Brain
//  struct brain
//  {
//    int weight;
//  } brain;
//  
//} body;
//
//
//// Function to init the patient fields
//body initialize (int a, int b)
//{
//  body init = {
//    // heart
//    {
//      {a, b},
//      {1.24, 2.65},
//      59.5},
//      // brain
//      4
//  };
//  
//  return init;
//};
//
//
//int
//  main ()
//  {
//    body patient_1 = initialize(4, 5);
//    
//    double HR;
//    HR = patient_1.heart.HR;
//    
//    printf ("Heart rate (BPM): %lf", HR);
//    
//    return 0;
//  }//