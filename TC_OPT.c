/* Co-based SA, V1.1 12-20-2014

  < SRO model is impelmented >
  < First kinetic model is implemented >

  Shengyen.li@nist.gov

  *** ^ ^ ***
  Inputs:
          chem_composition: nominal chemical composition of the alloy (mole fraction)
          temperature: temperature of the heat treatment (Kelvin)

  Intermidiate variables:

    Equilibrium Chemical Composition:
          chem_comp_gamma: chemical composition of gamma (mole fraction)
          chem_comp_gp: chemical composition of gamma_prime phase (mole fraction)

    Transient composition:
          chem_comp_trans: chemical composition of gamma during g-g' transformation

  Outputs:
          (1)VF: equilibrium volume fraction of gamma_prime
          (2)gp_mean_radius: optimum average radius of gamma_prime
          (3)density_all: density of the alloy
          (4)SS_stress+HP_stress: the "stress" form Peierls-Nabarro, solid solution, and Hell-Petch effect
          (5)Prec_stress: the "stress" increase from gamma_prime precipitates

  ** the output stresses are NOT shear stress

  ********** Adjust Model Parameter ************
  (1) no_ini_nu_site
  (2) E_intfac
  (3) dtime

*/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "tcapi.h"
#include "tcutils.h"

double *API_TC(int generation,int no_samples,int start_samples,int no_variables,double *XT,double M_para,int no_outputs)
{
  double *results;

// ++++ Variables for TCAPI ++++
  char *conditions,*conditions2;
  char error[80];
  char *log_file_directory;
  char *tc_installation_directory;

  TC_INT ierr;
  double value;

// ++++ Variables for process ++++ 
  int i, j, ii, jj, ierror;
  int sam, ind_XT;
  const int wstep=500,winterval=1E7,no_results=no_samples*no_outputs;
  int count_results;
  const double ST_error=0.00001;

// ++++ Variables for phsical quantities ++++ 
  int no_elements, no_elem_input;
  char *chem_elements[TC_MAX_NR_OF_ELEMENTS];    // max=40
  double chem_composition[TC_MAX_NR_OF_ELEMENTS];
  double tot_dg;
// Two temperatures: (1) for heat treatment (2) service temperature
  double temperature;
  double T_service;

  double VF[10], Vf_liq, density_all;
  int ind_L12, no_ordered;                    // phase index, high volume fraction for FCC_L12 and number of ordered phases
  double HVf;
//
// To avoid the misability gap in FCC_L12 phase:
// (1) I choose the high volume fraction one
// (2) the sample is rejected while the misability gap is shown 

  const double no_SL=2;
  double site_f_Ni[3],site_f_Al[3];
  double site_f_dis_Ni[3],site_f_dis_Al[3];
  double site_f_ord_Ni[3];
  double site_f_dis_solute[TC_MAX_NR_OF_ELEMENTS];
  double error_diff_Ni,error_diff_Al;
  double tot_Vf,tot_Part;

// ############ APBE ################
  double xNI, xsolution;
  double HM_FCCL12,HM_FCCL12_DA;
  double HM_FCCA1,HM_FCCA1_DA,HM_FCCA1_55;
  double dHM_int,dHM_int_L12,dHM_int_A1;
  double HM_PURE[TC_MAX_NR_OF_ELEMENTS];
  double V_eff_13,V_1,V_2,V_3;
  double EAPB;

  const int SRO_model=0;
// 0: Crudden et al., Acta 75 (2014) 356-370
// 1: 
// 2: Ansara J Alloys and Compounds 247 1997

//---- APBE, Crudden et al., Acta 75 (2014) 356-370 ----
  const double SRO_effect=1.1;

  double HM_SRO;
  double SRO_parameter[TC_MAX_NR_OF_ELEMENTS];
  const double Z_coordinate=4.;
  double PTOT[TC_MAX_NR_OF_ELEMENTS];
  double W2_AB;
  const double Avogadro=6.02E+23;
  const double unit_k=1.38E-23;

  double xal,xcr;
  double g_LP,gp_LP,misfit_LP,temp_LP;
  double APB_gp_LP;


// ------------- Kintics of Precipitate ------------
  int time_step;
  double tot_time,tot_time_p,time_incu,dtime;

  const int max_time_step=1E8;
//  const int max_time_step=2;

  const double c_LSW=8./9.;
  double E_intfac, alpha_Eint;
  double mvolume_gp,mvolume_gamma;
  double eq_chemsq_dM;
  double mobility[TC_MAX_NR_OF_ELEMENTS];
  double K_LSW;

  double dGM_NU,dGM_EL,GM_L12,GM_A1;  // these two are the equilibrium molar Gibbs free energies
  const int model_dGMNU=2;
  // =0, my model, tangent of the free energy curve
  // =1, Bhadeshia, Bainite in Steels, 2001, pp 131
  // =2, Wagner, Homogeneous Second-Phase Precipitation, 2001

  double GM_eq;  // Gibbs free energy of the system at equilibrium at each time step
  double GM_0;  // molar Gibbs Free energy of the initial composition (at each time step) in FCC_A1
  double dGM_critical;
  double G0_A1_AL, G0_A1_NI;
  double no_ini_nu_site;   //2E+26;   //2008B: model2_5E+26/model3_2.5E+26
  double no_nu_site;
  double nu_site;
  double gp_radius0; // critical nucleation size
  double gp_mean_radius,gp_mean_radius_p;  // mean radius of the precipitates
  double gp_critical_radius, dG_GibbsThomson;  // Maximum radius of precipitates, minimum driving force
  double dG_int, dG_grow;  // two energies to evaluate the stop of the growth
  const int model_dGgrow=0;
  // =0, the driving force is calculated at the annealing temperature
  // =1, Lupis, Chemical Thermodynamics of Materials, 1983, pp89-94

//---------------------------------------------------------------------------------------------
  double Zeldfactor, beta_trans;
  double chem_poten[TC_MAX_NR_OF_ELEMENTS], thermo_factor[TC_MAX_NR_OF_ELEMENTS];
  double activity[2][TC_MAX_NR_OF_ELEMENTS];
  double diffusivity[2][TC_MAX_NR_OF_ELEMENTS], diffusivity_off[TC_MAX_NR_OF_ELEMENTS];
  double diff[TC_MAX_NR_OF_ELEMENTS];
  double Velocity_int;  // for Avrami eq
  double chem_comp_trans[TC_MAX_NR_OF_ELEMENTS];
  double chem_comp_matrixint[TC_MAX_NR_OF_ELEMENTS];
  double k_chem;

// Qing et al., Analytical treatment..., Acta, 2008 
  double omega,lambda,C_growth,func_lambda;
  double lambda_max,lambda_min;
  double error_beta_fl;

  int onnucleation,oncoarsening;
  int onnustep;
  int no_prec_steps,no_calc_prec_steps;

//  double **tot_precipitate[max_time_step], **precipitates[4][max_time_step];
  double *tot_precipitate, **precipitates= (double **)malloc(4 * sizeof(double *));
  double swap;
// information of precipitates [feature][time step]
// [1][:] density of nucleus
// [2][:] radius of precipitates
// [3][:] volume fraction of precipitates

  const int kinetics_model=1;
// Two kinetics models are implemented:
// (0) use Avrami equation to calculate the volume fraction of precipitation
//     use LSW model to calculate the size of the precipitates
// (1) use nucleation model to get the density of nuclei
//     use LSW model to calculate the size of the precipitates
//     the volume fraction of precipitation is (density of nuclei*precipitate volume)

  const int model_velocity=3;
// =0, the same as Panda
// =1, Qing, Analytical treatment of diffusion..., 2008
// =2, M&P with diffusion distance
// =3, Rougier, Numerical simulation of ..., 2013

// ************* YS **************
  double YS,YS_p;
  double chem_comp_gamma[TC_MAX_NR_OF_ELEMENTS],chem_comp_gp[TC_MAX_NR_OF_ELEMENTS];
  double Vf_gp, Vf_gp_eq;    // mean characters of the microstructure
  const double grain_size=20E-6;
  double weak_stress,strong_stress,strong_stress_p;
  double Prec_stress,dislocation_stress,Back_stress;
  double Orwan_stress,Orwan_stress_p;
  double LT,burgers;

// ++++ Pollock et al, Creep resistance of ..., Acta Metal, 40, 1, 1992 ++++
// const double shear_modulus=4.82*pow(10,10);
// const double burgers=2.54*pow(10,-10);
// ++++ RC Reed, Superalloy - Fundamentals and applications, 2006 (Chap 2) ++++
// ++ To calculate shear stress ++
//  double shear_modulus=8*pow(10,10);
//  const double burgers=2.5*pow(10,-10);
// ++++ RC Reed, Superalloy - Fundamentals and applications, 2006 (Chap 3) ++++
// ++ To calculate Orowan Stress ++
//  double shear_modulus=5*pow(10,10); 
//  const double burgers=2.5*pow(10,-10);
// ++++ Thomas et al., J Mat. Pro. Tech, 177, 2006, 469 ++++
// E = 2 ( 1 + 1/3 ) Mu, poission ration=1/3
  double shear_modulus;
  double Taylor_factor=3.06;

  const double pi=acos(-1.),R_gas=8.314,KB=1.3806488E-23;
  const double C1=1.91, C2=0.22, cw=1.;

// ++++ Ahmadi et al, Yield strength ..., Materials science & Engineering A, 608, 114, 2014 ++++
  const double YS_0=21.8;
  const double C_HPeffect=0.158;   // MPa*m^(0.5)
  double HP_stress;
  double YS_SS[TC_MAX_NR_OF_ELEMENTS], SS_stress;
  double YieldStress;  

  double temp;
 
//  Two models are implemented for calculating maximum yield stress
//    0) iterative calculation to get the cross stress of weak pair and strong pair dislocation
//    1) derivative of the strong pair dislocation model (Crudden et al, Acta Mat, 2014)

  const int max_stress_model=0;


// ////////// Model parameters N0 //////////
  no_ini_nu_site=M_para;


//printf("-----------------%d %d %d \n",generation,no_variables,no_outputs);

//printf("%d, %d, %d, %d; %E; %E, %d \n",generation,no_samples,start_samples,no_variables,XT[0],M_para,no_outputs);

  results=calloc(no_results,sizeof(double));
  for(i=0;i<no_results;i++){
    results[i]=0.;
  }

  tot_precipitate=calloc(max_time_step,sizeof(double));


  for(i=0;i<4;i++){
    precipitates[i]=(double *)malloc(max_time_step*sizeof(double));
  }

//  printf("%E \n",precipitates[3][1]);




// +++ initialize the physical system +++
//  no_elements=no_variables-1;  // +(1)Ni -[ (1)temperature for heat treatment + (2)service temperature ]
  no_elements=no_variables;
  no_elem_input=1;       //the start of the inputs composition: must >= 1
  chem_elements[0]="NI";
  chem_elements[1]="AL";
  chem_elements[2]="CR";

// ++++ Ahmadi et al, Yield strength ..., Materials science & Engineering A, 608, 114, 2014 ++++
// AL:225; CO:39.4; CR:337; NB:1183; C:1061; FE:153; MO:1015; TI:775; W:997
// ++++ RC Reed, Superalloy - Fundamentals and applications, 2006 (Chap 2) ++++
// AL:225; CO:39; CR:337; FE:153; MO:1015; TI:775; W:997; ZR:2359; HF:1401; MN:448
// ++++ Roth et al, Modeling solid solution ..., Metal. & Materials Tran. A, 28A, 1329, 1997 ++++
// AL:225; SI:275; ZN:386; GA:310; GE:332; IN:985; SN:1225; SB:960; TI:775; V:408; ZR:2359; 
// HF:1401; NB:1183; TA:1191; CR:337; MO:1015; W:977; MN:448; FE:153; RU:1068; CO:39.4; RH:520; 
// CU:86.7; C:1061; PD:492
//
// Unit: MPa per (at fraction)^(-0.5)
  YS_SS[1]=225.0;
  YS_SS[2]=337.0;

  tc_deinit();
  log_file_directory = ".";
  tc_installation_directory = ".";



  for(sam=start_samples;sam<no_samples;sam++){
    count_results=sam*no_outputs;

    printf("========== SAMPLE:  %i ============\n",sam);
// ***** Calculate Enthalpy of FCC_L12 *****
// ++++ RE-Initialize the system ++++
//    tc_deinit();
//    log_file_directory = ".";
//    tc_installation_directory = ".";
    tc_init_root3(log_file_directory, tc_installation_directory);
    tc_open_database("TCNI6");
    for(i=0;i<no_elements;i++){
      tc_element_select(chem_elements[i]);
    }

    tc_phase_reject("*");
    tc_phase_select("FCC_A1");
    tc_phase_select("FCC_L12");
    tc_phase_select("LIQ");

    tc_get_data();

// ++++ Setup invariant conditions ++++
    tc_set_condition("N",1);
    tc_set_condition("P",101325);

    ind_XT=sam*no_variables;

    for(i=no_elem_input;i<no_elements;i++){
      chem_composition[i]=XT[ind_XT];
      ind_XT += 1;
    }
    temperature=XT[ind_XT];
    ind_XT += 1;
    tc_set_condition("T",temperature);

    T_service=1123.;  //XT[ind_XT];
//    ind_XT += 1;

//    printf("inputs: %E, %E, %E, %E \n",chem_composition[1],chem_composition[2],temperature,T_service);

    shear_modulus=(298.*3./8.)*(1.-0.5*(T_service-300.)/1673.)*1E+9;


// the excess enthalpy for Short rang ordering
// NIST database
    SRO_parameter[1]=6*(-13415.515+2.0819247*T_service);
    SRO_parameter[2]=6*(-4300.);

    for(i=no_elem_input;i<no_elements;i++)
    {
      conditions=NULL;
      asprintf(&conditions,"%s%s%s","X(",chem_elements[i],")");

      tc_set_condition(conditions,chem_composition[i]);
      printf("%s:  %f\n",conditions,chem_composition[i]);
    }

    ierr=0;
    if ( tc_error(&ierr,error,sizeof(error)) ) {
      printf("ERROR in setup conditions");
//      fprintf(stdout,"error: %d : %s\n",ierr,error);
      fprintf(stdout,"error: %s\n",error);
      tc_reset_error();
      continue;
    }


// ++++ Compute the equilibrium ++++
    tc_compute_equilibrium();

    ierr=0;
    if ( tc_error(&ierr,error,sizeof(error)) ) {
      fprintf(stdout,"error while calculating FCC_L12 equilibrium: %s\n",error);
      tc_reset_error();
      for(ierror=0;ierror<no_outputs;ierror++){
        results[count_results+ierror]=0.;
      }
      continue;
    }

    Vf_liq=tc_get_value("VPV(LIQ)");
    if(Vf_liq>0.0){
      printf("Higher than solidus temperature: %f, %f \n",temperature,Vf_liq);
      for(ierror=0;ierror<no_outputs;ierror++){
        results[count_results+ierror]=0.;
      }
      continue;
    }

// ---- Retrieve the desired value ----
    tot_Vf=0.0;
    HVf=0.0;
    no_ordered=0;
    ind_L12=-1;

    for(i=1;i<=5;i++){

      VF[i]=0.0;
      asprintf(&conditions,"%s%i%s","VPV(FCC_L12#",i,")");
      VF[i]=tc_get_value(conditions);

//      printf("%s:  %f \n",conditions,VF[i]);

      ierr=0;
      if ( tc_error(&ierr,error,sizeof(error)) ) {
        printf("ERROR while getting results \n");
        fprintf(stdout,"error: %s\n",error);
        tc_reset_error();
        continue;
      }

      for(j=1;j<=no_SL;j++){
        site_f_Ni[j]=0.0;
        asprintf(&conditions,"%s%i%s%i%s","y(FCC_L12#",i,",NI#",j,")");
        site_f_Ni[j]=tc_get_value(conditions);

        site_f_Al[j]=0.0;
        asprintf(&conditions,"%s%i%s%i%s","y(FCC_L12#",i,",AL#",j,")");
        site_f_Al[j]=tc_get_value(conditions);
      }

      error_diff_Ni = 0.;
      error_diff_Al = 0.;

      error_diff_Ni = fabs(site_f_Ni[1]-site_f_Ni[2]);
      error_diff_Al = fabs(site_f_Al[1]-site_f_Al[2]);


// ------ decide the phase index of order phase, FCC_L12 (and?) -----------
      if( (error_diff_Ni > ST_error) && (error_diff_Al > ST_error) && VF[i] > 0. ){
        no_ordered+=1;

        if(VF[i] > HVf){
          ind_L12 = i;
          HVf=VF[i];

//          PTOT=0.;
//          PTOT=site_f_Ni[1]*(1.-site_f_Ni[1])*site_f_Ni[2]*(1.-site_f_Ni[2]);
//          printf("PTOT: %f, %f, %f, %f \n",site_f_Ni[1],(1.-site_f_Ni[1]),site_f_Ni[2],1.-site_f_Ni[2]);

          for(j=0;j<no_elements;j++){
            chem_comp_gp[j]=0.0;
            asprintf(&conditions,"%s%i%s%s%s","x(FCC_L12#",i,",",chem_elements[j],")");
            chem_comp_gp[j]=tc_get_value(conditions);
//            printf("Composition of L12: %d, %E \n",j,chem_comp_gp[j]);

            activity[1][j]=0.;
            asprintf(&conditions,"%s%s%s%i%s","ACR(",chem_elements[j],",FCC_L12#",i,")");
            activity[1][j]=tc_get_value(conditions);
          }

          mvolume_gp=0.0;
          asprintf(&conditions,"%s%i%s","VM(FCC_L12#",i,")");
          mvolume_gp=tc_get_value(conditions);

          gp_LP=pow(mvolume_gp/(Avogadro/4.),1./3.);
//          printf("Molar Volume of FCC_L12, m^3: %E, %E \n",mvolume_gp,gp_LP);

          site_f_ord_Ni[1]=site_f_Ni[1];
          site_f_ord_Ni[2]=site_f_Ni[2];

          GM_L12=0.0;
          asprintf(&conditions,"%s%i%s","GM(FCC_L12#",i,")");
          GM_L12=tc_get_value(conditions);
//          printf("Molar Energy of L12: %E \n",GM_L12);

          dHM_int_L12=0.;
          asprintf(&conditions,"%s%i%s","HM(FCC_L12#",i,")");
          dHM_int_L12=tc_get_value(conditions);
        }
      }
      else if(i==1){
        for(j=0;j<no_elements;j++){
          chem_comp_gamma[j]=0.0;
          asprintf(&conditions,"%s%i%s%s%s","x(FCC_L12#",i,",",chem_elements[j],")");
          chem_comp_gamma[j]=tc_get_value(conditions);
//          printf("Composition of A1: %d, %E \n",j,chem_comp_gamma[j]);
        }

        GM_A1=0.0;
        asprintf(&conditions,"%s%i%s","GM(FCC_L12#",i,")");
        GM_A1=tc_get_value(conditions);
        dHM_int_A1=tc_get_value("HM");

        mvolume_gamma=0.0;
        asprintf(&conditions,"%s%i%s","VM(FCC_L12#",i,")");
        mvolume_gamma=tc_get_value(conditions);

        g_LP=pow(mvolume_gamma/(Avogadro/4.),1./3.);
//        printf("Molar Volume of FCC_A1, m^3: %E, %E \n",mvolume_gamma,g_LP);

      }

      tot_Vf+=VF[i];
//      printf("total Volume fraction, index of FCC_L12: %f, %i\n",tot_Vf,ind_L12);


      if(tot_Vf>=0.9999999){
//        printf("total Volume fraction, index of FCC_L12: %f, %i, %i\n",tot_Vf,ind_L12,no_ordered);
        break;
      }

    }

// ---- Checking points of the microstructure ----
    if(tot_Vf<0.9999999 || VF[ind_L12]>0.99){
      printf("Not in gamma/gamma_prime, two phase field: %E, %E \n",tot_Vf,VF[ind_L12]);
      for(ierror=0;ierror<no_outputs;ierror++){
        results[count_results+ierror]=0.;
      }
      continue;
    }

    printf("Ordered number %d \n",ind_L12);
    if(no_ordered > 1){
      printf("There are more than ONE order phase: %d \n",no_ordered);
      for(ierror=0;ierror<no_outputs;ierror++){
        results[count_results+ierror]=0.;
      }
      continue;
    }

    if(ind_L12==-1){
      printf("**** FCC_L12 does not exist **** \n");
      for(ierror=0;ierror<no_outputs;ierror++){
        results[count_results+ierror]=0.;
      }
      continue;
    }

// **************************** for Ni 
    if((VF[ind_L12]<0.4) || (VF[ind_L12]>0.95)){
      printf("Vf_L12 is too small: %E, %E \n",tot_Vf,VF[ind_L12]);
      for(ierror=0;ierror<no_outputs;ierror++){
        results[count_results+ierror]=0.;
      }
      continue;
    }


    printf("Equilibrium Vf_L12: %E \n",VF[ind_L12]);
// ---------------------------------------------------------------------
    tc_set_condition("T",T_service);
    tc_compute_equilibrium();

// ++++ Density of the alloy ++++
    density_all=0.0;
    conditions="BV";
    density_all=tc_get_value(conditions)/1E+6;
//    printf("Density,   g/cm^3: %f\n",density_all);


// ------------- dH_L12 -------------
    tc_init_root3(log_file_directory, tc_installation_directory);
    tc_open_database("TCNI6");

    for(i=0;i<no_elements;i++){
      tc_element_select(chem_elements[i]);
    }

    tc_phase_reject("*");
    tc_phase_select("FCC_L12");
    tc_get_data();

    for(i=no_elem_input;i<no_elements;i++)
    {
      conditions=NULL;
      asprintf(&conditions,"%s%s%s","X(",chem_elements[i],")");
      tc_set_condition(conditions,chem_comp_gp[i]);
    }

    tc_set_condition("N",1.);
    tc_set_condition("P",101325);
    tc_set_condition("T",T_service);

    tc_compute_equilibrium();

//    asprintf(&conditions,"%s%i%s","HM(FCC_L12#",ind_L12,")");
//    HM_FCCL12_DA=tc_get_value(conditions);
    APB_gp_LP=pow(tc_get_value("VM")/(Avogadro/4.),1./3.);

    HM_FCCL12_DA=tc_get_value("HM");

//  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ The model to calculate SRO Energy ++++++
// +++++++++ Calculating Enthalpy of FCC_A1 +++++++++++++++
//    tc_deinit();
//    log_file_directory = ".";
//    tc_installation_directory = ".";
    tc_init_root3(log_file_directory, tc_installation_directory);
    tc_open_database("TCNI6");

    for(i=0;i<no_elements;i++){
      tc_element_select(chem_elements[i]);
    }

    tc_phase_reject("*");
    tc_phase_select("FCC_A1");

    tc_get_data();

// ----- Setup diffusion database -----
    tc_append_database("NIMOB");

    for(i=0;i<no_elements;i++){
      tc_element_select(chem_elements[i]);
    }

    tc_phase_reject("*");
    tc_phase_select("FCC_A1");
    tc_get_data();


// ++++ Setup invariant conditions ++++ 
    tc_set_condition("N",1.);
    tc_set_condition("P",101325);
    tc_set_condition("T",T_service);

    for(i=no_elem_input;i<no_elements;i++)
    {
      conditions=NULL;
      asprintf(&conditions,"%s%s%s","X(",chem_elements[i],")");
      tc_set_condition(conditions,chem_comp_gp[i]);
    }

// ++++ Compute the equilibrium ++++
    tc_compute_equilibrium();

    ierr=0;
    if ( tc_error(&ierr,error,sizeof(error)) ) {
      fprintf(stdout,"error while calculating FCC_A1 equilibrium: %s\n",error);
      tc_reset_error();
      for(ierror=0;ierror<no_outputs;ierror++){
        results[count_results+ierror]=0.;
      }
      continue;
    }

    HM_FCCA1_DA=tc_get_value("HM");
    HM_FCCL12_DA-=HM_FCCA1_DA;

//    printf("HM:  %f,  %f \n",HM_FCCA1_DA,HM_FCCL12_DA);


// ========================================================================= SRO Model ========
//  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ (1) Miodownik 1995 $$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    if (SRO_model == 0){

//printf("sum:::: %E \n",HM_FCCA1_DA+HM_FCCL12_DA);
      HM_FCCA1=HM_FCCA1_DA/SRO_effect;
      HM_FCCL12=HM_FCCL12_DA+(1.-1./SRO_effect)*HM_FCCA1_DA;
//printf("sum:::: %E \n",HM_FCCA1+HM_FCCL12);
//      printf("H:  %f, %f, %f, %f \n",HM_FCCA1_DA,HM_FCCA1,HM_FCCL12_DA,HM_FCCL12);
    }
//  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ (2) Sundman CALPHAD 22 1998 $$$$$$$$$$$$$$$$$$$$
    else if (SRO_model == 1){

      asprintf(&conditions,"%s","y(FCC_A1#1,NI#1)");
      site_f_dis_Ni[1]=tc_get_value(conditions);

      for(i=1;i<no_elements;i++){
        PTOT[i]=0.;
        asprintf(&conditions,"%s%s%s","y(FCC_A#1,",chem_elements[i],"#1)");
        site_f_dis_solute[i]=tc_get_value(conditions);

        PTOT[i]=site_f_dis_Ni[1]*site_f_dis_solute[i]*site_f_dis_Ni[1]*site_f_dis_solute[i];
//        printf("PTOT: %f, %f \n",PTOT[i],site_f_dis_Ni[1]);
      }

      tc_set_condition("N",1.);
      tc_set_condition("P",101325);
      tc_set_condition("T",T_service);

      for(i=no_elem_input;i<no_elements;i++)
      {
        conditions=NULL;
        asprintf(&conditions,"%s%s%s","X(",chem_elements[i],")");
        tc_set_condition(conditions,0.0);
      }

      ierr=0;
      if ( tc_error(&ierr,error,sizeof(error)) ) {
        printf("ERROR in setup conditions");
        fprintf(stdout,"error: %s\n",error);
        tc_reset_error();
      }

// ++++ Compute the equilibrium ++++
      tc_compute_equilibrium();

      ierr=0;
      if ( tc_error(&ierr,error,sizeof(error)) ) {
        fprintf(stdout,"error while calculating FCC_A1 equilibrium: %s\n",error);
        tc_reset_error();
        for(ierror=0;ierror<no_outputs;ierror++){
          results[count_results+ierror]=0.;
        }
        continue;
      }

// --- pure element of NI ---
      HM_PURE[0] = tc_get_value("HM(FCC_A1)");

      tc_set_condition("X(AL)",1.0);
      ierr=0;
      if ( tc_error(&ierr,error,sizeof(error)) ) {
        printf("ERROR in setup conditions");
        fprintf(stdout,"error: %s\n",error);
        tc_reset_error();
     }

// ++++ Compute the equilibrium ++++
      tc_compute_equilibrium();

      ierr=0;
      if ( tc_error(&ierr,error,sizeof(error)) ) {
        fprintf(stdout,"error while calculating FCC_A1 equilibrium: %s\n",error);
        tc_reset_error();
        for(ierror=0;ierror<no_outputs;ierror++){
          results[count_results+ierror]=0.;
        }
        continue;
      }    

      HM_PURE[1] = tc_get_value("HM(FCC_A1)");

// --- 0.5NI-0.5AL ---
      tc_set_condition("X(AL)",0.5);
      ierr=0;
      if ( tc_error(&ierr,error,sizeof(error)) ) {
        printf("ERROR in setup conditions");
        fprintf(stdout,"error: %s\n",error);
        tc_reset_error();
     }

// ++++ Compute the equilibrium ++++
      tc_compute_equilibrium();

      ierr=0;
      if ( tc_error(&ierr,error,sizeof(error)) ) {
        fprintf(stdout,"error while calculating FCC_A1 equilibrium: %s\n",error);
        tc_reset_error();
        for(ierror=0;ierror<no_outputs;ierror++){
          results[count_results+ierror]=0.;
        }
        continue;
      }

      HM_FCCA1_55 = tc_get_value("HM(FCC_A1)");


// Sundman, CALPHAD, 22, 1998, 335  
//    W2_AB=(HM_FCCA1_DA-HM_PURE[0])*(HM_FCCA1_DA-HM_PURE[0]);
// Abe, Computer Coupling of Phase Diagrams and Tech, 27, 2003, 403.
//    W2_AB=(HM_FCCA1_DA+(1.0-chem_comp_gp[1])*HM_PURE[0]+chem_comp_gp[1]*HM_PURE[1]);

// --- My model ---
      W2_AB=(HM_FCCA1_55-2.*(HM_PURE[0]+HM_PURE[1]));

//    printf("W2AB: %f, %f, %f, %f \n",W2_AB,HM_FCCA1_DA,HM_PURE[0],HM_PURE[1]);
// ()()()()()()()()()()()()()()()()()()()()()()()() //

//    Z_coordinate=4.;
//    printf("for SRO: %f, %f, %f \n",Z_coordinate,PTOT,temperature);

// Sundman, CALPHAD, 22, 1998, 335
//    HM_SRO=-(Z_coordinate*PTOT*W2_AB)/(R_gas*temperature);
// Abe, Computer Coupling of Phase Diagrams and Tech, 27, 2003, 403.

      HM_SRO=0.;
      for(i=1;i<no_elements;i++){
        HM_SRO+=PTOT[i]*W2_AB;
      }

//      printf("SRO: ***** %f, %f, %f ***** \n",HM_SRO,HM_FCCA1_DA,(HM_SRO/HM_FCCA1_DA));


      HM_FCCA1=HM_FCCA1_DA-HM_SRO;
      HM_FCCL12=HM_FCCL12_DA+HM_SRO;
    }
//  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ (3) Ansara J Alloys and Compounds 247 1997 $$$$$$$$
    else{

      asprintf(&conditions,"%s","y(FCC_A1#1,NI#1)");
      site_f_dis_Ni[1]=tc_get_value(conditions);

//      PTOT=0.;
//      PTOT=site_f_dis_Ni[1]*(1.-site_f_dis_Ni[1])*site_f_dis_Ni[1]*(1.-site_f_dis_Ni[1]);
//      PTOT=site_f_ord_Ni[1]*(1.-site_f_ord_Ni[1])*site_f_ord_Ni[2]*(1.-site_f_ord_Ni[2]);
//      printf("PTOT: %E, %f, %f \n",PTOT,site_f_ord_Ni[1],site_f_dis_Ni[1]);

      for(i=1;i<no_elements;i++){
        PTOT[i]=0.;
        asprintf(&conditions,"%s%s%s","y(FCC_A#1,",chem_elements[i],"#1)");
        site_f_dis_solute[i]=tc_get_value(conditions);
        
        PTOT[i]=site_f_dis_Ni[1]*site_f_dis_solute[i]*site_f_dis_Ni[1]*site_f_dis_solute[i];
//        printf("PTOT: %f, %f \n",PTOT[i],site_f_dis_Ni[1]);
      }

      HM_SRO=0.;
      for(i=1;i<no_elements;i++){
//        HM_SRO+=PTOT*pow((SRO_parameter[i]-G0_A1_AL-G0_A1_NI),2.)/(Z_coordinate*R_gas*temperature);
        HM_SRO+=PTOT[i]*SRO_parameter[i];
      }

      printf("SRO: ***** %f, %f, %f ***** \n",HM_SRO,HM_FCCA1_DA,(HM_SRO/HM_FCCA1_DA));

      HM_FCCA1=HM_FCCA1_DA-HM_SRO;
      HM_FCCL12=HM_FCCL12_DA+HM_SRO;
    }
//    printf("H:  %f, %f, %f, %f \n",HM_FCCA1_DA,HM_FCCA1,HM_FCCL12_DA,HM_FCCL12);

    xNI=chem_comp_gp[0];
//    printf("Nickle content in gamma prime: %f\n",xNI);
    xsolution=1.-xNI;

// ^^^^^^^^^^^^^^^^^^^^^^^^^ Collins, 2014 ^^^^^^^^^^^^^^^^^^^^^^^^^  
//    V_eff_13=-( 3.*HM_FCCA1+HM_FCCL12*(1.-xsolution)/xsolution )/( 24.*Avogadro*xsolution*(1.-xsolution) );   // Collins, 2014
//    V_2=( HM_FCCL12*(1.-xsolution) - HM_FCCA1*xsolution )/ (12.*Avogadro*pow(xsolution,2.)*(1.-xsolution) );  // Collins, 2014


// ^^^^^^^^^^^^^^^^^^^^^^^  Crudden, 2014 ^^^^^^^^^^^^^^^^^^^^^^
    V_eff_13=-( 3.*HM_FCCA1+HM_FCCL12*(1.-xsolution)/xsolution )/( 24.*R_gas*xsolution*(1.-xsolution) );
    V_2=( HM_FCCL12*(1.-xsolution) - HM_FCCA1*xsolution )/ (12.*R_gas*pow(xsolution,2.)*(1.-xsolution) );


//    printf("%E, %E, %E, %E, %E \n",HM_FCCL12,(1.-xsolution), HM_FCCA1,xsolution, Avogadro,pow(xsolution,2.));

    V_1=0.75*V_eff_13;
    V_3=0.125*V_eff_13;

    EAPB=(V_1 - 3.*V_2 + 4.*V_3)/(sqrt(3.)*pow(APB_gp_LP*1E+10,2.)*1E+3);

    if ( EAPB<0.01 ) {
      EAPB=0.01;
    }


    printf("EAPB: %E  J/m^2\n",EAPB);


// =================================================== Kinetic simulation of gamma_prime precipitation =====
// ++++ Setup invariant conditions ++++
    tc_set_condition("N",1.);
    tc_set_condition("P",101325);
    tc_set_condition("T",temperature);

// Equilibrium diffusivity
    for(i=no_elem_input;i<no_elements;i++)
    {
      conditions=NULL;
      asprintf(&conditions,"%s%s%s","X(",chem_elements[i],")");
      tc_set_condition(conditions,chem_comp_gamma[i]);
    }

    tc_compute_equilibrium();
    ierr=0;
    if ( tc_error(&ierr,error,sizeof(error)) ) {
      printf("ERROR in setup conditions");
      fprintf(stdout,"error: %s\n",error);
      tc_reset_error();
    }

    for(i=1;i<no_elements;i++){
      asprintf(&conditions,"%s%s%s%s%s%s%s","DC(FCC_A1,",chem_elements[i],",",chem_elements[i],",",chem_elements[0],")");
      diffusivity[1][i]=tc_get_value(conditions);  // Diffusivity of equilibrium composition
    }


    for(i=no_elem_input;i<no_elements;i++)
    {
      conditions=NULL;
      asprintf(&conditions,"%s%s%s","X(",chem_elements[i],")");
      tc_set_condition(conditions,chem_composition[i]);
    }

    ierr=0;
    if ( tc_error(&ierr,error,sizeof(error)) ) {
      printf("ERROR in setup conditions");
      fprintf(stdout,"error: %s\n",error);
      tc_reset_error();
    }

// ++++ Compute the equilibrium ++++
    tc_compute_equilibrium();

    ierr=0;
    if ( tc_error(&ierr,error,sizeof(error)) ) {
      results[count_results]=-1.0;
      fprintf(stdout,"error while calculating FCC_A1 equilibrium: %s\n",error);
      tc_reset_error();
      for(ierror=0;ierror<no_outputs;ierror++){
        results[count_results+ierror]=0.;
      }
      continue;
    }

    GM_0=0.0;   // total fre energy of the system before phase transformation
    asprintf(&conditions,"%s","GM");
    GM_0=tc_get_value(conditions);
//    printf("Total Free Energy: %E \n",GM_0);

    GM_eq=0.;
    i=1;
//    for(i=1;i<no_elements;i++){
      GM_eq+=(GM_A1-GM_L12)*(chem_composition[i]-chem_comp_gamma[i])/(chem_comp_gp[i]-chem_comp_gamma[i]);
//    }
    GM_eq=GM_A1-GM_eq;
    dG_grow=(GM_0-GM_eq);
    if(dG_grow<0.1){dG_grow=0.; }
//  printf("dG growth: %E \n",dG_grow,GM_0,GM_eq);

    dGM_NU=0.0;
    if(model_dGMNU==0){
      for(i=1;i<no_elements;i++){

    // retrive the chemical potential of the initial alloy
        asprintf(&conditions,"%s%s%s","mu(",chem_elements[i],")");
        chem_poten[i]=tc_get_value(conditions);
        dGM_NU+=chem_poten[i]*(chem_comp_gp[i]-chem_composition[i]);

//        printf("Chemical Potential: %E, %E, %E \n",chem_poten[i],chem_comp_gp[i],chem_composition[i]);
      }

//printf("dG NU: %E %E %E \n", dGM_NU,GM_L12,GM_0);
     dGM_NU+=(GM_0-GM_L12);  // This is driving force, NOT dG ==> This value is positive
//   printf("dG NU: %E %E %E \n", dGM_NU,GM_0,GM_L12);
    }
    else if(model_dGMNU==1){
      dGM_NU=dG_grow*(chem_comp_gp[1]-chem_comp_gamma[1])/(chem_comp_gp[1]-chem_composition[1]);
    }
    else{

      chem_composition[0]=1.;
      for(i=0;i<no_elements;i++){
        activity[0][i]=0.;
        asprintf(&conditions,"%s%s%s","ACR(",chem_elements[i],",FCC_A1)");
        activity[0][i]=tc_get_value(conditions);

        dGM_NU+=chem_comp_gp[i]*log(activity[0][i]/activity[1][i]);
        chem_composition[0]-=chem_composition[1];
      }
      dGM_NU=8.314*temperature*dGM_NU; //J/mol  //mvolume_gp;

    }
//    printf("dG Nucleation: %E %E \n",dGM_NU,dGM_NU/mvolume_gp );   

// Wagner, Homogeneous Second-Phase Precipitation, 200
    dGM_EL=0.;
    misfit_LP=2.*(g_LP-gp_LP)/(g_LP+gp_LP);
    dGM_EL=2.*shear_modulus*(1.3/0.7)*misfit_LP*misfit_LP;   // poission ratio is taken as 0.3
    dGM_EL=-dGM_EL*mvolume_gp;   //2.67e+6*mvolume_gp;
//    printf("dG Elastic: %E %E %E %E %E \n",dGM_EL,dGM_EL/mvolume_gp,g_LP,gp_LP,misfit_LP );

    if ( dGM_NU+dGM_EL <=0.0 ) {
      fprintf(stdout,"The chemical driving force for nucleation is insufficient: %f \n", dGM_NU);
      for(ierror=0;ierror<no_outputs;ierror++){
        results[count_results+ierror]=0.;
      }
      continue;
    }

    for(i=1;i<no_elements;i++){

      asprintf(&conditions,"%s%s%s","mu(",chem_elements[i],")");
      chem_poten[i]=tc_get_value(conditions);
//      printf("Chemical Potential: %E, %E, %E \n",chem_poten[i],chem_comp_gp[i],chem_composition[i]);

      asprintf(&conditions,"%s%s%s","M(FCC_A1,",chem_elements[i],")");
      mobility[i]=tc_get_value(conditions);

      asprintf(&conditions,"%s%s%s%s%s%s%s","DC(FCC_A1,",chem_elements[i],",",chem_elements[i],",",chem_elements[0],")");
      diffusivity[0][i]=tc_get_value(conditions);

//      printf("Mobility & Diffusivity: %E, %E \n",mobility[i],diffusivity[i]);
    }
//    printf("Shear Modulus: %E \n",shear_modulus);

    burgers=g_LP/sqrt(2.);
    LT=shear_modulus*burgers*burgers/2.;
//    printf("Burgers and line tension: %E %E %E \n",burgers,LT,shear_modulus);


// Li et al, Met Mat Tran A, 2002, 33A, 3367
// For calculating EAPB, HM_FCCL12_DA -= HM_FCCA1_DA;
//    alpha_Eint=0.022/fabs(dHM_int);

//    E_intfac=0.020;    //alpha_Eint*fabs(dHM_int);
//    E_intfac=0.018;
//    dHM_int_A1=tc_get_value("HM");

    dHM_int=0.0;
    dHM_int=dHM_int_L12-dHM_int_A1;
//    alpha_Eint=E_intfac/fabs(dHM_int);

    E_intfac=(3.7475E-3*temperature-2.226)*1E-6*fabs(dHM_int);
    printf("Interfacial Energy: %E, %E, %E \n",E_intfac,fabs(dHM_int),alpha_Eint);

//    if(E_intfac>0.026 || (E_intfac!=E_intfac)){
    if((E_intfac!=E_intfac)){
      fprintf(stdout,"High/Error Interfacial Energy: %E \n", E_intfac);
      for(ierror=0;ierror<no_outputs;ierror++){
        results[count_results+ierror]=0.;
      }
      continue;
    }

    k_chem=2*E_intfac*mvolume_gp/(R_gas*temperature);

    Vf_gp_eq=VF[ind_L12];
    Vf_gp=0.0;
//    printf("Total volume fraction of L12: %E \n", Vf_gp_eq);

    dGM_critical=16.*pi*pow(E_intfac,3)/(3.*pow(((dGM_NU+dGM_EL)/mvolume_gp),2));
// Svoboda et al., 2004 & TC-Presima
//    printf("Chemical Driving Force: %E \n", dGM_critical);


    gp_radius0=(2.*E_intfac*mvolume_gp/(dGM_NU+dGM_EL) );   // critical nucleation size
//    printf("Critical Nucleation Size, meter: %E \n",gp_radius0);

// ---- Svoboda et al., Mat Sci Eng A, 385, 2004, 166 ----
// ---- from TC_presima ----
    Zeldfactor=pow((dGM_NU+dGM_EL),2)/(8.*pi*Avogadro*pow(E_intfac,2)*mvolume_gp);
    Zeldfactor=Zeldfactor*sqrt(E_intfac/(KB*temperature));
//    printf("Zfactor:  %E \n", Zeldfactor);

    beta_trans=0;
    eq_chemsq_dM=0.;

//    i=1;
    for(i=1;i<no_elements;i++){

// **** Qing ****
      error_beta_fl=10000.;

      chem_comp_matrixint[i]=chem_comp_gamma[i]*exp(k_chem/gp_radius0);
      printf("COMPOSITION %d: %E, %E %E \n",i,chem_comp_matrixint[i],chem_comp_gamma[i],chem_composition[i]);

      omega=(chem_composition[i]-chem_comp_gamma[i])/(chem_comp_gp[i]-chem_comp_gamma[i]);
//      omega=(chem_composition[i]-chem_comp_matrixint[i])/(chem_comp_gp[i]-chem_comp_matrixint[i]);

      lambda_max=30.;
      lambda_min=0.;
      while(true){
        lambda=(lambda_max+lambda_min)/2.;
        C_growth=1.-lambda*sqrt(pi)*exp(pow(lambda,2.))*erfc(lambda);
        func_lambda=2.*pow(lambda,2.)*C_growth;

        error_beta_fl=func_lambda-omega;

        if( fabs(error_beta_fl)<=ST_error*100.){
          break;
        }

        if(error_beta_fl>0.){
          lambda_max=lambda;
        }
        else{
          lambda_min=lambda;
        }
      }

//printf("Lumpda: %E %E %E %E %E \n",C_growth,omega,chem_composition[i],chem_comp_gamma[i],chem_comp_gp[i]);

    // **** Morral et al., Scripta, 30, 1994, 905 ****
//      eq_chemsq_dM+=pow((chem_comp_gp[i]-chem_comp_gamma[i]),2)/mobility[i];

      j=i;
//      for(j=1;j<no_elements;j++){
        if(model_velocity==0){
          eq_chemsq_dM+=(chem_comp_gp[i]-chem_comp_gamma[i])*(chem_comp_gp[j]-chem_comp_gamma[j])/mobility[i];
        }
        else if(model_velocity==1){
          eq_chemsq_dM+=C_growth*(chem_comp_gp[i]-chem_comp_gamma[i])*(chem_comp_gp[j]-chem_comp_gamma[j])/(mobility[i]*chem_comp_gamma[i] );
        }
        else if(model_velocity==2){
          eq_chemsq_dM+=C_growth*(chem_comp_gp[i]-chem_comp_gamma[i])*(chem_comp_gp[j]-chem_comp_gamma[j])/mobility[i];
        }
        else{
          diff[i]=0.;
          diff[i]=diffusivity[1][i];  //(diffusivity[0][i]+diffusivity[1][i])/2.;

          eq_chemsq_dM+=diff[i]*(chem_composition[j]-chem_comp_gamma[j])/(C_growth*(chem_comp_gp[i]-chem_comp_gamma[i]));
//          eq_chemsq_dM+=diff*(chem_composition[j]-chem_comp_matrixint[j])/(C_growth*(chem_comp_gp[i]-chem_comp_matrixint[i]));
        }
//      }

//      eq_chemsq_dM+=C_growth*(chem_comp_gp[i]-chem_comp_gamma[i])*(chem_comp_gp[i]-chem_comp_gamma[i])/mobility[i];

//      printf("eq_chemsq: %E %E %E %E \n",eq_chemsq_dM,C_growth,(chem_comp_gp[i]-chem_comp_gamma[i]),mobility[i]);

      // ---- from TC_presima ----
      beta_trans+=pow((chem_comp_gp[i]-chem_comp_gamma[i]),2)/(chem_comp_gamma[i]*diff[i]);
//      beta_trans+=pow((chem_comp_gp[i]-chem_composition[i]),2)/(chem_composition[i]*diffusivity[i]);
    }

    if(eq_chemsq_dM<0.){
      fprintf(stdout,"Cr escaping rate is higher than Al increasing rate: %f \n", eq_chemsq_dM);
      for(ierror=0;ierror<no_outputs;ierror++){
        results[count_results+ierror]=0.;
      }
      continue;
    }

    beta_trans=(4.*pi*pow(gp_radius0,2)/pow(gp_LP,4))/beta_trans;
//    printf("Beta: %E, %E, %E, %E, %E \n",pow(gp_radius0,2),pow(gp_LP,4),beta_trans,chem_comp_gp[1],chem_comp_gamma[1]);

    no_nu_site=no_ini_nu_site*Zeldfactor*beta_trans*exp(-dGM_critical/(KB*temperature));
    printf("# of initial-steady nucleation site: %E, %E, %E, %E \n",no_nu_site,no_ini_nu_site,Zeldfactor,beta_trans);
    if(no_nu_site<1E5){
      fprintf(stdout,"# of initial-steady nucleation site is too small: %E \n", no_nu_site);
      for(ierror=0;ierror<no_outputs;ierror++){
        results[count_results+ierror]=0.;
      }
      continue;
    }

    if(model_velocity==0 || model_velocity==2){
      K_LSW=c_LSW*mvolume_gp*E_intfac/eq_chemsq_dM;
    }
    else if(model_velocity==1){
      K_LSW=1./eq_chemsq_dM;
    }
    else{
      K_LSW=eq_chemsq_dM/(no_elements-1);
    }
//    printf("K_LSW: %E, %E, %E \n",K_LSW,mvolume_gp,eq_chemsq_dM);


//    time_incu=1./(2.*pow(Zeldfactor,2)*beta_trans);
//      time_incu=2./(pi*pow(Zeldfactor,2)*beta_trans);
    time_incu=1.5/(pow(Zeldfactor,2)*beta_trans);
    printf("Incubation time: %E, %E %E \n", time_incu,Zeldfactor,beta_trans);
    if(time_incu<0.){
      fprintf(stdout,"Incubation time is negative: %E \n", time_incu);
      for(ierror=0;ierror<no_outputs;ierror++){
        results[count_results+ierror]=0.;
      }
      continue;
    }


//    Velocity_int=( pow((pow(gp_radius0,3)+K_LSW*dtime),(1./3.)) - gp_radius0)/dtime;
//    printf("Interface Velocity: %E \n",Velocity_int);

    weak_stress=0.;
    strong_stress=0.;
    strong_stress_p=0;
    Orwan_stress_p=0.;
    YS=0.;
    YS_p=0.;
    gp_mean_radius=gp_radius0;
    gp_mean_radius_p=gp_radius0;
    tot_time=0.;
    tot_time_p=0.;
    no_prec_steps=max_time_step;
    onnucleation=0;
    oncoarsening=0;
    tot_precipitate[0]=0.;

// linear function to calculate the Gibbs free energy of the system
//printf("%E, %E, %E, %E \n",GM_eq,GM_L12,GM_A1,GM_A1-(GM_A1-GM_L12)*(chem_composition[1]-chem_comp_gamma[1])/(chem_comp_gp[1]-chem_comp_gamma[1]));

    for(i=0;i<max_time_step+1;i++){
      tot_precipitate[i]=0.;
      precipitates[1][i]=0.;
      precipitates[2][i]=0.;
      precipitates[3][i]=0.;
    }


    if(max_stress_model==0){
      for(time_step=1;time_step<=max_time_step;time_step++){
//printf("time step: %d \n",time_step);

        if(onnucleation==0 && oncoarsening==0){   // Nucleation
          dtime=time_incu*5E-4;
          tot_time+=dtime;
        }
        else if(onnucleation==1 && oncoarsening==0){  //Growth
          dtime=time_incu*1E-3;  //310*exp(-10/(time_step-onnustep));
          tot_time+=dtime;
        }
        else{
          if((time_step-onnustep)<50000){
            dtime=time_incu*0.1;
          }
          else{
            dtime=time_incu;
          }

          tot_time+=dtime;  // Coarsening
        }

//if(time_step%10==1 ){printf("(1) time: %d  %d: %d \n",generation,sam,time_step); }
//printf("(1) time: %d, %E, %E \n",time_step,dtime,tot_time);
    // ----------------------------------------- Model 1: Avrami type of simulation; the number of particles are not counted ---
        if(kinetics_model==0){
          gp_mean_radius=pow((pow(gp_radius0,3)+K_LSW*tot_time),(1./3.));
          Vf_gp=Vf_gp_eq*(1.-exp(-pi*no_nu_site*exp(-time_incu/tot_time)*pow(Velocity_int,3)*pow(tot_time,4)/3.));
          tot_Vf=Vf_gp;

	  weak_stress=3.06*EAPB*( ( sqrt(C1*EAPB*Vf_gp*gp_mean_radius/LT) - Vf_gp )/(2.*burgers))/1E+6;
	  strong_stress=3.06*C2*sqrt(Vf_gp*cw*cw*(pi*gp_mean_radius*EAPB/(cw*LT)-1.))*(shear_modulus*burgers/gp_mean_radius)/1E+6;

	}

    // -------------------- Model 2: volume fraction of gamma prime = (#nuclei site * one particle volume / total volume) -----
	else{

		tot_precipitate[time_step]=tot_precipitate[time_step-1];
		if(onnucleation==0){
		  tot_precipitate[time_step]=no_nu_site*exp(-time_incu/tot_time)*(tot_time-tot_time_p);
		  precipitates[1][time_step]=tot_precipitate[time_step]-tot_precipitate[time_step-1];

//		precipitates[1][time_step]=no_nu_site*exp(-time_incu/tot_time)*dtime;
//              if(precipitates[1][time_step]<1.){precipitates[1][time_step]=0.; }
//		tot_precipitate[time_step]=tot_precipitate[time_step-1]+precipitates[1][time_step];

		  precipitates[2][time_step]=gp_radius0;


		  if((precipitates[1][time_step]/tot_precipitate[time_step])<1E-4 || oncoarsening==1){
//		    printf("**********@ %d Nucleation -> Growth ********** \n",time_step);
		    no_prec_steps=time_step-1;
		    onnucleation=1;
		    onnustep=time_step;
		  }

		}

		weak_stress=0.0;
		strong_stress=0.0;
		tot_Vf=0.0;
		gp_mean_radius=0.;
		tot_Part=0.;


		no_calc_prec_steps=time_step;
		if(no_prec_steps<time_step){
			no_calc_prec_steps=no_prec_steps;
		}

		//printf("------- %d, %E, %E, %d, %d \n",time_step,precipitates[2][1],gp_mean_radius_p,oncoarsening,onnucleation);
//printf("(2) *************  \n");

		for(i=1;i<=no_calc_prec_steps;i++){

			dG_int=2.*mvolume_gp*E_intfac/precipitates[2][i];

                        if(model_velocity<3){
			  precipitates[2][i]=pow((pow(precipitates[2][i],3.)+K_LSW*dtime),(1./3.));   // radius of precipitates
                        }
                        else{
                          if(oncoarsening==0){
                            Velocity_int=K_LSW/precipitates[2][i];
                            precipitates[2][i]+=Velocity_int*dtime;
                          }
                          else{
                            Velocity_int=K_LSW/(gp_mean_radius_p*gp_mean_radius_p);
                            precipitates[2][i]+=Velocity_int*dtime;
                          }
                        }

			if(precipitates[2][i]<=0.){
				precipitates[2][i]=0.;
				precipitates[1][i]=0.;
			}

			precipitates[3][i]=4.*pi*pow(precipitates[2][i],3)/3.;  // volume of the precipitates in m^3
//			printf("nuclei & size & fraction %d: %E, %E, %E \n",i,precipitates[1][i],precipitates[2][i],precipitates[3][i]);

			Vf_gp=precipitates[1][i]*precipitates[3][i];
//			printf("Volume Fraction: %d, %E, %E, %E, %E \n",i,Vf_gp,precipitates[1][i],precipitates[2][i],precipitates[3][i]);

			tot_Vf+=Vf_gp;
			tot_Part+=precipitates[1][i];
			// =========================================================== For Mean Radius of GP ==== 
			//            gp_mean_radius+=precipitates[2][i]*Vf_gp;
			gp_mean_radius+=precipitates[2][i]*precipitates[1][i];

			if(i>=2){
				if(precipitates[2][i]>precipitates[2][i-1]){
					swap=0.;
					swap=precipitates[1][i];
					precipitates[1][i]=precipitates[1][i-1];
					precipitates[1][i-1]=swap;

					swap=0.;
					swap=precipitates[2][i];
					precipitates[2][i]=precipitates[2][i-1];
					precipitates[2][i-1]=swap;
				}
			}

			if(tot_Vf>Vf_gp_eq){

                          precipitates[1][i]-=(tot_Vf-Vf_gp_eq)/precipitates[3][i];
                          if(precipitates[1][i]<=0.){
                            precipitates[1][i]=0.;
                            no_prec_steps=i-1;
                          }
                          else{
                            no_prec_steps=i;
                          }
                          tot_Vf=Vf_gp_eq;
                        }
		}

		gp_mean_radius=gp_mean_radius/tot_Part;

		weak_stress=3.06*EAPB*( ( sqrt(C1*EAPB*tot_Vf*gp_mean_radius/LT) - tot_Vf )/(2.*burgers))/1E+6;
		strong_stress=3.06*C2*cw*sqrt(tot_Vf*(pi*gp_mean_radius*EAPB/(cw*LT)-1.))*(shear_modulus*burgers/gp_mean_radius)/1E+6;

                //Kozar et al, Strengthing Mechanisms..., MMT A, 2009
                  // 3D sphere
//                Orwan_stress=3.06*shear_modulus*burgers/(2.*(1./pow(tot_Vf,1./3.)-1.)*gp_mean_radius*1E+6);
                  // 2D plane: the radius of particle is ignored because of the gemortric error- d<0
                  Orwan_stress=3.06*shear_modulus*burgers/(sqrt(pi/tot_Vf)*gp_mean_radius*1E+6);

//               Orwan_stress=3.06*shear_modulus*burgers

//printf("Orowan: %E %E %E %E %E \n",Orwan_stress,shear_modulus,burgers,pow((4.*pi*pow(gp_mean_radius,3.)/(3.*tot_Vf)),1./3.),2.*gp_mean_radius);

                if ( gp_mean_radius != gp_mean_radius ){gp_mean_radius=precipitates[2][1];  }
		if ( (strong_stress != strong_stress) || (strong_stress<0.) ){strong_stress=0.; }
		if ( (weak_stress != weak_stress) || (weak_stress<0.) ){weak_stress=0.; }
                if ( (Orwan_stress != Orwan_stress) || (Orwan_stress<0.) ){Orwan_stress=0.; }

//printf("YS: %E, %E, %E, %E \n",weak_stress,strong_stress,Orwan_stress,Prec_stress);


                if( gp_mean_radius>=gp_mean_radius_p && (weak_stress*strong_stress*Orwan_stress!=0) ){
// *** Weak-Strong ***
                  if( (strong_stress<weak_stress) && ((strong_stress-strong_stress_p)<10) && (strong_stress>200.) ){
                    printf("Stress Maximized EXIT (1) \n");
                    printf("*-*-*-*-* Max Yield Stress: %E, %E, %E *-*-*-*-*\n",YS,YS_0+SS_stress+HP_stress+Back_stress,Prec_stress);
                    printf("Volume fraction, Mean radius & Processing time: %E, %E, %E \n",tot_Vf,gp_mean_radius,(tot_time+time_incu)/60.);
                    printf("P_stress: %E, %E, %E \n",weak_stress,strong_stress,Orwan_stress);
                    printf("=============================================================\n");
                    break;
                  }

// *** Orowan ***
                  if( (Orwan_stress<strong_stress) && (Orwan_stress<weak_stress) ){
                    printf("Stress Maximized EXIT (2) \n");
                    printf("*-*-*-*-* Max Yield Stress: %E, %E, %E *-*-*-*-*\n",YS,YS_0+SS_stress+HP_stress+Back_stress,Prec_stress);
                    printf("Volume fraction & Mean radius: %E, %E \n",tot_Vf,gp_mean_radius);
                    printf("P_stress: %E, %E, %E \n",weak_stress,strong_stress,Orwan_stress);
                    break;
                  }
                }

                gp_mean_radius_p=gp_mean_radius;
                strong_stress_p=strong_stress;
                Prec_stress=weak_stress;
                if( (strong_stress<weak_stress) && (strong_stress<Orwan_stress) ){
                  if(strong_stress>0.){ Prec_stress=strong_stress;}
                }
                if( (Orwan_stress<weak_stress) && (Orwan_stress<strong_stress) ){
                  if(Orwan_stress>0.){ Prec_stress=Orwan_stress;}
                }

                Orwan_stress_p=Orwan_stress;

//                HP_stress=C_HPeffect/sqrt(grain_size);
                HP_stress=0.;
                Back_stress=(shear_modulus*burgers*4./grain_size)/1E+6;

                SS_stress=0.;
                for(i=1;i<no_elements;i++){
                  SS_stress+=pow(YS_SS[i],2.)*chem_comp_trans[i];
                }

                dislocation_stress=0.25*Taylor_factor*shear_modulus*burgers*sqrt(1E+13)/1E+6;   // MPa

                SS_stress=sqrt(SS_stress);
                YS=YS_0+SS_stress+HP_stress+Back_stress+sqrt(pow(dislocation_stress,2)+pow(Prec_stress,2));
//                YS=(1.-tot_Vf)*(YS_0+SS_stress+HP_stress+Back_stress+sqrt(pow(dislocation_stress,2)+pow(Prec_stress,2)));
//                YS=YS_0+SS_stress+HP_stress+Back_stress+(1.-tot_Vf)*sqrt(pow(dislocation_stress,2)+pow(Prec_stress,2))+tot_Vf*dislocation_stress;
//                YS_p=YS;

		// ************************************************************************ Update Model Parameters ****
		// ------------------------------- After phase transformation, the chemical composition of gamma is updated

                if ( Vf_gp_eq-tot_Vf > 0.001 ) {
//printf("Vf: %E, %E, %E, %E \n",Vf_gp_eq,tot_Vf,Vf_gp_eq-tot_Vf,K_LSW );

		  for(i=no_elem_input;i<no_elements;i++)
		  {
		    chem_comp_trans[i]=chem_composition[i]-(mvolume_gamma*(chem_comp_gp[i]-chem_composition[i])*tot_Vf)/(mvolume_gp*(1.-tot_Vf));

//printf("composition %d: %E %E %E %E \n",i,chem_comp_trans[i],chem_comp_gamma[i],chem_composition[i],chem_comp_gp[i]);

           // If the molar volume is not obtained from database, chem_comp_trans[i] may be negative!!
                    if(chem_comp_trans[i]<chem_comp_gamma[i]){chem_comp_trans[i]=chem_comp_gamma[i];  }
		  }

                  beta_trans=0;
                  for(i=1;i<no_elements;i++){
                    beta_trans+=pow((chem_comp_gp[i]-chem_comp_gamma[i]),2)/(chem_comp_gamma[i]*diff[i]);
                  }

                  beta_trans=(4.*pi*pow(gp_radius0,2)/pow(gp_LP,4))/beta_trans;

                }
               else{
                 if(oncoarsening==0){
                   oncoarsening=1;
//                   printf("**************** COARSENING: %d ****************** \n",time_step);

                   eq_chemsq_dM=0.;
   K_LSW=4.*chem_comp_gamma[1]*diff[1]*(2*E_intfac*(mvolume_gp/Avogadro))/(27.*(chem_comp_gp[1]-chem_comp_gamma[1])*KB*temperature);
//                   printf("K_LSW: %E \n",K_LSW);
                 }
               }



        if(oncoarsening==0){

          eq_chemsq_dM=0.;
//          i=1;
          for(i=1;i<no_elements;i++){

          // **** Qing ****
            error_beta_fl=10000.;

            chem_comp_matrixint[i]=chem_comp_gamma[i]*exp(k_chem/gp_mean_radius);
//            printf("COMPOSITION %d: %E, %E %E \n",i,chem_comp_matrixint[i],chem_comp_gamma[i],chem_composition[i]);

            omega=(chem_comp_trans[i]-chem_comp_gamma[i])/(chem_comp_gp[i]-chem_comp_gamma[i]);
//            omega=(chem_comp_trans[i]-chem_comp_matrixint[i])/(chem_comp_gp[i]-chem_comp_matrixint[i]);

            if(oncoarsening==0){
              lambda_max=30.;
              lambda_min=0.;
              while(true){
                lambda=(lambda_max+lambda_min)/2.;
                C_growth=1.-lambda*sqrt(pi)*exp(pow(lambda,2.))*erfc(lambda);
                func_lambda=2.*pow(lambda,2.)*C_growth;

                error_beta_fl=func_lambda-omega;

                if( fabs(error_beta_fl)<=ST_error*100.){
                  break;
                }

                if(error_beta_fl>0.){
                  lambda_max=lambda;
                }
                else{
                  lambda_min=lambda;
                }
                if(lambda_max==lambda_min){
                  C_growth==1.;
                  break;
                }
              }
            }
            else{
              C_growth==1.;
            }

            j=i;
//            for(j=1;j<no_elements;j++){
              if(model_velocity==0){
                eq_chemsq_dM+=(chem_comp_gp[i]-chem_comp_gamma[i])*(chem_comp_gp[j]-chem_comp_gamma[j])/mobility[i];
              }
              else if(model_velocity==1){
                eq_chemsq_dM+=C_growth*(chem_comp_gp[i]-chem_comp_gamma[i])*(chem_comp_gp[j]-chem_comp_gamma[j])/(mobility[i]*chem_comp_gamma[i] );
              }
              else if(model_velocity==2){
                eq_chemsq_dM+=C_growth*(chem_comp_gp[i]-chem_comp_gamma[i])*(chem_comp_gp[j]-chem_comp_gamma[j])/mobility[i];
              }
              else{
                eq_chemsq_dM+=diff[i]*(chem_comp_trans[j]-chem_comp_gamma[j])/(C_growth*(chem_comp_gp[i]-chem_comp_gamma[i]));
//                eq_chemsq_dM+=diff*(chem_comp_trans[j]-chem_comp_matrixint[j])/(C_growth*(chem_comp_gp[i]-chem_comp_matrixint[i]));
//                printf("eq_chemsq_dM %d: %E %E %E \n",i,eq_chemsq_dM,C_growth,diff);
//                if(eq_chemsq_dM<0){
//                  printf("negative growth rate \n");
//                 continue;
//                }
              }
//            }
          }

          if(model_velocity==0 || model_velocity==2){
            K_LSW=c_LSW*mvolume_gp*E_intfac/eq_chemsq_dM;
          }
          else if(model_velocity==1){
            K_LSW=1./eq_chemsq_dM;
          }
          else{
            K_LSW=eq_chemsq_dM/(no_elements-1);
          }
        }

          tot_time_p=tot_time;
        }
      }
    }
    else{

    // -+-+-+ The strong precipitate stress is the maximum value of the strong coupling model -+-+-+
    // Crudden et al., Acta, 2014
    // under the conditions (RC Reed, P81), model (0)strong_stress=381; (1)strong_stress=336

      strong_stress=3.06*0.5*EAPB*sqrt(Vf_gp)/burgers/1E+6;
      Prec_stress=strong_stress;
    }


    if ( strong_stress != strong_stress ) {
      printf("The Yield Stress calculation is not accurate!! \n");
      for(ierror=0;ierror<no_outputs;ierror++){
        results[count_results+ierror]=0.;
      }
      continue;
    }

// %%%% Output of Results %%%%
    if ( ((tot_time+time_incu)/60.) < 1.  ) {
      printf("Processing time is too short!!  \n");
      for(ierror=0;ierror<no_outputs;ierror++){
        results[count_results+ierror]=0.;
      }
    }
    else{
      results[count_results]=VF[ind_L12];
      results[count_results+1]=gp_mean_radius;
      results[count_results+2]=density_all;
      results[count_results+3]=YS_0+SS_stress+HP_stress+Back_stress;
      results[count_results+4]=Prec_stress;

      results[count_results+5]=YS;
      results[count_results+6]=Vf_gp_eq;
      results[count_results+7]=E_intfac;
      results[count_results+8]=EAPB;
      results[count_results+9]=(tot_time+time_incu)/60.;
    }

  }

  return results;

}



