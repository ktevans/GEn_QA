#include <iostream>

#include "../include/Analysis.h"
#include "../include/FF_param.h"


// This script calculates the variables for GE/GM extraction.
void UpdateAverageKinematics(analyzed_tree *T){
 
  SetHe3PolAngle();  // Sets the He3 polarization value/direction

  TLorentzVector Pe(0,0,T->ebeam,T->ebeam);   // incoming e-
  TLorentzVector Peprime(T->trPx,   // scattered e-
			 T->trPy,
			 T->trPz,
			 T->trP);
  
  // Calculate the scattering plane, which we will use to get angles
  TLorentzVector q = Pe - Peprime;
  TVector3 Pe_vect = Pe.Vect();
  TVector3 Peprime_vect = Peprime.Vect();
  TVector3 q_vect = q.Vect();
  TVector3 normal = Pe_vect.Cross(Peprime_vect);
  normal = normal.Unit();
  q_vect = q_vect.Unit();
  
  // Get our kinematics for the event
  double m = constant::Mp;
  double tau = T->Q2 / (4 * m * m);   // Defenition of tau
  double tan_theta = tan(T->etheta/2);
  double epsilon = 1.0 / (1 + 2*(1 + tau)*pow(tan_theta,2)); // Defenition

  // These two give the polarization components 
  double Px = normal.Dot(q_vect.Cross(TargetPolDirection)); 
  double Pz = q_vect.Dot(TargetPolDirection);

  // These variables are for fifth order Taylor expansion
  double B = -2 * sqrt( tau * (1 + tau) ) * tan_theta * Px;
  double C = -2 * tau * sqrt(1 + tau + pow((1 + tau) * tan_theta,2) ) * tan_theta * Pz; 
  double D = tau / epsilon;

  double T_0 = C / D;
  double T_1 = B / D;
  double T_2 = -1*C / (D*D);
  double T_3 = -1*B / (D*D);
  double T_4 = C / (D*D*D);
  double T_5 = B / (D*D*D);

  // Now we update the average for each value
  count_avg_exp++;
  T_avg[0] += (T_0 - T_avg[0]) / count_avg_exp;
  T_avg[1] += (T_1 - T_avg[1]) / count_avg_exp;
  T_avg[2] += (T_2 - T_avg[2]) / count_avg_exp;
  T_avg[3] += (T_3 - T_avg[3]) / count_avg_exp;
  T_avg[4] += (T_4 - T_avg[4]) / count_avg_exp;
  T_avg[5] += (T_5 - T_avg[5]) / count_avg_exp;
  T_1_Q2_avg += (T_1*T->Q2 - T_1_Q2_avg) / count_avg_exp;
  //Q2_avg = T_1_Q2_avg / T_avg[1];

  // We also get the average kinematic variables for later
  Q2_avg += (T->Q2 - Q2_avg) / count_avg_exp;
  tau_avg += (tau - tau_avg) / count_avg_exp;
  epsilon_avg += (epsilon - epsilon_avg) / count_avg_exp;
  Px_avg += (Px - Px_avg) / count_avg_exp;
  Pz_avg += (Pz - Pz_avg) / count_avg_exp;

}


// This gets the GE/GM value from parameterizations from a Q2 value
// This can be found in literature or Seamus/Freddy's theses 
void GetGEGMFromTheory(bool is_neutron, double Q2, double &R, double &R_err){

  double m = constant::Mp;
  double tau = Q2 / (4 * m * m);
  double GD = pow(1.0 + Q2/(0.71), -2.0);  // Dipole approximation
  double GE,GM,GE_err,GM_err;

  if(is_neutron){ // Neutron parameterizations
    // Seamus Fit
    //GE = (1.520*tau + 2.629*tau*tau + 3.055*tau*tau*tau)*GD/(1.0+5.222*tau+0.040*tau*tau+11.438*tau*tau*tau);
    // Kelly
    //GM = -1.913*(1.0+2.33*tau)/(1.0 + 14.72*tau + 24.20*tau*tau + 84.1*tau*tau*tau );
    // Instead of Seamus or Kelly we use Ye global fit
    GetGENFit(Q2,GE,GE_err);
    GetGMNFit(Q2,GM,GM_err);
  }
  else{ // Proton parameterizations
    // Kelly
    //GE = (1.0-0.24*tau)/(1.0 + 10.98*tau + 12.82*tau*tau + 21.97*tau*tau*tau );
    // Kelly
    //GM = 2.79*(1.0+0.12*tau)/(1.0 + 10.97*tau + 18.86*tau*tau + 6.55*tau*tau*tau );
    // Instead of Seamus or Kelly we use Ye global fit
    GetGEPFit(Q2,GE,GE_err);
    GetGMPFit(Q2,GM,GM_err);
  }

  R = GE / GM;
  R_err = R*sqrt( pow(GE_err/GE,2) + pow(GM_err/GM,2) ); 
}

// Given the average kinematic values get the asymmetry from this Q2 value
// by following the GE/GM parameterizations
void GetAFromQ2(bool is_neutron, double Q2, double &A, double &A_err){

  if(epsilon_avg == 0 || tau_avg == 0){
    cout<<"Error: [Analysis::GetAFromQ2] epsilon and tau values have not been set!!!"<<endl;
    exit(0);
  }

  double R, R_err;
  GetGEGMFromTheory(is_neutron, Q2_avg, R, R_err); // Get GE/GM from parameterization
  double A_par = epsilon_avg / tau_avg;
  double B_par = sqrt(2*epsilon_avg*(1-epsilon_avg)/tau_avg)*Px_avg;
  double C_par = sqrt(1 - epsilon_avg*epsilon_avg)*Pz_avg;

  if(!is_neutron){
    GEGMp = R;
    GEGMp_err = R_err;
  }
  
  A = -1.0 / ( 1 + A_par * R* R)* ( B_par * R + C_par );
  A_err = (2*A_par*R*A + B_par)/(1 + A_par*R*R) * R_err;
}

// From the average kinematic values and the asymmetry measurement calculate
// GE/GM
double GetGEGMFromA(double A_phys, double A_stat_err, double A_sys_err){

  if(epsilon_avg == 0 || tau_avg == 0){
    cout<<"Error: [Analysis::GetGEGMFromA] epsilon and tau values have not been set!!!"<<endl;
    exit(0);
  }
  
  double e = epsilon_avg;
  double t = tau_avg;

  double A = e / t * A_phys;
  double B = Px_avg*sqrt(2*e*(1 - e) / t);
  double C = A_phys + Pz_avg*sqrt(1 - e*e);

  A_quad = A;
  B_quad = B;
  C_quad = C;
  
  // We are solving the quadratic equation to get GE/GM
  double R1 = (-B + sqrt(B*B - 4*A*C)) / (2*A);
  double R2 = (-B - sqrt(B*B - 4*A*C)) / (2*A);
  GEGMn = R1;  // I think this root is correct
  
  // Now we calculate the error on the quadratic formula above
  // See thesis for derivation of this calculation
  GEGMn_stat_err = sqrt( pow(-1*C / (A*sqrt(B*B - 4*A*C)) + (B - sqrt(B*B - 4*A*C)) / (2*A*A),2) * e*e/(t*t) + 1.0 / (B*B - 4*A*C) ) * A_stat_err;
  GEGMn_sys_err = sqrt( pow(-1*C / (A*sqrt(B*B - 4*A*C)) + (B - sqrt(B*B - 4*A*C)) / (2*A*A),2) * e*e/(t*t) + 1.0 / (B*B - 4*A*C) ) * A_sys_err;

  return GEGMn;
}


// Here we use Newtons method to find the 0 point of the asymmetry expansion
void GetGEGMFromA_old(double A){

  double gamma = (A - T_avg[0]) / T_avg[1];  //  Start with first order apporximation
  int niter = 1000;
  bool converge = false;

  for(int inum = 0; inum < niter; inum++){
    double fsum = 0;
    double fderivsum = 0;
    
    for(int i=0; i < nexp; i++){
      fsum += T_avg[i]*pow(gamma,i);
      if(i > 0) fderivsum += i*T_avg[i]*pow(gamma,i-1);
    }
    
    double f = A - fsum;
    double fderiv = -1*fderivsum;
    double gamma_old = gamma;
    gamma = gamma_old - f / fderiv;
    
    converge = abs((gamma - gamma_old) / gamma) < 10e-6;
    if(converge) break;
  }

  if(!converge){
    cout<<"Error: [Analysis::GetGEGMFromA] Newtons method did not converge!!!"<<endl;
    exit(0);
  }

  GEGMn = gamma;   //Save the form factor result
}


// This function is used as a test to see what GE/GM we expect from A just from
// the theory
double GetAFromGEGM_theory(double GEGM,double etheta, double Q2){

  double m = constant::Mp;
  double tau = Q2 / (4 * m * m);

  double num = -2 * GEGM * sqrt(tau*(tau + 1)) * tan(etheta/2);
  double denom = GEGM*GEGM + tau + 2*tau*(tau + 1)*pow(tan(etheta/2),2);
  double A = num / denom;

  return A;
}


double distribution_fits::fitbg_pol4( double *x, double *par ){
  double dx = x[0];

  // Let's use 4th order polynomial for the background:
  double bg = 0.0;
  for( int i = 0; i<5; i++ ){
    bg += par[i]*pow(dx,i);
  }

  return bg; 
}

double distribution_fits::fitbg_pol3( double *x, double *par ){
  double dx = x[0];

  // Let's use 3rd order polynomial for the background:
  double bg = 0.0;
  for( int i = 0; i<4; i++ ){
    bg += par[i]*pow(dx,i);
  }

  return bg; 
}

double distribution_fits::fitbg_pol2( double *x, double *par ){
  double dx = x[0];

  // Let's use 3rd order polynomial for the background:
  double bg = 0.0;
  for( int i = 0; i<3; i++ ){
    bg += par[i]*pow(dx,i);
  }

  return bg; 
}

double distribution_fits::fitbg_gaus( double *x, double *par ){
  double dx = x[0];

  // Let's use gaussian for the background:
  double bg = par[0]*exp(-0.5*pow((dx-par[1])/par[2],2));

  return bg; 
}

// Custom fit for the He3 distributions. This assumes we have a proton,
// neutron, and background distribution. These data is then fitted to these
// three distributions
double distribution_fits::fitsim( double *x, double *par){
  double dx = x[0];

  double Norm_overall = par[0]; // Normalization for the spectrum
  double R_pn = par[1];         // Ratio between proton/neutron peak
  double Bg_norm = par[2];      // Normalization of the background

  double bg = 0.0;
  
  // Get the background function:
  if(bg_shape_option == "pol4") bg = fitbg_pol4(x,&par[3]);
  else if(bg_shape_option == "pol3") bg = fitbg_pol3(x,&par[3]);
  else if(bg_shape_option == "pol2") bg = fitbg_pol2(x,&par[3]);
  else if(bg_shape_option == "gaus") bg = fitbg_gaus(x,&par[3]);
  else if(bg_shape_option == "from data") bg = hdx_bg_data->Interpolate(dx);
  
  // The fit is now the sum of these three distributions
  // fit = N * (p_dist + R * n_dist + N_bg * bg_dist)
  double simu = Norm_overall * (hdx_sim_p->Interpolate(dx) + R_pn * hdx_sim_n->Interpolate(dx) + Bg_norm * bg);
  return simu;   
}


// This function is used to fit the He3 data
void distribution_fits::He3_fit_dists(){

  ////////////First do some checks to make sure that things are set///////
  ////////////////////////////////////////////////////////////////////////
  // Check bg shape is set
  if(!bg_shape_set){
    cout<<"Error: [distribution_fits::He3_sim_fit] bg shape has not been set!!!"<<endl;
    exit(0);
  }
  else {
    if(bg_shape_option != "pol2" && bg_shape_option != "pol3" && bg_shape_option != "pol4" && bg_shape_option != "gaus" && bg_shape_option != "from data"){
      cout<<"Error: [distribution_fits::He3_sim_fit] bg shape option " + bg_shape_option + " is not supported!!!"<<endl;
      exit(0);
    }
  }

  //Check if all the histograms are set
  if(hdx_data == NULL){
    cout<<"Error: [distribution_fits::He3_sim_fit] hdx_data has not been set!!!"<<endl;
    exit(0);
  }
  if(hdx_sim_p == NULL){
    cout<<"Error: [distribution_fits::He3_sim_fit] hdx_sim_p has not been set!!!"<<endl;
    exit(0);
  }
  if(hdx_sim_n == NULL){
    cout<<"Error: [distribution_fits::He3_sim_fit] hdx_sim_n has not been set!!!"<<endl;
    exit(0);
  }
  if(hdx_bg_data == NULL){
    cout<<"Error: [distribution_fits::He3_sim_fit] hdx_bg_data has not been set!!!"<<endl;
    exit(0);
  }

  ////////////////////////////////////////////////////////////////////////
  
  //Normalize the histograms so they are all of a similar scale
  double scale = hdx_data->Integral();
  hdx_data->Scale(1.0/hdx_data->Integral());
  hdx_sim_p->Scale(1.0/hdx_sim_p->Integral());
  hdx_sim_n->Scale(1.0/hdx_sim_n->Integral());
  hdx_bg_data->Scale(1.0/hdx_bg_data->Integral());

  //Get histogram data so all future histograms match
  int nbins = hdx_data->GetNbinsX();
  double xmin = hdx_data->GetXaxis()->GetBinLowEdge(1);
  double xmax = hdx_data->GetXaxis()->GetBinUpEdge(nbins);
  
  //Set background type and check if the option exists
  int npar = -1;
  if(bg_shape_option == "pol4") npar = 3 + 5;
  else if(bg_shape_option == "pol3") npar = 3 + 4;
  else if(bg_shape_option == "pol2") npar = 3 + 3;
  else if(bg_shape_option == "gaus") npar = 3 + 3;
  else if(bg_shape_option == "from data") npar = 3;  //Our only bg par is the scale of the bg

  //Set fit to function fitsim
  TF1 *FitFunc = new TF1( "FitFunc", this, &distribution_fits::fitsim,xmin,xmax,npar);
      
  //Set some arbitrary starting values
  FitFunc->SetNpx(1000);
  double startpar[] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
  FitFunc->SetParameters(startpar);
  FitFunc->SetParLimits(0,0.1,100);
  FitFunc->SetParLimits(1,0.1,100);
  FitFunc->SetParLimits(2,0,100);
    
  hdx_data->Fit(FitFunc,"q0","",xmin,xmax);
    
  //Scale the p/n functions
  //This assumes our fit was like N*(p_dist + R*n_dist + N_bg * bg_dist)
  hdx_sim_p->Scale( FitFunc->GetParameter(0) );
  hdx_sim_n->Scale( FitFunc->GetParameter(0)*FitFunc->GetParameter(1) );

  hdx_bg_fit = new TH1F("hdx_bg","",nbins,xmin,xmax);
  hdx_total_fit = new TH1F("hdx_total_fit","",nbins,xmin,xmax);

  //Set the backgrond function using the fit we already found
  TF1 *fit_bg;

  if(bg_shape_option == "pol4") 
    fit_bg = new TF1("fit_bg",this, &distribution_fits::fitbg_pol4,xmin,xmax,npar - 3);
  else if(bg_shape_option == "gaus") 
    fit_bg = new TF1("fit_bg",this, &distribution_fits::fitbg_gaus,xmin,xmax,npar - 3);
  else if(bg_shape_option == "pol3")
    fit_bg = new TF1("fit_bg",this, &distribution_fits::fitbg_pol3,xmin,xmax,npar - 3);
  else 
    fit_bg = new TF1("fit_bg",this, &distribution_fits::fitbg_pol2,xmin,xmax,npar - 3);

  for(int i=0; i<npar - 3; i++)
    fit_bg->SetParameter(i,FitFunc->GetParameter(i+3));
    
  //Fill the bg hist and total hist
  for(int ibin = 0; ibin < nbins;ibin++){
    //Again normalize bg assuming our fit is N*(p_dist + R*n_dist + N_bg * bg_dist)

    if(bg_shape_option == "from data"){
      hdx_bg_fit->SetBinContent(ibin,FitFunc->GetParameter(0) * FitFunc->GetParameter(2) * hdx_bg_data->GetBinContent(ibin));
    }
    else {
      hdx_bg_fit->SetBinContent(ibin,FitFunc->GetParameter(0) * FitFunc->GetParameter(2) * fit_bg->Eval(hdx_total_fit->GetBinCenter(ibin)));
    }
      
	
    hdx_total_fit->SetBinContent(ibin,hdx_sim_n->GetBinContent(ibin) + hdx_sim_p->GetBinContent(ibin) + hdx_bg_fit->GetBinContent(ibin));
  }


  // Last step is to scale everything back up to the original data scale
  hdx_data->Scale(scale);
  hdx_sim_p->Scale(scale);
  hdx_sim_n->Scale(scale);
  hdx_bg_fit->Scale(scale);
  hdx_total_fit->Scale(scale);

}


// This function sets the beam polarization for 
void analyzed_info::SetPolBeam(int key, TDatime *evtime, vector<vector<TDatime*>> PolTime, vector<vector<double>> PolVal) {
  
  if(!Asym_runs[key].run_init){
    
    time_t time = evtime->Convert();
    
    int pol_index = -1;
    for(int ipol=0; ipol < PolTime.size(); ipol++){
      time_t time1 = PolTime[ipol][0]->Convert();
      time_t time2 = PolTime[ipol][1]->Convert();
      time_t diff1 = time - time1;
      time_t diff2 = time - time2;
      
      if(diff1 > 0 && diff2 < 0) pol_index = ipol;
      
    }
    
    Asym_runs[key].P_beam[0] = PolVal[pol_index][0] / 100;  //conver to decimal
    Asym_runs[key].P_beam[1] = PolVal[pol_index][1] / 100;  //conver to decimal
    
    Asym_runs[key].run_init = true;
  } 
}


// This function sets the beam polarization for 
void analyzed_info::SetAvgPol(){
    double weightHe3vals = 0;
    double weightBeamvals = 0;
    double weightHe3errs = 0;
    double weightBeamerrs = 0;
    double weightsum = 0;
    
    for (const auto& pair : Asym_runs) {
      run_info Asym_runs = pair.second;
      double weight = Asym_runs.N_raw_p + Asym_runs.N_raw_m;
      
      weightsum += weight;

      weightHe3vals += weight*Asym_runs.P_He3;
      weightBeamvals += weight*Asym_runs.P_beam[0];

      weightHe3errs += weight*Asym_runs.P_He3_err;
      weightBeamerrs += weight*Asym_runs.P_beam[1];
    }
    
    P_He3_avg = weightHe3vals / weightsum;
    P_beam_avg = weightBeamvals / weightsum;
    P_He3_avg_err = weightHe3errs / weightsum;
    P_beam_avg_err = weightBeamerrs / weightsum;

  }



// This function calculates the asymmetry values. This formalism follows a more
// standard formalism of calculating each asymmetry contribution and scaling 
// it by the fraction inside the quasielastic cuts
void analyzed_info::CalcAsymVals(){

  ////////////////// Make sure all the variables have been loaded /////////////
  if(A_acc < -100 || A_acc_err < -100){
    cout<<"Warning: [analyzed_info::CalcAsymVals] A_acc or A_acc_err are not set!!!"<<endl;
  }
  if(f_acc < -100 || f_acc_err < -100){
    cout<<"Warning: [analyzed_info::CalcAsymVals] f_acc or f_acc_err are not set!!!"<<endl;
  }
  if(A_pion < -100 || A_pion_err < -100){
    cout<<"Warning: [analyzed_info::CalcAsymVals] A_pion or A_pion_err are not set!!!"<<endl;
  }
  if(f_pion < -100 || f_pion_err < -100){
    cout<<"Warning: [analyzed_info::CalcAsymVals] f_pion or f_pion_err are not set!!!"<<endl;
  }
  if(A_in < -100 || A_in_err < -100){
    cout<<"Warning: [analyzed_info::CalcAsymVals] A_in or A_in_err are not set!!!"<<endl;
  }
  if(f_in < -100 || f_in_err < -100){
    cout<<"Warning: [analyzed_info::CalcAsymVals] f_in or f_in_err are not set!!!"<<endl;
  }
  if(f_N2 < -100 || f_N2_err < -100){
    cout<<"Warning: [analyzed_info::CalcAsymVals] f_N2 or f_N2_err are not set!!!"<<endl;
  }
  if(N_proton == 0){
    cout<<"Warning: [analyzed_info::CalcAsymVals] N_proton is zero!!!"<<endl;
  }
  /////////////////////////////////////////////////////////////////////////


  double Sigma_raw_tot = 0;
  
  for (const auto& pair : Asym_runs) {
    run_info Asym_runs = pair.second;

    Sigma_raw_tot += Asym_runs.N_raw_p;
    Sigma_raw_tot += Asym_runs.N_raw_m;
    N_raw_p += Asym_runs.N_raw_p;
    N_raw_m += Asym_runs.N_raw_m;
  }
  
  // Get the values from the raw asymmetry
  double Sigma_raw = N_raw_p + N_raw_m;
  double Delta_raw = N_raw_p - N_raw_m;
  
  // Calculate the fractions
  f_p = 1.0*N_proton / Sigma_raw_tot;
  f_FSI = 0.0287; // This value comes from GEN-I
  f_FSI_err = 0.0026;   // This value comes from GEN-I
  // The rest are calculated in previous scripts

  f_n = 1 - f_p - f_acc - f_N2 - f_pion - f_in - f_FSI;
  
  // Calculate the asymmetry values
  A_raw = Delta_raw / Sigma_raw;
  A_raw_err = CalcAsymErr(N_raw_p, N_raw_m);
  A_FSI = 0.0003; // Number comes from GEN-I
  A_FSI_err = 0.0005; // Number comes from GEN-I
  
  SetAvgPol();  // Get avg polarization for proton asym

  // Calculate the proton asymmetry from the expansion
  //A_p = 0;
  //for(int i=0; i < nexp; i++)
  //A_p += T_avg[i] * pow(GetGEGMFromTheory(0,Q2_avg),i);

  // Calculate proton asym from averages
  GetAFromQ2(0, Q2_avg, A_p_phys, A_p_phys_err);

  A_p = P_He3_avg*P_beam_avg*P_p*A_p_phys;
  A_p_err = A_p * sqrt( pow(A_p_phys_err/A_p_phys,2) + pow(P_beam_avg_err/P_beam_avg,2) + pow(P_He3_avg_err/P_He3_avg,2) );
  
  f_p_err = CalcFractionErr(N_proton, Sigma_raw_tot);
  //!!!! The rest are calculated in previous scripts !!!!!!!!!///////
  
  // Run summation to get total asymmetry
  double A_phys_num = 0;
  double A_phys_stat_sum = 0;
  double A_phys_sys_sum = 0;

  // Loop over all runs
  for (const auto& pair : Asym_runs) {
    run_info Asym_runs = pair.second;
    
    // Skip if run has no events in it
    if(Asym_runs.N_raw_p == 0 || Asym_runs.N_raw_m == 0) continue;

    // Calculate the raw asymmetry and the error
    double A_raw_run = 1.0*(Asym_runs.N_raw_p - Asym_runs.N_raw_m) / (Asym_runs.N_raw_p + Asym_runs.N_raw_m);
    double A_raw_err_run = CalcAsymErr(Asym_runs.N_raw_p, Asym_runs.N_raw_m);

    // Calculate the physical asymmetry and statistical error
    double denom = Asym_runs.P_He3 * Asym_runs.P_beam[0] * P_n * f_n;
    double A_phys_run = ( A_raw_run - f_acc*A_acc - f_p*A_p - f_pion*A_pion - f_in*A_in - f_FSI*A_FSI ) / denom;
    double A_stat_err_run = A_raw_err_run / denom;

    // These sums are needed to get the total asym from all runs
    A_phys_num += A_phys_run / (A_stat_err_run*A_stat_err_run);
    A_phys_stat_sum += 1.0 / (A_stat_err_run*A_stat_err_run);
    
  }
  
  // From all runs combined we now calculate the total
  A_phys = A_phys_num / A_phys_stat_sum;

  // This is the statistical error
  A_phys_stat_err = 1 / sqrt(A_phys_stat_sum);

  // Here we separate parts of the systematic error because they are long //////

  // Asymmetry systematic error squared
  double A_sys_err = (f_acc*f_acc*A_acc_err*A_acc_err + f_pion*f_pion*A_pion_err*A_pion_err + f_in*f_in*A_in_err*A_in_err + f_p*f_p*A_p_err*A_p_err + f_FSI*f_FSI*A_FSI_err*A_FSI_err) / pow(P_He3_avg*P_beam_avg*P_n*f_n,2);
  
  // Fraction systematic error squared
  double f_sys_err = pow(A_phys/f_n*f_N2_err,2) + pow(CalcAsymFractionErr(A_acc, f_acc_err),2) + pow(CalcAsymFractionErr(A_pion, f_pion_err),2) + pow(CalcAsymFractionErr(A_in, f_in_err),2) + pow(CalcAsymFractionErr(A_p, f_p_err),2) + pow(CalcAsymFractionErr(A_FSI, f_FSI_err),2);

  // Polarization systematic error squared
  double P_sys_err = A_phys*A_phys * pow(P_He3_avg_err/P_He3_avg,2) + pow(P_beam_avg_err/P_beam_avg,2);

  // Total systematic error adding up the three parts above
  A_phys_sys_err = sqrt(A_sys_err + f_sys_err + P_sys_err);
  
  // The final error is the statistical + systematic
  A_phys_tot_err = sqrt(A_phys_stat_err*A_phys_stat_err + A_phys_sys_err*A_phys_sys_err);

  // Get final GE/GM and GE
  GetGEGMFromA(A_phys, A_phys_stat_err, A_phys_sys_err); // Get form factor result

  GetGMNFit(Q2_avg,GMn,GMn_err);  // Get GMn from parameterization

  GEn = GEGMn * GMn;
  GEn_stat_err = sqrt(pow(GEGMn*GMn_err,2) + pow(GMn*GEGMn_stat_err,2));
  GEn_sys_err = sqrt(pow(GEGMn*GMn_err,2) + pow(GMn*GEGMn_sys_err,2));
}


