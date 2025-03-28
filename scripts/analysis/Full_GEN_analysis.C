//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//   Created by Sean Jeffas, sj9ry@virginia.edu
//   Last Modified August 8, 2024
//
//
//   The purpose of this script is to complete the full GEN
//   analysis in a single script. It sets a few kinematic cuts
//   and then runs all the correction scripts in order and
//   writes the output into a text file.
//
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

#include "../../include/gen-ana.h"

// Input types: GEN2, GEN3, GEN4
void Full_GEN_analysis(TString cfg){

  TString DB_dir = DBparse::DB_corr_dir;

  double W2step = 0.1;
  double dystep = 0.1;

  double W2min = 1.2;
  double W2max = 2.0;

  double dymin = 0.1;
  double dymax = 0.5;

  double W2stable;
  double dystable;

  if(cfg == "GEN2"){
    W2min = 1.2;
    W2max = 2.0;
    W2stable = 1.6;

    dymin = 0.3;
    dymax = 1.0;
    dystable = 0.5;
  }
  else if(cfg == "GEN3"){
    W2min = 1.6;
    W2max = 2.2;
    W2stable = 2.0;

    dymin = 0.1;
    dymax = 0.5;
    dystable = 0.4;
  }
  else if(cfg == "GEN4"){
    W2min = 1.6;
    W2max = 2.2;
    W2stable = 2.0;

    dymin = 0.1;
    dymax = 0.5;
    dystable = 0.3;
  }
  else{
    cout<<"Error: cfg value is not accepted!!!"<<endl;
    exit(0);
  }

  int ndy = (dymax - dymin) / dystep + 1;
  int nW2 = (W2max - W2min) / W2step + 1;
  
  for(int iW2 = 0; iW2 < nW2; iW2++){
    double W2 = W2min + iW2*W2step;

    TString file_old = cfg + "_corr.csv";
    TString file = Form(cfg + "_corr_W2_%g_%g_dy_%g.csv",-2.0,W2,dystable);

    TString command = "cp " + DB_dir + file_old + " " + DB_dir + file;
    gSystem->Exec(command); // Create the DB file to start

    TString accidental = Form("root -b -q 'accidental_contamination.C(\"" + cfg + "\",%g,%g,%g)'",-2.0,W2,dystable);
    TString nitrogen= Form("root -b -q 'Nitrogen_contamination.C(\"" + cfg + "\",%g,%g,%g)'",-2.0,W2,dystable);
    TString pion = Form("root -b -q 'pion_contamination.C(\"" + cfg + "\",%g,%g,%g)'",-2.0,W2,dystable);
    TString inelastic = Form("root -b -q 'Inelastic_contamination.C(\"" + cfg + "\",%g,%g,%g)'",-2.0,W2,dystable);
    TString GEN_calc = Form("root -b -q 'GEN_Extraction.C(\"" + cfg + "\",%g,%g,%g)'",-2.0,W2,dystable);

    
    gSystem->Exec(accidental); // Run the accidental contamination
    gSystem->Exec(nitrogen); // Run the nitrogen contamination
    gSystem->Exec(pion); // Run the pion contamination
    gSystem->Exec(inelastic); // Run the inelastic contamination
    gSystem->Exec(GEN_calc); // Run the GEN_calc contamination
  }
  
  
  for(int idy = 0; idy < ndy; idy++){
    double dy = dymin + idy*dystep;

    TString file_old = cfg + "_corr.csv";
    TString file = Form(cfg + "_corr_W2_%g_%g_dy_%g.csv",-2.0,W2stable,dy);

    TString command = "cp " + DB_dir + file_old + " " + DB_dir + file;
    gSystem->Exec(command); // Create the DB file to start

    TString accidental = Form("root -b -q 'accidental_contamination.C(\"" + cfg + "\",%g,%g,%g)'",-2.0,W2stable,dy);
    TString nitrogen= Form("root -b -q 'Nitrogen_contamination.C(\"" + cfg + "\",%g,%g,%g)'",-2.0,W2stable,dy);
    TString pion = Form("root -b -q 'pion_contamination.C(\"" + cfg + "\",%g,%g,%g)'",-2.0,W2stable,dy);
    TString inelastic = Form("root -b -q 'Inelastic_contamination.C(\"" + cfg + "\",%g,%g,%g)'",-2.0,W2stable,dy);
    TString GEN_calc = Form("root -b -q 'GEN_Extraction.C(\"" + cfg + "\",%g,%g,%g)'",-2.0,W2stable,dy);

    
    gSystem->Exec(accidental); // Run the accidental contamination
    gSystem->Exec(nitrogen); // Run the nitrogen contamination
    gSystem->Exec(pion); // Run the pion contamination
    gSystem->Exec(inelastic); // Run the inelastic contamination
    gSystem->Exec(GEN_calc); // Run the GEN_calc contamination
  }
  

}
