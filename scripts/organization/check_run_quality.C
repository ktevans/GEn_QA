//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//   Created by Sean Jeffas, sj9ry@virginia.edu
//   Last Modified February 9, 2024
//
//
//   The purpose of this script check runs for good helicity
//   readout and create a csv file with the results for 
//   easier recovery.
//
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

#include "../../include/gen-ana.h"
#include "../../dflay/src/JSONManager.cxx"


void check_run_quality(){

  const int ncfg = 4;
  TString cfgs[ncfg] = {"GEN2","GEN3","GEN4","GEN4b"};

  std::vector<int> runnums;
  TChain *T = new TChain("Tout");

  for(int icfg=0; icfg < ncfg; icfg++){
    
    TString jmgr_file = "../../config/" + cfgs[icfg] + "_He3.cfg";
    JSONManager *jmgr = new JSONManager(jmgr_file);

    if(icfg == 0)
      jmgr->GetVectorFromKey<int>("runnums",runnums);
    else{
      std::vector<int> runnums_temp; 
      jmgr->GetVectorFromKey<int>("runnums",runnums_temp);
      runnums.insert( runnums.end(), runnums_temp.begin(), runnums_temp.end() );
    }


    T->Add("../outfiles/QE_data_" + cfgs[icfg] + "_sbs100p_nucleon_np_model2.root");

  }

  T->SetBranchStatus("*",0);

  // This is the variables we need from the analyzed root file
  int runnum;   setrootvar::setbranch(T,"runnum","",&runnum);
  int helicity;   setrootvar::setbranch(T,"helicity","",&helicity);

  TString outfilename = "../../DB/Helicity_quality.csv";
  ofstream outfile;
  outfile.open (outfilename);

  outfile<<"Run Number,Good Helicity"<<endl;

  for(int irun = 0; irun < runnums.size(); irun++){
    cout << 1.0*irun / runnums.size() << "% \r";
    std::cout.flush();
    
    bool is_good = true;
    TH1F *h = new TH1F("h","",3,-1.5,1.5);
    T->Draw("helicity>>h",Form("runnum == %i",runnums[irun]));
      
    if(h->GetBinContent(2) / h->GetEntries() > 0.03 || T->GetEntries(Form("runnum == %i",runnums[irun])) == 0)
      is_good = false;
      
    outfile<<runnums[irun]<<","<<is_good<<endl;
    h->Delete();
  }

}
