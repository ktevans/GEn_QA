#include "../../include/gen-ana.h"
#include "../../dflay/src/JSONManager.cxx"

void organize_runs(const char *configfilename){

  JSONManager *jmgr = new JSONManager(configfilename);

  string volatile_dir = "/lustre19/expphy/volatile/halla/sbs/jeffas/GEN_root/Rootfiles";

  std::string rootfile_dir = jmgr->GetValueFromKey_str("rootfile_dir");
  std::vector<int> runnums; jmgr->GetVectorFromKey<int>("runnums",runnums);
  int nruns = jmgr->GetValueFromKey<int>("Nruns_to_ana"); // # runs to analyze
  TChain *C = new TChain("T");
  if (nruns < 1 || nruns > runnums.size()) nruns = runnums.size();
  for (int i=0; i<nruns; i++) {
    
    string runnum_string = Form("%i",runnums[i]);
    string command = "mv " + volatile_dir + "/*" + runnum_string + "* " + rootfile_dir;
    gSystem->Exec(command.c_str());
	       
  }
  
}
