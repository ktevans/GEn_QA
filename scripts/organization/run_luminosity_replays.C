


void run_luminosity_replays(){


  ifstream infile("../../../GEM_analysis/scripts/GEM_luminosity_good_events.txt");

  
  TString first, second;
  int run = -1;
  int first_event, last_event;
  string line;
  while(getline(infile, line)){ 
    
    if(line.empty() || line[0] == '#') continue;

    std::istringstream iss(line);
    iss >> first >> second;

    if(first == "run") run = stoi(second.Data());
    else{
      first_event = stoi(first.Data());
      last_event = stoi(second.Data());
      
      gSystem->Exec(Form("run_luminosity_replay.sh %i %i %i",run,first_event,last_event));
    }

  }





}
