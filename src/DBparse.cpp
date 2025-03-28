#include "../include/DBparse.h"
#include "../include/Analysis.h"

namespace DBparse {

  void DB_load(DBInfo &request){

    for(int ivar = 0; ivar < request.var_req.size(); ivar++){
      DBrequest var_info = request.var_req[ivar];
      bool found_var = false;

      TString DIR = DB_dir;

      auto it = DBFileMap.find(var_info.var_names);
      TString file = it->second;
      if(it->first == "Asymmetry Correction"){
	file = Form(request.cfg + "_" + file + "_W2_%g_%g_dy_%g.csv",request.W2min,request.W2max,request.dymax);
	DIR = DB_corr_dir;
      }
      
      fstream file_csv; file_csv.open(DIR + file);
      
      // Here we check that the DB files are there
      if (it == DBFileMap.end()) { // This is not in the list of DB files
        if(!var_info.mandatory){
	  cout<<"DB variable NOT found but NOT mandatory: "<<var_info.var_names<<endl;
	  continue;
	}
	if(var_info.mandatory){
	  cout<<"!!!!WARNING, ERROR, DBparse::DB_load: DB variable NOT found and IS mandatory: "<<var_info.var_names<<endl;
	  exit(0);
	}
      }
      else if (!file_csv.is_open()) {// This is in the list but the file is not there
	if(!var_info.mandatory){
	  cout<<"DB file NOT found but NOT mandatory: "<<file<<endl;
	  continue;
	}
	if(var_info.mandatory){
	  cout<<"!!!!WARNING, ERROR, DBparse::DB_load: File "<<file<<" is not found"<<endl;
	  exit(0);
	}
      }
      else {// The file is found and everything is fine
	cout<<"DB file found: "<<file<<endl;
      }
      
      string line;
      int iline = 0;
      
      while (getline(file_csv, line)) {
	// Create a stringstream to parse the line
	stringstream ss(line);
	string cell;
	iline++;
	
	vector<string> val;

	// Split the line into cells using a comma as a delimiter
	while (getline(ss, cell, ',')) {
	  val.push_back(cell);  // Put one line into vectros
	}
	if(it->first == "He3 Polarization"){// He3 polarization DB
	  if(iline == 1) continue;
	  request.He3Pol.insert(std::make_pair(Utilities::SetTime(val[0]),stod(val[1])));
	}// end He3 Pol
	else if(it->first == "Beam Polarization"){// Beam polarization DB
	  if(iline == 1) continue;	 
	  vector<TDatime*> beam_time;
	  vector<double> beam_val;
	  
	  for(int irow = 0; irow < val.size(); irow++){ 
	    if(irow < 2) { // Time numbers
	      TDatime *time = new TDatime(val[irow].c_str());
	      beam_time.push_back(time);
	    }
	    if(irow >= 2) { // Values
	      beam_val.push_back(stod(val[irow]));
	    }
	  }
	  
	  request.BeamPolTime.push_back(beam_time);
	  request.BeamPolValue.push_back(beam_val); 
	} // End Beam Pol
	else if(it->first == "Helicity Quality"){// Helicity quality DB
	  if(iline == 1) continue;	  
	  request.GoodHel.insert(std::make_pair(stoi(val[0]),stod(val[1])));
	}
	else if(it->first == "Moller Quality"){// Moller quality DB
	  if(iline == 1) continue;	 
	  request.GoodMoller.insert(std::make_pair(stoi(val[0]),stod(val[1])));
	}
	else if(it->first == "Field Measurement"){// Field Measurement DB
	  if(iline == 1) continue;	  
	  if(request.cfg == val[0]){
	    field_hor_ang = stod(val[1]);
	    field_vert_ang = stod(val[2]);
	  }
	}// End Field Measurement
	else if(it->first == "Field Measurement"){// Field Measurement DB
	  if(iline == 1) continue;	  
	  if(request.cfg == val[0]){
	    field_hor_ang = stod(val[1]);
	    field_vert_ang = stod(val[2]);
	  }
	}// End Field Measurement
	else if(it->first == "Asymmetry Correction"){// Asym Correction DB
	  if(val.size() == 0) continue;

	  if(val[0] == "pion"){// Pion Corrections
	    request.PionAsymmetry = stod(val[1]);
	    request.PionAsymmetryErr = stod(val[2]);
	    request.PionFraction = stod(val[3]);
	    request.PionFractionErr = stod(val[4]);
	  } 
	  else if(val[0] == "Inelastic"){// Inelastic Corrections
	    request.InelasticAsymmetry = stod(val[1]);
	    request.InelasticAsymmetryErr = stod(val[2]);
	    request.InelasticFraction = stod(val[3]);
	    request.InelasticFractionErr = stod(val[4]);
	  } 
	  else if(val[0] == "Nitrogen"){// Nitrogen Corrections
	    request.NitrogenFraction = stod(val[1]);
	    request.NitrogenFractionErr = stod(val[2]);
	  } 
	  else if(val[0] == "time accidentals"){// Accidental Corrections
	    request.AccidentalAsymmetry = stod(val[1]);
	    request.AccidentalAsymmetryErr = stod(val[2]);
	    request.AccidentalFraction = stod(val[3]);
	    request.AccidentalFractionErr = stod(val[4]);
	  } 

	}// End Asym Correction 
      }
      
    }
  }


  void DB_SetCorrections(DBInfo &DBInfo){

    TString file = Form(DBInfo.cfg + "_corr_W2_%g_%g_dy_%g.csv",DBInfo.W2min,DBInfo.W2max,DBInfo.dymax);

    ofstream file_csv; file_csv.open(DB_corr_dir + file);

    file_csv<<"#time accidentals,Asymmetry, Asymmetry Err, Fraction, Fraction Err"<<endl;
    file_csv<<Form("time accidentals,%.4f,%.4f,%.4f,%.4f",DBInfo.AccidentalAsymmetry,DBInfo.AccidentalAsymmetryErr,DBInfo.AccidentalFraction,DBInfo.AccidentalFractionErr)<<endl;
    file_csv<<endl;

    file_csv<<"#Nitrogen, Fraction, Fraction Err"<<endl;
    file_csv<<Form("Nitrogen,%.4f,%.4f",DBInfo.NitrogenFraction,DBInfo.NitrogenFractionErr)<<endl;
    file_csv<<endl;
 
    file_csv<<"#pion,Asymmetry, Asymmetry Err, Fraction, Fraction Err"<<endl;
    file_csv<<Form("pion,%.4f,%.4f,%.4f,%.4f",DBInfo.PionAsymmetry,DBInfo.PionAsymmetryErr,DBInfo.PionFraction,DBInfo.PionFractionErr)<<endl;
    file_csv<<endl;
 
    file_csv<<"#Inelastic,Asymmetry, Asymmetry Err, Fraction, Fraction Err"<<endl;
    file_csv<<Form("Inelastic,%.4f,%.4f,%.4f,%.4f",DBInfo.InelasticAsymmetry,DBInfo.InelasticAsymmetryErr,DBInfo.InelasticFraction,DBInfo.InelasticFractionErr)<<endl;
  }
  
}
