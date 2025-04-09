#include <iostream>
#include <iomanip>
#include <fstream>
#include "TString.h"
#include "TObjString.h"
#include <filesystem>
namespace fs = std::filesystem;

using namespace std;

void phoebus_default_ecal_hv(int use_this_hv_crate, const char* read_this_backup_file){

  //for arg read_this_backup_file just want the name of the backup file from phoebus in aslow/EPICS/HV/backup-dir/ without the entire file path
  
  int itrip = -1250.0;
  int vmax = -1950.0;
  int rmpup = 61.0;
  int rmpdn = 61.0;

  int cells;

  if(use_this_hv_crate == 31 || use_this_hv_crate == 32 || use_this_hv_crate == 36){
    cells = 192;
  }else if(use_this_hv_crate == 33 || use_this_hv_crate == 40){
    cells = 180;
  }else if(use_this_hv_crate == 34 || use_this_hv_crate == 38){
    cells = 168;
  }else if(use_this_hv_crate == 35){
    cells = 186;
  }else if(use_this_hv_crate == 39){
    cells = 126;
  }else if(use_this_hv_crate == 41){
    cells = 72;
  }
  
  if(use_this_hv_crate < 31 || use_this_hv_crate > 41){
    cout << use_this_hv_crate << " is not a valid crate number, should be 31-41. Exiting. Note that 37 is not used, and 40 is used in it's place." << endl;
    return;
  }
 
  int ind1_arr[cells];
  int ind2_arr[cells];
  int ind3_arr[cells];
  int ind4_arr[cells];
  int ind5_arr[cells];
  int ind6_arr[cells];
  int V0Set_arr[cells];

  TString ecr_arr[cells];
  // these are the crate slot and channel numbers defined by EPICS/phoebus for each chanell. These are mapped to row col indices in HV.hvc
  
  int arr_index2 = 0;

  ifstream epicsMap(read_this_backup_file);//this should be the name of the phoebus backup file

  TString find_crate_start;
  find_crate_start.Form("# Detector: ECAL Crate %i", use_this_hv_crate);

  TString currentline_block1;

  bool found_crate = false;
  int find_end = 0;

  while( currentline_block1.ReadLine(epicsMap) ){
    
    //this breaks the loop once all channels of the crate we want have been read in
    if(found_crate == true && currentline_block1.BeginsWith( "#" ) && find_end == 1){
      break;
    }

    //find the line in the backup file which starts the ecal crate we want to read
    if(currentline_block1 == find_crate_start){
      found_crate = true;
      continue;

      //next if is to skip the first comment describing the format of the text below
      //and increment find_end so that next comment found breaks the loop over the backup file
    }else if(found_crate == true && currentline_block1.BeginsWith( "#" )){
      find_end = 1;
      continue;
    }
      
    if( found_crate == true && !currentline_block1.BeginsWith( "#" )){

      TObjArray *tokens = ( (TObjArray*) currentline_block1.Tokenize("\t") );
      int ntokens = tokens->GetEntries();
      
      
      if( ntokens == 12 ){
	TString ecr = ( (TObjString*) (*tokens)[0] )->GetString();
	
	int ind1 = ( (TObjString*) (*tokens)[1] )->GetString().Atoi();
	int ind2 = ( (TObjString*) (*tokens)[2] )->GetString().Atoi();
	int ind3 = ( (TObjString*) (*tokens)[3] )->GetString().Atoi();
	int ind4 = ( (TObjString*) (*tokens)[4] )->GetString().Atoi();
	int ind5 = ( (TObjString*) (*tokens)[5] )->GetString().Atoi();
	int ind6 = ( (TObjString*) (*tokens)[6] )->GetString().Atoi();
	int V0Set = ( (TObjString*) (*tokens)[7] )->GetString().Atoi();
	
	ind1_arr[arr_index2] = ind1;//hv crate
	ind2_arr[arr_index2] = ind2;//hv slot
	ind3_arr[arr_index2] = ind3;//hv channel
	ind4_arr[arr_index2] = ind4;//group??
	ind5_arr[arr_index2] = ind5;//row 1-69
	ind6_arr[arr_index2] = ind6;//col 1-28

	if(V0Set != 0){
	
	  //V0Set_arr[arr_index2] = V0Set;
	  V0Set_arr[arr_index2] = -1300;

	}else{
	  V0Set_arr[arr_index2] = 0;
	}
	ecr_arr[arr_index2] = ecr;//string of form ECR(row)C(col)_(crate)_(slot)_(channel)
	
	arr_index2++;
	
      }
    }    
  }

  //done parsing the phoebus backup file, now we can make the shell scripts with all of the relevant caput commands for the desired hv crate

  //fs::path(path/to/file.txt).stem() gets the filename from the directory path without the extension, .filename would get the filename with the extension

  //Form wants a char* form rather than a path, so need the .stem().c_str()

  ofstream hv_inputs(Form("HV_inputs_hv_crate%i_%s.sh", use_this_hv_crate, fs::path(read_this_backup_file).stem().c_str() ));
 
  for(int cell = 0; cell < cells; cell++){

    TString ecr_check;
    
    if(ind1_arr[cell] == use_this_hv_crate){


      TString hv_string;
      TString itrip_string;
      TString vmax_string;
      TString rmpup_string;
      TString rmpdn_string;

      //%02i means it will format an integer with two leading 0's

      hv_string.Form("caput HAHV%i:%02i:%03i:V0Set %i", ind1_arr[cell], ind2_arr[cell], ind3_arr[cell], V0Set_arr[cell]);
      
      itrip_string.Form("caput HAHV%i:%02i:%03i:I0Set %i", ind1_arr[cell], ind2_arr[cell], ind3_arr[cell], itrip);

      vmax_string.Form("caput HAHV%i:%02i:%03i:SVMax %i", ind1_arr[cell], ind2_arr[cell], ind3_arr[cell], vmax);

      rmpup_string.Form("caput HAHV%i:%02i:%03i:RUp %i", ind1_arr[cell], ind2_arr[cell], ind3_arr[cell], rmpup);

      rmpdn_string.Form("caput HAHV%i:%02i:%03i:RDWn %i", ind1_arr[cell], ind2_arr[cell], ind3_arr[cell], rmpdn);
      
      hv_inputs << hv_string << endl;
      hv_inputs << itrip_string << endl;
      hv_inputs << rmpup_string << endl;
      hv_inputs << rmpdn_string << endl;

      TString hv_check;

      hv_check.Form("caget HAHV%i:%02i:%03i:V0Setr", ind1_arr[cell], ind2_arr[cell], ind3_arr[cell]);

      hv_inputs << hv_check << " >> ecal_caget_test.txt" << endl;

    }
  }


  //hv_inputs << "if [" << hv_check <<


  ofstream power_on_crate_channels(Form("power_on_HV_channels_crate%i_%s.sh", use_this_hv_crate, fs::path(read_this_backup_file).stem().c_str() ));

  for(int cell = 0; cell < cells; cell++){
    if(ind1_arr[cell] == use_this_hv_crate){
      TString power_on_hv_string;
      int power_toggle = 0;

      if(V0Set_arr[cell] < -100 && V0Set_arr[cell] > -1860){
	power_toggle = 1;
      }

      power_on_hv_string.Form("caput HAHV%i:%02i:%03i:Pw %i", ind1_arr[cell], ind2_arr[cell], ind3_arr[cell], power_toggle);

      power_on_crate_channels << power_on_hv_string << endl;
    }
  }
  

}
