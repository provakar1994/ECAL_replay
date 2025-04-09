#include <iostream>
#include <iomanip>
#include <fstream>
#include "TString.h"
#include "TObjString.h"
#include <filesystem>

using namespace std;

void bad_ch_check(const *char caget_list_textfile){

  int cells = 1656;

  int arr_index = 0;

  int row_arr[cells];
  int col_arr[cells];
  double xpos_arr[cells];
  double ypos_arr[cells];
  
  ifstream ecalfile("maps/ECAL_r_c_x_y_cpr.csv");

  TString currentline;

  int check = 0;
  
  while( currentline.ReadLine(ecalfile) ){

    //cout << currentline(0,2) << endl;

    if( check != 0  ){ //new cell information:
      TObjArray *tokens = ( (TObjArray*) currentline.Tokenize(",") );
      int ntokens = tokens->GetEntries();

      //cout << currentline << " " << ntokens << endl;

      if( ntokens == 6 ){
	int cell = ( (TObjString*) (*tokens)[0] )->GetString().Atoi();
	int row = ( (TObjString*) (*tokens)[1] )->GetString().Atoi();
	int col = ( (TObjString*) (*tokens)[2] )->GetString().Atoi();
	double pos_hori = ( (TObjString*) (*tokens)[3] )->GetString().Atof();
	double pos_vert = ( (TObjString*) (*tokens)[4] )->GetString().Atof();

	row_arr[arr_index] = row;
	col_arr[arr_index] = col;
	xpos_arr[arr_index] = pos_hori;
	ypos_arr[arr_index] = pos_vert;
	arr_index++;

	//cout << cell << " " << row << " " << col << " " << pos_hori << " " << pos_vert << endl;

      }
    }
    check = 1;
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

  TString read_this_backup_file;
  read_this_backup_file = Form("/adaqfs/home/aslow/EPICS/HV/earm_ecal/HV.hvc");

  ifstream epicsMap(read_this_backup_file);//this should be the name of the phoebus backup file

  TString currentline_block1;

  while( currentline_block1.ReadLine(epicsMap) ){

    if( !currentline_block1.BeginsWith( "#" )){
      
      TObjArray *tokens = ( (TObjArray*) currentline_block1.Tokenize(" ") );
      int ntokens = tokens->GetEntries();
      
      
      if( ntokens == 7 ){
	TString ecr = ( (TObjString*) (*tokens)[0] )->GetString();
	
	int ind1 = ( (TObjString*) (*tokens)[1] )->GetString().Atoi();
	int ind2 = ( (TObjString*) (*tokens)[2] )->GetString().Atoi();
	int ind3 = ( (TObjString*) (*tokens)[3] )->GetString().Atoi();
	int ind4 = ( (TObjString*) (*tokens)[4] )->GetString().Atoi();
	int ind5 = ( (TObjString*) (*tokens)[5] )->GetString().Atoi();
	int ind6 = ( (TObjString*) (*tokens)[6] )->GetString().Atoi();
	
	ind1_arr[arr_index2] = ind1;//hv crate
	ind2_arr[arr_index2] = ind2;//hv slot
	ind3_arr[arr_index2] = ind3;//hv channel
	ind4_arr[arr_index2] = ind4;//group??
	ind5_arr[arr_index2] = ind5;//row 1-69
	ind6_arr[arr_index2] = ind6;//col 1-28
	
	ecr_arr[arr_index2] = ecr;//string of form ECR(row)C(col)_(crate)_(slot)_(channel)
	
	arr_index2++;
	
      }
    }
  }
    
  //ifstream caget_list("/adaqfs/home/aslow/nhunt/apr8_1200V_all_beam_comm_hv_control/test.txt");
  ifstream caget_list(caget_list_textfile);
  
  TString currentline_caget;
  int caget_ind = 0;

  vector<string> store_hv_epics_if_0;
  vector<int> hv_crate_arr_caget;
  vector<int> hv_slot_arr_caget;
  vector<int> hv_ch_arr_caget;

  while( currentline_caget.ReadLine(caget_list) ){

    if( !currentline_caget.BeginsWith( "#" )){

      TObjArray *tokens = ( (TObjArray*) currentline_caget.Tokenize(" ") );
      int ntokens = tokens->GetEntries();
      
      //cout << currentline_caget << " " << ntokens << endl;
      
      if( ntokens ==2){
	int ind1 = ( (TObjString*) (*tokens)[1] )->GetString().Atoi();
	TString ind0 = ( (TObjString*) (*tokens)[0] )->GetString();
      //TString zero = "0";
	if(ind1 == 0){
	
	//cout << currentline_caget << endl;
	
	//store_hv_epics_if_0.push_back(int (currentline_caget(0,13)));
	//hv_crate_arr_caget.push_back(stoi(currentline_caget(4,6)));
	//hv_slot_arr_caget.push_back(stoi(currentline_caget(7,9)));
	//hv_ch_arr_caget.push_back(stoi(currentline_caget(11,13)));

	hv_crate_arr_caget.push_back(stoi(ind0(4,6)));
	hv_slot_arr_caget.push_back(stoi(ind0(7,9)));
	hv_ch_arr_caget.push_back(stoi(ind0(11,13)));
	}
      }
    }
  }
  
  ifstream bad_ch_list_JY("maps/bad_ecal_channels_04_08_25.csv");

  TString currentline_bad;
  vector<int> bad_chan_vector;

  while( currentline_bad.ReadLine(bad_ch_list_JY) ){
    TObjArray *tokens = ( (TObjArray*) currentline_bad.Tokenize(" ") );
    int bad_ch = ( (TObjString*) (*tokens)[0] )->GetString().Atoi();
    bad_chan_vector.push_back(bad_ch);
    //cout << bad_ch << endl;
  }

  for(int i = 0; i < hv_crate_arr_caget.size(); i++){

    int caget_ch;

    for(int ch = 0; ch < 1656; ch++){
      if( ind1_arr[ch] == hv_crate_arr_caget[i] && ind2_arr[ch] == hv_slot_arr_caget[i] && ind3_arr[ch] == hv_ch_arr_caget[i] ){
	caget_ch = ch;
      }
    }
    
    int check = 0;
    for(int j = 0; j < bad_chan_vector.size(); j++){
      //cout << hv_crate_arr_caget[i] << " " << hv_slot_arr_caget[i] << " " << hv_ch_arr_caget[i] << endl;
      if(caget_ch == bad_chan_vector[j]){
	check = 1;
      }
      
    }
    if(check == 0){
      cout << "ch: " << caget_ch << " is not on the bad channel list" << endl;
    }
  }


  
  for(int i = 0; i < bad_chan_vector.size(); i++){

    int check = 0;
    for(int j = 0; j < hv_crate_arr_caget.size(); j++){

      int caget_ch;

      for(int ch = 0; ch < 1656; ch++){
	if( ind1_arr[ch] == hv_crate_arr_caget[j] && ind2_arr[ch] == hv_slot_arr_caget[j] && ind3_arr[ch] == hv_ch_arr_caget[j] ){
	  caget_ch = ch;
	}
      }

      if(caget_ch == bad_chan_vector[i]){
	check = 1;
      }
      
    }
    if(check == 0){
      cout << "ch: " << bad_chan_vector[i] << " is not set to 0V according to the caget list" << endl;
    }
  }
    
}

/*











  
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
      
  }


}
*/
