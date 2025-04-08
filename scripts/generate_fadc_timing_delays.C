#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <TFile.h>

using namespace std;


void generate_fadc_timing_delays(){

  int add_time_arr[1656];
  int arr_index1 = 0;
  int track_largest = 0;
  /*
  ifstream times("run1914/run_1914_timing_peak_multi1.txt");

  TString currentlinetime;
  while( currentlinetime.ReadLine(times) ){


      TObjArray *tokens = ( (TObjArray*) currentlinetime.Tokenize(" ") );
      int ntokens = tokens->GetEntries();

      if( ntokens == 2 ){
	int time = ( (TObjString*) (*tokens)[1] )->GetString().Atoi();

	if(time < 200){
	  time = 325;
	}
	if(time > track_largest){
	  track_largest = time;
	}
//	if(time > 325){//using 325 largest peak for now
//	  time = 325;
//	}
	add_time_arr[arr_index1] = 8;//scintillator time minus the time, since scint is max
	arr_index1++;
      }
    }
  */


  
  
  ifstream ecalfile("/adaqfs/home/a-onl/sbs/ECAL_replay/git-repo/sbs_devel/SBS-replay/DB/db_earm.ecal.dat");

  TString currentline;

  bool in_det_map = false;
  bool in_chan_map = false;

  int arr_index = 0;

  int tot_slots = 108;

  int rocid_arr[tot_slots];
  int slot_arr[tot_slots];
  int start_ch_arr[tot_slots];
  int end_ch_arr[tot_slots];

  while( currentline.ReadLine(ecalfile) ){


    if( currentline == "## crate roc id,  slot,  start_ch,  end_ch,  ref ch" ){
      in_det_map = true;
      in_chan_map = false;
      continue;
    }



    if( currentline == "earm.ecal.chanmap =" ){
      in_chan_map = true;
      in_det_map = false;
      continue;
    }




    if(in_det_map == true && in_chan_map == false){


      TObjArray *tokens = ( (TObjArray*) currentline.Tokenize("\t") );
      int ntokens = tokens->GetEntries();

      if( ntokens == 5 ){
	int rocid = ( (TObjString*) (*tokens)[0] )->GetString().Atoi();
	int slot = ( (TObjString*) (*tokens)[1] )->GetString().Atoi();
	int start_ch = ( (TObjString*) (*tokens)[2] )->GetString().Atoi();
	int end_ch = ( (TObjString*) (*tokens)[3] )->GetString().Atoi();

	rocid_arr[arr_index] = rocid - 37;//rocid is starting at 38 but want to start at 1 and go to 7
	slot_arr[arr_index] = slot;
	start_ch_arr[arr_index] = start_ch;
	end_ch_arr[arr_index] = end_ch;
	arr_index++;

      }
      
      
    }


  }

  ifstream cnffile("/adaqfs/home/sbs-onl/ecal/cfg/fadc250/ecal.cnf");

  TString currentlinecnf;

  bool start_peds = 0;
  int crate_iter;
  int slot_iter;
  int global_slot_iter = 0;
  int ecal_ch_iter = 0;
  int ecal_delay_cnf[1656];
  int ch_ind_new_line[1656];//values should be 1 - 108// 108 total slots, just want to list the global slot number for every channel
  int ch_ind_slot_val[1656];
  int ch_ind_crate_val[1656];
  
  while( currentlinecnf.ReadLine(cnffile) ){

    if(currentlinecnf == "#### PEDS only"){
      start_peds = 1;
    }
    if(start_peds == 0){
      continue;
    }else{
      for(int i = 1; i < 8; i++){
	if(currentlinecnf == Form("FADC250_CRATE sbsecal%i.jlab.org", i)){
            crate_iter = i;
          }
      }
      if( currentlinecnf(0,12) == "FADC250_SLOT"){ //look for the FADC250_SLOT # line
        TObjArray *tokens = ( (TObjArray*) currentlinecnf.Tokenize(" ") );
        int ntokens = tokens->GetEntries();
        int slot_id = ( (TObjString*) (*tokens)[1] )->GetString().Atoi();
        slot_iter = slot_id;
        //cout << currentlinecnf << endl;
      }

      int ch_amount;

      for(int i = 0; i < tot_slots; i++){
        if(rocid_arr[i] == crate_iter && slot_arr[i] == slot_iter){
          ch_amount = end_ch_arr[i] - start_ch_arr[i] + 1;
        }
      }
      TString string_start = Form("FADC250_ALLCH_DELAY");
      TString string_check;

      if( currentlinecnf(0,19) == "FADC250_ALLCH_DELAY" ){
        TObjArray *tokens = ( (TObjArray*) currentlinecnf.Tokenize(" ") );
        int ntokens = tokens->GetEntries();
	for(int i = 1; i <= ch_amount; i++){
	  int delay = ( (TObjString*) (*tokens)[i] )->GetString().Atof();
	  ecal_delay_cnf[ecal_ch_iter] = delay;
	  ch_ind_slot_val[ecal_ch_iter] = slot_iter;
	  ch_ind_crate_val[ecal_ch_iter] = crate_iter;
	  ecal_ch_iter++;
	}
	global_slot_iter++;
	cout << currentlinecnf << endl;
      }
    }
  }

  vector<int> new_delay;
  vector<int> delay_entries;
  vector<int> delay_diff;

  int ch_index = 0;
  
  ifstream time_outlier("time_outliers.txt");
  TString time_currentline;

  while( time_currentline.ReadLine(time_outlier) ){

    if(time_currentline.BeginsWith("#")){
	continue;
    }

    TObjArray *tokens = ( (TObjArray*) time_currentline.Tokenize(" ") );
    int ntokens = tokens->GetEntries();
    
    int delay_ch = ( (TObjString*) (*tokens)[0] )->GetString().Atoi();
    int delay_entries_ns = ( (TObjString*) (*tokens)[1] )->GetString().Atoi();
    int delay_diff_ns = ( (TObjString*) (*tokens)[2] )->GetString().Atoi();

    new_delay.push_back(delay_ch);
    delay_entries.push_back(delay_entries_ns);
    delay_diff.push_back(delay_diff_ns);
    
    ch_index++;
  }
  
  for(int i = 0; i < 1656; i++){
    for(int j = 0; j < new_delay.size(); j++){
      if(i == new_delay[j]){
	ecal_delay_cnf[i] = 99999;
      }
    }
  }
  
  
  //bool new_delay_slot = 1;
  for(int i = 0; i < 108; i++){
    
    if(i == 0 || rocid_arr[i] != rocid_arr[i-1]){
      cout << endl << "FADC250 Crate " << rocid_arr[i] << endl << endl;
    }
    
    //cout << rocid_arr[i] << " " << slot_arr[i] << endl;

    if(i == 0 || slot_arr[i] != slot_arr[i-1]){
      cout << "FADC250 Slot " << slot_arr[i] << endl;
    }
    
    bool new_delay_slot = 1;
    for(int j = 0; j <= end_ch_arr[i]; j++){

      int k_good;
      bool yes = 0;

      for(int k = 0; k < 1656; k++){
	if(ch_ind_crate_val[k] == rocid_arr[i] && ch_ind_slot_val[k] == slot_arr[i]){
	  k_good = k;
	  break;
	}
      }
      //cout << ch_ind_crate_val[k_good] << " " << ch_ind_slot_val[k_good] << " " << k_good << endl;
      //cout << k_good+j << endl;
      
      if(new_delay_slot == 1){
	//cout << "Crate: " << ch_ind_crate_val[k_good+j] << ", Slot: " << ch_ind_slot_val[k_good+j] << "    FADC250_ALLCH_DELAY ";
	cout << "FADC250_ALLCH_DELAY ";
      }
      new_delay_slot = 0;
      
      cout << ecal_delay_cnf[k_good+j] << " ";
    }

    if( end_ch_arr[i] < 15 ){
      for(int z = 0; z < 15 - end_ch_arr[i]; z++){
	cout << 0 << " ";
      }
    }
    
    cout << endl;
  }


  /*
  
  TString new_delay_line;

  for(int i = 0; i < 109; i++){
    bool new_delay_slot = 1;
    for(int j = 0; j < 1656; j++){
      
      //for(int k = 1; k < (new_delay.size()); k++){
	//if(ch_ind_crate_val[new_delay[i]] == ch_ind_crate_val[new_delay[i+k]] && ch_ind_slot_val[new_delay[i+k]] == ch_ind_slot_val[new_delay[i+k]]){
	  
	//}
      //}
      
      if(new_delay[i] == j){
	if(new_delay_slot == 1){
	  cout << "Crate: " << ch_ind_crate_val[j] << ", Slot: " << ch_ind_slot_val[j] << "    FADC250_ALLCH_DELAY ";
	}
	new_delay_slot = 0;

	cout << "new ";
      }else if(ch_ind_slot_val[j] == ch_ind_slot_val[new_delay[i]] && ch_ind_crate_val[j] == ch_ind_crate_val[new_delay[i]]){
	if(new_delay_slot == 1){
	  cout << "Crate: " << ch_ind_crate_val[j] << ", Slot: " << ch_ind_slot_val[j] << "    FADC250_ALLCH_DELAY ";
	}
	new_delay_slot = 0;

	cout << ecal_delay_cnf[j] << " ";

      }
    }
    cout << endl;
  }
  */
  

  /*
   for(int npmt=0; npmt<1656;npmt++){
  new_delay[npmt]= add_time_arr[npmt]+ecal_delay_cnf[npmt];

   }
  ofstream fadc_list;
  fadc_list.open("fadc_list_test.txt");

  fadc_list << "#pmt, crate, slot, ch, delay" << endl;

  int pmtindx = 0;

  for(int all_slots = 0; all_slots < tot_slots; all_slots++){

    for(int ch = start_ch_arr[all_slots]; ch <= end_ch_arr[all_slots]; ch++ ){
      fadc_list << pmtindx << " " << rocid_arr[all_slots] << " " << slot_arr[all_slots] << " " << ch << " " << new_delay[pmtindx] << endl;
	pmtindx++;
    }
    
  }
  */

}
