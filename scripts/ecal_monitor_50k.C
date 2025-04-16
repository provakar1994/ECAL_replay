#include <TTree.h>
#include <TMatrixD.h>
#include <TROOT.h>
#include <TMatrixDBase.h>
#include <TMath.h>
#include <TChain.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TF1.h>
#include <cstdlib>
#include "TColor.h"
#include "TAttFill.h"
#include "stdlib.h"
#include "TStreamerInfo.h"
#include "TStreamerElement.h"
#include <vector>
#include <TChain.h>
#include "TCanvas.h"
#include "TGraphErrors.h"
#include <iostream>
#include <fstream>
#include "TSystem.h"
#include "TFile.h"
#include "TString.h"

void ecal_monitor_50k(int second_file){//compare second file to first file

  int first_file = 2379;//this should be the golden run

  //want first file to "Golden Run", so just need to input replayed run that you want to compare it to, run2159 is just the most recent cosmic run with decent timing, so we will want to change this once we have beam
  
  //want to look at two different runs and make comparisons of the peds, energy spectra, and dead channels. Update the root file to keep track of these. Make 2d heatmap of the differences between the two

  bool fit_arg = 0;

  //cout << "Type 0 for Gaus fit, type 1 for expo fit." << endl;
  //cin >> fit_arg;








  //Read ecal fadc info from database file in order to read in pedestals from ecal.cnf to get expected peds

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
  int ecal_ch_iter = 0;
  double ecal_ped_cnf[1656];
  
  while( currentlinecnf.ReadLine(cnffile) ){

    //cout << currentlinecnf << endl;

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
	//cout << ntokens << endl;
      }
      
      int ch_amount;
      
      for(int i = 0; i < tot_slots; i++){
	if(rocid_arr[i] == crate_iter && slot_arr[i] == slot_iter){
	  ch_amount = end_ch_arr[i] - start_ch_arr[i] + 1;
	}
      }
      

      TString string_start = Form("FADC250_ALLCH_PED");
      TString string_check;
      
      if( currentlinecnf(0,17) == "FADC250_ALLCH_PED" ){
	TObjArray *tokens = ( (TObjArray*) currentlinecnf.Tokenize("  ") );
	int ntokens = tokens->GetEntries();
	//cout << ntokens << endl;
	  for(int i = 1; i <= ch_amount; i++){
	    double ped = ( (TObjString*) (*tokens)[i] )->GetString().Atof();
	    ecal_ped_cnf[ecal_ch_iter] = ped/4.096;
	    //cout << ped << endl;
	    ecal_ch_iter++;
	  }
      }
      
    }

  }











  
  int arr_index1 = 0;
  int cells = 1656;

  double ex[cells], mean_ECAL[cells], sigma_ECAL[cells], dmean_ECAL[cells], dsigma_ECAL[cells];

  double cells_ECAL[cells];
  
  int row_arr[cells];
  int col_arr[cells];
  double xpos_arr[cells];
  double ypos_arr[cells];

  double mc_ep[2000];
  int ecal_adcID[1656];
  
  ifstream ecalfile1("maps/ECAL_r_c_x_y_cpr.csv");

  TString currentline1;

  int check = 0;

  while( currentline1.ReadLine(ecalfile1) ){

    if( check != 0 ) { //new cell information:
      TObjArray *tokens = ( (TObjArray*) currentline1.Tokenize(",") );
      int ntokens = tokens->GetEntries();

      if( ntokens == 7 ){
	int cell = ( (TObjString*) (*tokens)[0] )->GetString().Atoi();
	int row = ( (TObjString*) (*tokens)[1] )->GetString().Atoi();
	int col = ( (TObjString*) (*tokens)[2] )->GetString().Atoi();
	double pos_hori = ( (TObjString*) (*tokens)[3] )->GetString().Atof();
	double pos_vert = ( (TObjString*) (*tokens)[4] )->GetString().Atof();

	//index 0 corresponds to row0 col0
	row_arr[arr_index1] = row;
	col_arr[arr_index1] = col;
	xpos_arr[arr_index1] = pos_hori;
	ypos_arr[arr_index1] = pos_vert;
	arr_index1++;

      }
    }
    check = 1;
  }

  TH1F *htemp;

  int list_of_changed_ch[1656];

  //everything raw before any clustering is done
  Int_t ndata_adc1;
  Double_t adc_amp_mv1[1656];
  Double_t adc_ch_id1[1656];
  Double_t adc_int_pC1[1656];
  Double_t adc_row1[1656];
  Double_t adc_col1[1656];
  Double_t adc_ped1[1656];
  Double_t adc_xpos1[1656];
  Double_t adc_ypos1[1656];
  Double_t adc_time1[1656];

  //blocks in primary clusters
  Int_t ndata_clus_blk1;
  Double_t clus_blk_ch_id1[1656];
  Double_t clus_blk_e1[1656];
  Double_t clus_blk_row1[1656];
  Double_t clus_blk_col1[1656];
  Double_t clus_blk_x1[1656];
  Double_t clus_blk_y1[1656];

  //individual clusters, includes all clusters
  Int_t ndata_clus1;
  Double_t ndata_clus_arr1[1656];
  Double_t clus_ch_id1[1656];//this is the cluster seed ch id
  Double_t clus_e1[1656];
  Double_t clus_row1[1656];//also corresponds to the cluster seed (I think)
  Double_t clus_col1[1656];
  Double_t clus_x1[1656];
  Double_t clus_y1[1656];



  //everything raw before any clustering is done
  Int_t ndata_adc2;
  Double_t adc_amp_mv2[1656];
  Double_t adc_ch_id2[1656];
  Double_t adc_int_pC2[1656];
  Double_t adc_row2[1656];
  Double_t adc_col2[1656];
  Double_t adc_ped2[1656];
  Double_t adc_xpos2[1656];
  Double_t adc_ypos2[1656];
  Double_t adc_time2[1656];
  /*
  //blocks in primary clusters
  Int_t ndata_clus_blk2;
  Double_t clus_blk_ch_id2[1656];
  Double_t clus_blk_e2[1656];
  Double_t clus_blk_row2[1656];
  Double_t clus_blk_col2[1656];
  Double_t clus_blk_x2[1656];
  Double_t clus_blk_y2[1656];

  //individual clusters, includes all clusters
  Int_t ndata_clus2;
  Double_t ndata_clus_arr2[1656];
  Double_t clus_ch_id2[1656];//this is the cluster seed ch id
  Double_t clus_e2[1656];
  Double_t clus_row2[1656];//also corresponds to the cluster seed (I think)
  Double_t clus_col2[1656];
  Double_t clus_x2[1656];
  Double_t clus_y2[1656];
  */


  TH1F* adc_peds1_per_block[1656];
  TH1F* adc_peds2_per_block[1656];

  TH1F* adc_int1_per_block[1656];
  TH1F* adc_int2_per_block[1656];

  TH1F* adc_time1_per_block[1656];
  TH1F* adc_time2_per_block[1656];  

  TH2D *peds_diff_vs_xy = new TH2D("peds_diff_golden_run_vs_row_col", "", 27, 0.5, 27.5, 69, 0.5, 69.5);
  TH2D *int_diff_vs_xy = new TH2D("int_diff_vs_xy", "", 27, 0.5, 27.5, 69, 0.5, 69.5);
  TH2D *ch_diff_vs_xy = new TH2D("ch_diff__vs_xy", "", 27, 0.5, 27.5, 69, 0.5, 69.5);
  TH2D *peds_width_diff_vs_xy = new TH2D("peds_width_diff_golden_run_vs_row_col", "", 27, 0.5, 27.5, 69, 0.5, 69.5);
  TH2D *time_diff_vs_xy = new TH2D("time_diff_golden_run_vs_row_col", "", 27, 0.5, 27.5, 69, 0.5, 69.5);
  TH2D *time_diff_from_expect_vs_xy = new TH2D("time_diff_from_expect_vs_row_col", "", 27, 0.5, 27.5, 69, 0.5, 69.5);

  TH2D *peds_diff_vs_ch = new TH2D("peds_diff_golden_run_vs_ch", "", 1656, 0.5, 1656.5, 100, 0, 5);
  TH2D *peds_width_diff_vs_ch = new TH2D("peds_width_diff_golden_run_vs_ch", "", 1656, 0.5, 1656.5, 200, 0, 1.0);
  TH2D *int_diff_vs_ch = new TH2D("int_diff_golden_run_vs_ch", "", 1656, 0.5, 1656.5, 200, 0, 20);
  TH2D *ch_diff_vs_ch = new TH2D("ch_diff_golden_run_vs_ch", "", 1656, 0.5, 1656.5, 10, 0, 10);
  TH2D *peds_cnf_diff_vs_ch = new TH2D("peds_diff_cnf_vs_ch", "peds_diff_cnf_vs_ch", 1656, 0.5, 1656.5, 250, 0, 5);
  TH2D *time_diff_vs_ch = new TH2D("time_diff_golden_run_vs_ch", "", 1656, 0.5, 1656.5, 200, -50, 50);
  TH2D *time_diff_from_expect_vs_ch = new TH2D("time_diff_from_expect_vs_ch", "", 1656, 0.5, 1656.5, 100, -50, 50);

  TH2D *peds_vs_xy = new TH2D("peds_vs_row_col", "peds_vs_row_col", 27, 0.5, 27.5, 69, 0.5, 69.5);
  TH2D *int_vs_xy = new TH2D("int_vs_xy", "int_vs_row_col", 27, 0.5, 27.5, 69, 0.5, 69.5);
  TH2D *ch_vs_xy = new TH2D("ch_vs_xy", "ch_vs_xy", 27, 0.5, 27.5, 69, 0.5, 69.5);
  TH2D *peds_width_vs_xy = new TH2D("peds_width_vs_row_col", "peds_width_vs_row_col", 27, 0.5, 27.5, 69, 0.5, 69.5);
  TH2D *time_vs_xy = new TH2D("time_vs_row_col", "time_vs_row_col", 27, 0.5, 27.5, 69, 0.5, 69.5);
  
  TH2D *peds_vs_ch = new TH2D("peds_vs_ch", "peds_vs_ch", 1656, 0.5, 1656.5, 400, 60, 140);
  TH2D *peds_width_vs_ch = new TH2D("peds_width_vs_ch", "peds_width_vs_ch", 1656, 0.5, 1656.5, 200, 0, 4.0);
  TH2D *int_vs_ch = new TH2D("int_vs_ch", "int_vs_ch", 1656, 0.5, 1656.5, 200, 0, 20);
  TH2D *time_vs_ch = new TH2D("time_vs_ch", "time_vs_ch", 1656, 0.5, 1656.5, 200, 150, 250);
  TH2D *ch_vs_ch = new TH2D("ch_vs_ch", "", 1656, 0.5, 1656.5, 10, 0, 10);
  
  
  for(int i = 0; i < 1656; i++){
    adc_peds1_per_block[i] = new TH1F("", "", 300, 50, 150);
    adc_peds2_per_block[i] = new TH1F("", "", 300, 50, 150);

    adc_int1_per_block[i] = new TH1F("", "", 100, 0, 80);
    adc_int2_per_block[i] = new TH1F("", "", 100, 0, 80);
    
    adc_time1_per_block[i] = new TH1F("", "", 100, 150, 250);
    adc_time2_per_block[i] = new TH1F("", "", 100, 150, 250);
  }
  
  TH1F *adc_peds1 = new TH1F("adc_peds1", "", 300, 50, 150);
  TH1F *adc_peds2 = new TH1F("adc_peds2", "", 300, 50, 150);

  TH1F *adc_energy1 = new TH1F("adc_energy1", "", 100, 130, 160);
  TH1F *adc_energy2 = new TH1F("adc_energy2", "", 100, 130, 160);

  TH1F *adc_hits1 = new TH1F("adc_hits1", "", 100, 130, 160);
  TH1F *adc_hits2 = new TH1F("adc_hits2", "", 100, 130, 160);


























     //
     TH2D* h_vnclus_enclus;
     h_vnclus_enclus = new TH2D("h_vnclus_enclus"," VTP DETID =12; VTP nclus ; ECAL nclus",20,0,20,20,0,20);
       HList.Add(h_vnclus_enclus);
    TH2D* h_hvnclus_hnclus;
     h_hvnclus_hnclus = new TH2D("h_hvnclus_hnclus"," hVTP ; hVTP nclus ; HCAL nclus",20,0,20,20,0,20);
       HList.Add(h_hvnclus_hnclus);
    TH2D* h_hvtp_x_y;
     h_hvtp_x_y = new TH2D("h_hvtp_x_y"," hVTP ; Col ; Row ",20,0,20,20,0,20);
       HList.Add(h_hvtp_x_y);
    TH2D* h_hvtp_x_y_time;
     h_hvtp_x_y_time = new TH2D("h_hvtp_x_y_time"," hVTP ; Col ; Row ",20,0,20,20,0,20);
       HList.Add(h_hvtp_x_y_time);
     TH2D* h_vnclus_enclus_detidzero;
     h_vnclus_enclus_detidzero = new TH2D("h_vnclus_enclus_detidzero","VTp DetID = 0 ; VTP nclus ; ECAL nclus",20,0,20,20,0,20);
       HList.Add(h_vnclus_enclus_detidzero);
     TH2D* h_ecal_id_e;
     h_ecal_id_e = new TH2D("h_ecal_id_e"," ; ID ; Energy",1656,0,1656,100,0,4.);
       HList.Add( h_ecal_id_e);
     TH2D* h_ecal_row_col;
      h_ecal_row_col = new TH2D("h_ecal_row_col"," ; Col ; Row",27,1,28,69,1,70);
       HList.Add( h_ecal_row_col);
     TH2D* h_vtp_row_col;
      h_vtp_row_col = new TH2D("h_vtp_row_col"," ; Col ; Row",27,1,28,69,1,70);
       HList.Add( h_vtp_row_col);
    TH2D* h_ecal_row_col_e;
      h_ecal_row_col_e = new TH2D("h_ecal_row_col_e","Energy ; Col ; Row",27,1,28,69,1,70);
       HList.Add( h_ecal_row_col_e);
    TH2D* h_ecal_row_col_t;
      h_ecal_row_col_t = new TH2D("h_ecal_row_col_t","Time ; Col ; Row",27,1,28,69,1,70);
       HList.Add( h_ecal_row_col_t);
     TH2D* h_vtpEclus_ecalEclus;
     h_vtpEclus_ecalEclus = new TH2D("h_vtpEclus_ecalEclus"," ; VTP E clus ; ECAL E clus",200,0,4000,200,0,4);
       HList.Add(h_vtpEclus_ecalEclus);
     TH2D* h_vtpEclus_ecalEclus_detidzero;
     h_vtpEclus_ecalEclus_detidzero = new TH2D("h_vtpEclus_ecalEclus_detidzero","VTP DEtid = 0 ; VTP E clus ; ECAL E clus",200,0,4000,200,0,4);
       HList.Add(h_vtpEclus_ecalEclus_detidzero);
      TH1D* h_vtpEclus;
       h_vtpEclus = new TH1D("h_vtpEclus"," ; VTP E clus ",200,0,4000);
       HList.Add(h_vtpEclus);
     TH1D* h_vtpSizeclus;
       h_vtpSizeclus = new TH1D("h_vtpSizeclus"," ; VTP Size clus ",10,0,10);
       HList.Add(h_vtpSizeclus);
      TH1D* h_vtpTclus;
       h_vtpTclus = new TH1D("h_vtpTclus"," ; VTP Time clus ",100,0,100);
       HList.Add(h_vtpTclus);
      TH1D* h_hvtpEclus;
       h_hvtpEclus = new TH1D("h_hvtpEclus"," ; HVTP E clus ",200,0,4000);
       HList.Add(h_hvtpEclus);
     TH1D* h_hvtpSizeclus;
       h_hvtpSizeclus = new TH1D("h_hvtpSizeclus"," ; HVTP Size clus ",10,0,10);
       HList.Add(h_hvtpSizeclus);
      TH1D* h_hvtpTclus;
       h_hvtpTclus = new TH1D("h_hvtpTclus"," ; HVTP Time clus ",100,0,100);
       HList.Add(h_hvtpTclus);
    TH1D* h_ecalEclus;
     h_ecalEclus = new TH1D("h_ecalEclus"," ; ECAL E clus",200,0,4);
       HList.Add(h_ecalEclus);
     TH1D* h_ecalSizeclus;
     h_ecalSizeclus = new TH1D("h_ecalSizeclus"," ; ECAL Size clus",10,0,10);
       HList.Add(h_ecalSizeclus);
      TH2D* h_ecalSize_ecalEclus;
     h_ecalSize_ecalEclus = new TH2D("h_ecalSize_ecalEclus"," ; ECal Size clus ; ECAL E clus",20,0,20,200,0,4);
       HList.Add(h_ecalSize_ecalEclus);
     TH1D* h_ecalTclus;
      h_ecalTclus = new TH1D("h_ecalTclus"," ; ECAL Time clus",400,0,400);
       HList.Add(h_ecalTclus);
     TH2D* h_hcalTime_ecaltime;
     h_hcalTime_ecaltime = new TH2D("h_hcalTime_ecaltime"," ; HCal time; ECAL Time",200,0,400,400,0,400);
       HList.Add(h_hcalTime_ecaltime);
     TH1D* h_hcalTime;
     h_hcalTime = new TH1D("h_hcalTime"," ; HCal time",200,0,400);
       HList.Add(h_hcalTime);
     TH1D* h_hcalE;
     h_hcalE = new TH1D("h_hcalE"," ; HCal E",200,0,4);
       HList.Add(h_hcalE);











  

  //peds_vs_xy->SetMarkerStyle(8);

  peds_diff_vs_xy->SetMarkerStyle(8);
  int_diff_vs_xy->SetMarkerStyle(8);
  //ch_diff_vs_xy->Write("ch_diff_vs_row_col");
  peds_width_diff_vs_xy->SetMarkerStyle(8);

  peds_diff_vs_ch->SetMarkerStyle(8);
  int_diff_vs_ch->SetMarkerStyle(8);
  //ch_diff_vs_ch->Write("ch_diff_vs_ch");
  peds_width_diff_vs_ch->SetMarkerStyle(8);
  peds_cnf_diff_vs_ch->SetMarkerStyle(8);
  time_diff_vs_ch->SetMarkerStyle(8);

  time_diff_from_expect_vs_ch->SetMarkerStyle(8);

  peds_vs_xy->SetMarkerStyle(8);
  int_vs_xy->SetMarkerStyle(8);
  time_vs_xy->SetMarkerStyle(8);
  //ch_diff_vs_xy->Write("ch_diff_vs_row_col");
  peds_width_vs_xy->SetMarkerStyle(8);

  peds_vs_ch->SetMarkerStyle(8);
  int_vs_ch->SetMarkerStyle(8);
  //ch_diff_vs_ch->Write("ch_diff_vs_ch");
  peds_width_vs_ch->SetMarkerStyle(8);
  time_vs_ch->SetMarkerStyle(8);



  TChain chain1("T");
  //chain1.Add(Form("/adaqfs/home/a-onl/sbs/Rootfiles/ecal_%i_-1.root", first_file));

  chain1.Add(Form("/adaqfs/home/a-onl/sbs/Rootfiles/gep5_replayed_%i_50k_events.root", first_file));
  
  chain1.SetBranchStatus("*",0);
  chain1.SetBranchStatus("*ecal*",1);
  //chain1.SetBranchStatus("*MC*",1);



  //primary clusters are clusters with *ideally* highest energy
  
  chain1.SetBranchAddress("Ndata.earm.ecal.a_p", &ndata_adc1);//ndata array defining number of adc events in an entry
  //chain1.SetBranchAddress("earm.ecal.a_amp_p", adc_amp_mv1);
  chain1.SetBranchAddress("earm.ecal.adcelemID", adc_ch_id1);
  chain1.SetBranchAddress("earm.ecal.a_p", adc_int_pC1);
  chain1.SetBranchAddress("earm.ecal.adcrow", adc_row1);
  chain1.SetBranchAddress("earm.ecal.adccol", adc_col1);
  chain1.SetBranchAddress("earm.ecal.ped", adc_ped1);
  chain1.SetBranchAddress("earm.ecal.adcxpos", adc_xpos1);
  chain1.SetBranchAddress("earm.ecal.adcypos", adc_ypos1);
  chain1.SetBranchAddress("earm.ecal.a_time", adc_time1);
  
  /*
    chain1.SetBranchAddress("Ndata.earm.ecal.clus_blk.id", &ndata_clus_blk1);//ndata array defining number of events for blocks found in primary clusters
    chain1.SetBranchAddress("earm.ecal.clus_blk.id", clus_blk_ch_id1);
    chain1.SetBranchAddress("earm.ecal.clus_blk.e", clus_blk_e1);
    chain1.SetBranchAddress("earm.ecal.clus_blk.row", clus_blk_row1);
    chain1.SetBranchAddress("earm.ecal.clus_blk.col", clus_blk_col1);
    chain1.SetBranchAddress("earm.ecal.clus_blk.x", clus_blk_x1);
    chain1.SetBranchAddress("earm.ecal.clus_blk.y", clus_blk_y1);

    chain1.SetBranchAddress("Ndata.earm.ecal.clus.id", &ndata_clus1);//ndata array defining number of events for all clusters
    //chain1.SetBranchAddress("Ndata.earm.ecal.clus.e", ndata_clus_arr1);
    chain1.SetBranchAddress("earm.ecal.clus.id", clus_ch_id1);
    chain1.SetBranchAddress("earm.ecal.clus.e", clus_e1);
    chain1.SetBranchAddress("earm.ecal.clus.row", clus_row1);
    chain1.SetBranchAddress("earm.ecal.clus.col", clus_col1);
    chain1.SetBranchAddress("earm.ecal.clus.x", clus_x1);
    chain1.SetBranchAddress("earm.ecal.clus.y", clus_y1);
    */
  Int_t nentries = chain1.GetEntries();
  
  for(Int_t ent=0; ent<nentries; ent++){
    chain1.GetEntry(ent);
    for(Int_t j=0; j<ndata_adc1; j++){//get info from the blk's in primary clusters for this entry
      
      adc_peds1_per_block[int (adc_ch_id1[j])]->Fill(adc_ped1[j]);
      adc_int1_per_block[int (adc_ch_id1[j])]->Fill(adc_int_pC1[j]);
      adc_time1_per_block[int (adc_ch_id1[j])]->Fill(adc_time1[j]);
      
      adc_peds1->Fill(adc_ped1[j]);
      adc_energy1->Fill(adc_int_pC1[j]);
      //adc_hits1->Fill(adc_1[j]);
      
    }
  }
  
  //chain1.Reset();
  

  
  
  
  
  
  TChain chain2("T");
  
  chain2.Add(Form("/adaqfs/home/a-onl/sbs/Rootfiles/ecal_%i_50000.root", second_file));
  //chain2.Add(Form("/adaqfs/home/a-onl/sbs/Rootfiles/gep5_replayed_%i_50k_events.root", second_file));
  //chain2.Add(second_file);
  
  chain2.SetBranchStatus("*",0);
  chain2.SetBranchStatus("*ecal*",1);
  //chain1.SetBranchStatus("*MC*",1);






  

  Int_t vtp_nclus;
  chain2->SetBranchAddress("Ndata.earm.ecal.vtp.clus.e",&vtp_nclus) ;
  chain2->SetBranchStatus("Ndata.earm.ecal.vtp.clus.e",1) ; 
  Double_t vtp_clus_detid;
  chain2->SetBranchAddress("earm.ecal.vtp.detid",&vtp_clus_detid) ;
  chain2->SetBranchStatus("earm.ecal.vtp.detid",1) ; 
  Double_t vtp_clus_e[1000];
  chain2->SetBranchAddress("earm.ecal.vtp.clus.e",&vtp_clus_e) ;
  chain2->SetBranchStatus("earm.ecal.vtp.clus.e",1) ; 
  Double_t vtp_clus_x[1000];
  chain2->SetBranchAddress("earm.ecal.vtp.clus.x",&vtp_clus_x) ;
  chain2->SetBranchStatus("earm.ecal.vtp.clus.x",1) ; 
  Double_t vtp_clus_y[1000];
  chain2->SetBranchAddress("earm.ecal.vtp.clus.y",&vtp_clus_y) ;
  chain2->SetBranchStatus("earm.ecal.vtp.clus.y",1) ; 
  Double_t vtp_clus_size[1000];
  chain2->SetBranchAddress("earm.ecal.vtp.clus.size",&vtp_clus_size) ;
  chain2->SetBranchStatus("earm.ecal.vtp.clus.size",1) ; 
  Double_t vtp_clus_time[1000];
  chain2->SetBranchAddress("earm.ecal.vtp.clus.time",&vtp_clus_time) ;
  chain2->SetBranchStatus("earm.ecal.vtp.clus.time",1) ;





  
  





  
  
  //primary clusters are clusters with *ideally* highest energy
  
  chain2.SetBranchAddress("Ndata.earm.ecal.a_p", &ndata_adc2);//ndata array defining number of adc events in an entry
  //chain2.SetBranchAddress("earm.ecal.a_amp_p", adc_amp_mv2);
  chain2.SetBranchAddress("earm.ecal.adcelemID", adc_ch_id2);
  chain2.SetBranchAddress("earm.ecal.a_p", adc_int_pC2);
  chain2.SetBranchAddress("earm.ecal.adcrow", adc_row2);
  chain2.SetBranchAddress("earm.ecal.adccol", adc_col2);
  chain2.SetBranchAddress("earm.ecal.ped", adc_ped2);
  chain2.SetBranchAddress("earm.ecal.adcxpos", adc_xpos2);
  chain2.SetBranchAddress("earm.ecal.adcypos", adc_ypos2);
  chain2.SetBranchAddress("earm.ecal.a_time", adc_time2);

  /*
    chain1.SetBranchAddress("Ndata.earm.ecal.clus_blk.id", &ndata_clus_blk1);//ndata array defining number of events for blocks found in primary clusters
    chain1.SetBranchAddress("earm.ecal.clus_blk.id", clus_blk_ch_id1);
    chain1.SetBranchAddress("earm.ecal.clus_blk.e", clus_blk_e1);
    chain1.SetBranchAddress("earm.ecal.clus_blk.row", clus_blk_row1);
    chain1.SetBranchAddress("earm.ecal.clus_blk.col", clus_blk_col1);
    chain1.SetBranchAddress("earm.ecal.clus_blk.x", clus_blk_x1);
    chain1.SetBranchAddress("earm.ecal.clus_blk.y", clus_blk_y1);
  */
    chain2.SetBranchAddress("Ndata.earm.ecal.clus.id", &ndata_clus1);//ndata array defining number of events for all clusters
    //chain2.SetBranchAddress("Ndata.earm.ecal.clus.e", ndata_clus_arr1);
    chain2.SetBranchAddress("earm.ecal.clus.id", clus_ch_id1);
    chain2.SetBranchAddress("earm.ecal.clus.e", clus_e1);
    chain2.SetBranchAddress("earm.ecal.clus.row", clus_row1);
    chain2.SetBranchAddress("earm.ecal.clus.col", clus_col1);
    chain2.SetBranchAddress("earm.ecal.clus.x", clus_x1);
    chain2.SetBranchAddress("earm.ecal.clus.y", clus_y1);
    
  Int_t nentries2 = chain2.GetEntries();
  
  for(Int_t ent=0; ent<nentries2; ent++){
    chain2.GetEntry(ent);
    for(Int_t j=0; j<ndata_adc2; j++){//get info from the blk's in primary clusters for this entry
      adc_peds2_per_block[int (adc_ch_id2[j])]->Fill(adc_ped2[j]);
      adc_int2_per_block[int (adc_ch_id2[j])]->Fill(adc_int_pC2[j]);
      adc_time2_per_block[int (adc_ch_id2[j])]->Fill(adc_time2[j]);
      
      adc_peds2->Fill(adc_ped2[j]);
      
    }
  }

  //chain2.Reset();
  
  double mean_ped1[1656];
  double mean_ped2[1656];
  
  double mean_time1[1656];
  double mean_time2[1656];

  double mean_int1[1656];
  double mean_int2[1656];

  bool ch_check1[1656];
  bool ch_check2[1656];

  double ped_width1[1656];
  double ped_width2[1656];
  
  for(int i = 0; i < 1656; i++){

    if(adc_int1_per_block[i]->GetEntries() <= 5){
      ch_check1[i] = 0;
    }else{
      ch_check1[i] = 1;
    }

    if(adc_int2_per_block[i]->GetEntries() <= 5){
      ch_check2[i] = 0;
    }else{
      ch_check2[i] = 1;
    }

    if(ch_check1[i] != ch_check2[i]){
      list_of_changed_ch[i] = 1;
    }else{
      list_of_changed_ch[i] = 0;
    }

    ch_diff_vs_xy->Fill( col_arr[i]+1, row_arr[i]+1, list_of_changed_ch[i] );
    ch_diff_vs_ch->Fill( i, list_of_changed_ch[i] );
    
    mean_ped1[i] = adc_peds1_per_block[i]->GetMean();
    mean_ped2[i] = adc_peds2_per_block[i]->GetMean();

    mean_time1[i] = adc_time1_per_block[i]->GetMean();
    mean_time2[i] = adc_time2_per_block[i]->GetMean();

    //cout << mean_ped1[i] << " vs " << mean_ped2[i] << endl;

    //fitting is being difficult, so maybe just try finding the first and last bins filled to define the peak
    //ped_width1[i] = abs(adc_peds1_per_block[i]->FindLastBinAbove() - adc_peds1_per_block[i]->FindFirstBinAbove());
    //ped_width2[i] = abs(adc_peds2_per_block[i]->FindLastBinAbove() - adc_peds2_per_block[i]->FindFirstBinAbove());

    //want to try GetRMS instead of fitting

    ped_width1[i] = adc_peds1_per_block[i]->GetRMS();
    ped_width2[i] = adc_peds2_per_block[i]->GetRMS();
    
    //peds_diff_vs_xy->Fill( xpos_arr[i]/100.0, ypos_arr[i]/100.0, abs(mean_ped2[i]-mean_ped1[i]) );
    /*
    TF1 *gaus_fit1 = new TF1("gaus_fit1","gaus");
    int maxBin1 = adc_peds1_per_block[i]->GetMaximumBin();
    double maxVal1 = adc_peds1_per_block[i]->GetMaximum();
    double binWidth1 = adc_peds1_per_block[i]->GetBinWidth(maxBin1);
    double stdDev1 = adc_peds1_per_block[i]->GetStdDev();
    double maxBinCenter1 = adc_peds1_per_block[i]->GetXaxis()->GetBinCenter(maxBin1);
    double lowBin1 = maxBin1*binWidth1 - 2*stdDev1;
    double upBin1 = maxBin1*binWidth1 + 2*stdDev1;

    gaus_fit1->SetParameters( maxVal1, maxBinCenter1, stdDev1);
    //gaus_fit1->SetRange(1.0, 4.0);

    adc_peds1_per_block[i]->Fit(gaus_fit1, "Q");
    
    TF1 *gaus_fit2 = new TF1("gaus_fit2","gaus");
    int maxBin2 = adc_peds2_per_block[i]->GetMaximumBin();
    double maxVal2 = adc_peds2_per_block[i]->GetMaximum();
    double binWidth2 = adc_peds2_per_block[i]->GetBinWidth(maxBin2);
    double stdDev2 = adc_peds2_per_block[i]->GetStdDev();
    double maxBinCenter2 = adc_peds2_per_block[i]->GetXaxis()->GetBinCenter(maxBin2);
    double lowBin2 = maxBin2*binWidth2 - 2*stdDev2;
    double upBin2 = maxBin2*binWidth2 + 2*stdDev2;

    gaus_fit2->SetParameters( maxVal2, maxBinCenter2, stdDev2);
    //gaus_fit2->SetRange(140, 450);

    adc_peds2_per_block[i]->Fit(gaus_fit2, "Q");
    */
    //ped_width2[i] = gaus_fit2->GetParameter(2);
    //cout << ped_width2[i] << endl;
    //ped_width1[i] = gaus_fit1->GetParameter(2);

    peds_diff_vs_xy->Fill( col_arr[i]+1, row_arr[i]+1, abs(mean_ped2[i]-mean_ped1[i]) );
    peds_width_diff_vs_xy->Fill( col_arr[i]+1, row_arr[i]+1, abs(ped_width2[i]-ped_width1[i]) );
    peds_diff_vs_ch->Fill( i, abs(mean_ped2[i]-mean_ped1[i]) );
    peds_width_diff_vs_ch->Fill( i, abs(ped_width2[i]-ped_width1[i]) );
    peds_cnf_diff_vs_ch->Fill( i, abs(mean_ped2[i]-ecal_ped_cnf[i]) );

    peds_vs_xy->Fill( col_arr[i]+1, row_arr[i]+1, mean_ped2[i] );
    peds_width_vs_xy->Fill( col_arr[i]+1, row_arr[i]+1, ped_width2[i] );
    peds_vs_ch->Fill( i, mean_ped2[i] );
    peds_width_vs_ch->Fill( i, ped_width2[i] );
    
    
    mean_int1[i] = adc_int1_per_block[i]->GetMean();
    mean_int2[i] = adc_int2_per_block[i]->GetMean();

    //mean_int1[i] = gaus_fit1->GetParameter(1);
    //mean_int2[i] = gaus_fit1->GetParameter(1);
    //int_diff_vs_xy->Fill( xpos_arr[i]/100.0, ypos_arr[i]/100.0, abs(mean_int2[i]-mean_int1[i]) );
    int_diff_vs_xy->Fill( col_arr[i]+1, row_arr[i]+1, abs(mean_int2[i]-mean_int1[i]) );
    int_diff_vs_ch->Fill( i, abs(mean_int2[i]-mean_int1[i]) );
    
    int_vs_xy->Fill( col_arr[i]+1, row_arr[i]+1, abs(mean_int2[i]) );
    int_vs_ch->Fill( i, abs(mean_int2[i]) );

    time_diff_vs_xy->Fill( col_arr[i]+1, row_arr[i]+1, (mean_time2[i]-mean_time1[i]) );
    time_diff_vs_ch->Fill( i, (mean_time2[i]-mean_time1[i]) );
    time_diff_from_expect_vs_ch->Fill( i, (mean_time2[i]-195) );
    
    time_vs_xy->Fill( col_arr[i]+1, row_arr[i]+1, abs(mean_time2[i]) );
    time_vs_ch->Fill( i, abs(mean_time2[i]) );
    
  }

  //TFile f(Form("/adaqfs/home/a-onl/sbs/Rootfiles/ecal_%i_-1.root", second_file), "UPDATE"); //use this when we want to write these monitor histos to the actual rootfile for use with panguin
  
  TFile f(Form("hist/test_monitor_output_%i.root",second_file), "RECREATE");
  //TFile f(Form("/adaqfs/home/a-onl/sbs/Rootfiles/ecal_%i_-1.root", second_file), "UPDATE");
  //TFile f(replayed_root_file, "UPDATE")
  peds_diff_vs_xy->SetMarkerStyle(8);
  int_diff_vs_xy->SetMarkerStyle(8);
  //ch_diff_vs_xy->Write("ch_diff_vs_row_col");
  peds_width_diff_vs_xy->SetMarkerStyle(8);

  peds_diff_vs_ch->SetMarkerStyle(8);
  int_diff_vs_ch->SetMarkerStyle(8);
  //ch_diff_vs_ch->Write("ch_diff_vs_ch");
  peds_width_diff_vs_ch->SetMarkerStyle(8);
  peds_cnf_diff_vs_ch->SetMarkerStyle(8);
  time_diff_vs_ch->SetMarkerStyle(8);

  time_diff_from_expect_vs_ch->SetMarkerStyle(8);

  peds_vs_xy->SetMarkerStyle(8);
  int_vs_xy->SetMarkerStyle(8);
  time_vs_xy->SetMarkerStyle(8);
  //ch_diff_vs_xy->Write("ch_diff_vs_row_col");
  peds_width_vs_xy->SetMarkerStyle(8);

  peds_vs_ch->SetMarkerStyle(8);
  int_vs_ch->SetMarkerStyle(8);
  //ch_diff_vs_ch->Write("ch_diff_vs_ch");
  peds_width_vs_ch->SetMarkerStyle(8);
  time_vs_ch->SetMarkerStyle(8);

  TH1D *y_proj_time = time_vs_ch->ProjectionY("y_proj_time", 0, 1656);
  y_proj_time->Write("y_proj_time");



  
  peds_diff_vs_xy->Write("peds_diff_vs_row_col");
  int_diff_vs_xy->Write("int_diff_vs_row_col");
  //ch_diff_vs_xy->Write("ch_diff_vs_row_col");
  peds_width_diff_vs_xy->Write("peds_width_diff_vs_xy");

  peds_diff_vs_ch->Write("peds_diff_vs_ch");
  int_diff_vs_ch->Write("int_diff_vs_ch");
  //ch_diff_vs_ch->Write("ch_diff_vs_ch");
  peds_width_diff_vs_ch->Write("peds_width_diff_vs_ch");
  peds_cnf_diff_vs_ch->Write("peds_cnf_diff_vs_ch");
  time_diff_vs_ch->Write("time_diff_vs_ch");

  time_diff_from_expect_vs_ch->Write("time_diff_from_expect_vs_ch");

  TH2D *test_vs_xy = new TH2D("test_vs_row_col", "peds_vs_row_col", 27, 0.5, 27.5, 69, 0.5, 69.5);

  for (int i = 1; i <= test_vs_xy->GetNbinsX(); ++i) {
        for (int j = 1; j <= test_vs_xy->GetNbinsY(); ++j) {
	  
        }
    }
  
  peds_vs_xy->Write("peds_vs_row_col");
  int_vs_xy->Write("int_vs_row_col");
  time_vs_xy->Write("time_vs_row_col");
  //ch_diff_vs_xy->Write("ch_diff_vs_row_col");
  peds_width_vs_xy->Write("peds_width_vs_row_col");
  
  peds_vs_ch->Write("peds_vs_ch");
  int_vs_ch->Write("int_vs_ch");
  //ch_diff_vs_ch->Write("ch_diff_vs_ch");
  peds_width_vs_ch->Write("peds_width_vs_ch");
  time_vs_ch->Write("time_vs_ch");

  //TH1D *y_proj_time = time_vs_ch->ProjectionY("y_proj_time", 0, 1656);
  //y_proj_time->Write("y_proj_time");
  
  /*
  for(int i = 0; i < 5; i++){
    adc_peds2_per_block[i]->Write(Form("peds_example%i",i));
  }
  */



  gApplication->Terminate(0);
}



