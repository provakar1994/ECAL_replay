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

void ecal_monitor(int first_file, int second_file){//compare second file to first file

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

      if( ntokens == 6 ){
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

  TH2D *peds_diff_vs_xy = new TH2D("peds_diff_vs_xy", "", 27, 0.5, 27.5, 69, 0.5, 69.5);
  TH2D *int_diff_vs_xy = new TH2D("int_diff_vs_xy", "", 27, 0.5, 27.5, 69, 0.5, 69.5);
  TH2D *ch_diff_vs_xy = new TH2D("ch_diff_vs_xy", "", 27, 0.5, 27.5, 69, 0.5, 69.5);
  TH2D *peds_width_diff_vs_xy = new TH2D("peds_width_diff_vs_xy", "", 27, 0.5, 27.5, 69, 0.5, 69.5);

  TH2D *peds_diff_vs_ch = new TH2D("peds_diff_vs_ch", "", 1656, 0.5, 1656.5, 10, 0, 5);
  TH2D *peds_width_diff_vs_ch = new TH2D("peds_width_diff_vs_ch", "", 1656, 0.5, 1656.5, 15, 0, 50);
  TH2D *int_diff_vs_ch = new TH2D("int_diff_vs_ch", "", 1656, 0.5, 1656.5, 800, 0, 20);
  TH2D *ch_diff_vs_ch = new TH2D("ch_diff_vs_ch", "", 1656, 0.5, 1656.5, 10, 0, 10);
  TH2D *peds_cnf_diff_vs_ch = new TH2D("peds_diff_vs_ch", "", 1656, 0.5, 1656.5, 80, 0, 10);
  
  for(int i = 0; i < 1656; i++){
    adc_peds1_per_block[i] = new TH1F("", "", 300, 50, 150);
    adc_peds2_per_block[i] = new TH1F("", "", 300, 50, 150);

    adc_int1_per_block[i] = new TH1F("", "", 100, 0, 80);
    adc_int2_per_block[i] = new TH1F("", "", 100, 0, 80);
  }
  
  TH1F *adc_peds1 = new TH1F("adc_peds1", "", 300, 50, 150);
  TH1F *adc_peds2 = new TH1F("adc_peds2", "", 300, 50, 150);

  TH1F *adc_energy1 = new TH1F("adc_energy1", "", 100, 130, 160);
  TH1F *adc_energy2 = new TH1F("adc_energy2", "", 100, 130, 160);

  TH1F *adc_hits1 = new TH1F("adc_hits1", "", 100, 130, 160);
  TH1F *adc_hits2 = new TH1F("adc_hits2", "", 100, 130, 160);

  TChain chain1("T");
  chain1.Add(Form("/adaqfs/home/a-onl/sbs/Rootfiles/ecal_%i_-1.root", first_file));
  
  chain1.SetBranchStatus("*",0);
  chain1.SetBranchStatus("*ecal*",1);
  //chain1.SetBranchStatus("*MC*",1);
  
  //primary clusters are clusters with *ideally* highest energy
  
  chain1.SetBranchAddress("Ndata.earm.ecal.a_amp_p", &ndata_adc1);//ndata array defining number of adc events in an entry
  chain1.SetBranchAddress("earm.ecal.a_amp_p", adc_amp_mv1);
  chain1.SetBranchAddress("earm.ecal.adcelemID", adc_ch_id1);
  chain1.SetBranchAddress("earm.ecal.a_p", adc_int_pC1);
  chain1.SetBranchAddress("earm.ecal.adcrow", adc_row1);
  chain1.SetBranchAddress("earm.ecal.adccol", adc_col1);
  chain1.SetBranchAddress("earm.ecal.ped", adc_ped1);
  chain1.SetBranchAddress("earm.ecal.adcxpos", adc_xpos1);
  chain1.SetBranchAddress("earm.ecal.adcypos", adc_ypos1);
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
      
      adc_peds1->Fill(adc_ped1[j]);
      adc_energy1->Fill(adc_int_pC1[j]);
      //adc_hits1->Fill(adc_1[j]);
      
    }
  }
  
  
  

  
  
  
  
  
  TChain chain2("T");
  
  chain2.Add(Form("/adaqfs/home/a-onl/sbs/Rootfiles/ecal_%i_-1.root", second_file));
  
  chain2.SetBranchStatus("*",0);
  chain2.SetBranchStatus("*ecal*",1);
  //chain1.SetBranchStatus("*MC*",1);
  
  //primary clusters are clusters with *ideally* highest energy
  
  chain2.SetBranchAddress("Ndata.earm.ecal.a_amp_p", &ndata_adc2);//ndata array defining number of adc events in an entry
  chain2.SetBranchAddress("earm.ecal.a_amp_p", adc_amp_mv2);
  chain2.SetBranchAddress("earm.ecal.adcelemID", adc_ch_id2);
  chain2.SetBranchAddress("earm.ecal.a_p", adc_int_pC2);
  chain2.SetBranchAddress("earm.ecal.adcrow", adc_row2);
  chain2.SetBranchAddress("earm.ecal.adccol", adc_col2);
  chain2.SetBranchAddress("earm.ecal.ped", adc_ped2);
  chain2.SetBranchAddress("earm.ecal.adcxpos", adc_xpos2);
  chain2.SetBranchAddress("earm.ecal.adcypos", adc_ypos2);
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
  Int_t nentries2 = chain2.GetEntries();
  
  for(Int_t ent=0; ent<nentries2; ent++){
    chain2.GetEntry(ent);
    for(Int_t j=0; j<ndata_adc2; j++){//get info from the blk's in primary clusters for this entry
      adc_peds2_per_block[int (adc_ch_id2[j])]->Fill(adc_ped2[j]);
      adc_int2_per_block[int (adc_ch_id2[j])]->Fill(adc_int_pC2[j]);
      
      adc_peds2->Fill(adc_ped2[j]);
      
    }
  }
  
  double mean_ped1[1656];
  double mean_ped2[1656];

  double mean_int1[1656];
  double mean_int2[1656];

  bool ch_check1[1656];
  bool ch_check2[1656];

  int ped_width1[1656];
  int ped_width2[1656];
  
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

    //cout << mean_ped1[i] << " vs " << mean_ped2[i] << endl;

    //fitting is being difficult, so maybe just try finding the first and last bins filled to define the peak
    ped_width1[i] = abs(adc_peds1_per_block[i]->FindLastBinAbove() - adc_peds1_per_block[i]->FindFirstBinAbove());
    ped_width2[i] = abs(adc_peds2_per_block[i]->FindLastBinAbove() - adc_peds2_per_block[i]->FindFirstBinAbove());
    
    //peds_diff_vs_xy->Fill( xpos_arr[i]/100.0, ypos_arr[i]/100.0, abs(mean_ped2[i]-mean_ped1[i]) );
    peds_diff_vs_xy->Fill( col_arr[i]+1, row_arr[i]+1, abs(mean_ped2[i]-mean_ped1[i]) );
    peds_width_diff_vs_xy->Fill( col_arr[i]+1, row_arr[i]+1, abs(ped_width2[i]-ped_width1[i]) );
    peds_diff_vs_ch->Fill( i, abs(mean_ped2[i]-mean_ped1[i]) );
    peds_width_diff_vs_ch->Fill( i, abs(ped_width2[i]-ped_width1[i]) );
    peds_cnf_diff_vs_ch->Fill( i, abs(mean_ped2[i]-ecal_ped_cnf[i]) );
    /*
    TF1 *gaus_fit1 = new TF1("gaus_fit1","gaus");
    int maxBin1 = adc_int1_per_block[i]->GetMaximumBin();
    double maxVal1 = adc_int1_per_block[i]->GetMaximum();
    double binWidth1 = adc_int1_per_block[i]->GetBinWidth(maxBin1);
    double stdDev1 = adc_int1_per_block[i]->GetStdDev();
    double maxBinCenter1 = adc_int1_per_block[i]->GetXaxis()->GetBinCenter(maxBin1);
    double lowBin1 = maxBin1*binWidth1 - 2*stdDev1;
    double upBin1 = maxBin1*binWidth1 + 2*stdDev1;

    gaus_fit1->SetParameters( maxVal1, maxBinCenter1, stdDev1);
    gaus_fit1->SetRange(1.0, 4.0);

    adc_int1_per_block[i]->Fit(gaus_fit1,"R");

    TF1 *gaus_fit2 = new TF1("gaus_fit2","gaus");
    int maxBin2 = adc_int2_per_block[i]->GetMaximumBin();
    double maxVal2 = adc_int2_per_block[i]->GetMaximum();
    double binWidth2 = adc_int2_per_block[i]->GetBinWidth(maxBin2);
    double stdDev2 = adc_int2_per_block[i]->GetStdDev();
    double maxBinCenter2 = adc_int2_per_block[i]->GetXaxis()->GetBinCenter(maxBin2);
    double lowBin2 = maxBin2*binWidth2 - 2*stdDev2;
    double upBin2 = maxBin2*binWidth2 + 2*stdDev2;

    gaus_fit2->SetParameters( maxVal2, maxBinCenter2, stdDev2);
    gaus_fit2->SetRange(140, 450);

    adc_int2_per_block[i]->Fit(gaus_fit2,"R");
    */
    mean_int1[i] = adc_int1_per_block[i]->GetMean();
    mean_int2[i] = adc_int2_per_block[i]->GetMean();

    //mean_int1[i] = gaus_fit1->GetParameter(1);
    //mean_int2[i] = gaus_fit1->GetParameter(1);
    //int_diff_vs_xy->Fill( xpos_arr[i]/100.0, ypos_arr[i]/100.0, abs(mean_int2[i]-mean_int1[i]) );
    int_diff_vs_xy->Fill( col_arr[i]+1, row_arr[i]+1, abs(mean_int2[i]-mean_int1[i]) );
    int_diff_vs_ch->Fill( i, abs(mean_int2[i]-mean_int1[i]) );
  }

  //TFile f(Form("/adaqfs/home/a-onl/sbs/Rootfiles/ecal_%i_-1.root", second_file), "UPDATE"); //use this when we want to write these monitor histos to the actual rootfile for use with panguin
  
  TFile f(Form("test_monitor_output.root"), "RECREATE");
  peds_diff_vs_xy->Write("peds_diff_vs_row_col");
  int_diff_vs_xy->Write("int_diff_vs_row_col");
  ch_diff_vs_xy->Write("ch_diff_vs_row_col");
  peds_width_diff_vs_xy->Write("peds_width_diff_vs_xy");

  peds_diff_vs_ch->Write("peds_diff_vs_ch");
  int_diff_vs_ch->Write("int_diff_vs_ch");
  ch_diff_vs_ch->Write("ch_diff_vs_ch");
  peds_width_diff_vs_ch->Write("peds_width_diff_vs_ch");
  peds_cnf_diff_vs_ch->Write("peds_cnf_diff_vs_ch");
}



