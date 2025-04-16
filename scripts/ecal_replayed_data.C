#include <iostream>
#include <fstream>
#include <string>
#include <TFile.h>
#include <TDirectory.h>
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

void ecal_replayed_data(int runis){

  int arr_index = 0;
  int cells = 1656;


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

      if( ntokens == 7 ){
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
  
  Int_t ndata_adc;
  Double_t adc_amp_mv[1656];
  Double_t adc_ch_id[1656];
  Double_t adc_int_pC[1656];
  Double_t adc_row[1656];
  Double_t adc_col[1656];

  Int_t ndata_clus_blk;
  Double_t clus_blk_ch_id[1656];
  Double_t clus_blk_e[1656];
  Double_t clus_blk_row[1656];
  Double_t clus_blk_col[1656];
  Double_t clus_blk_x[1656];
  Double_t clus_blk_y[1656];

  Int_t ndata_clus;
  Double_t ndata_clus_arr[1656];
  Double_t clus_ch_id[1656];
  Double_t clus_e[1656];
  Double_t clus_row[1656];
  Double_t clus_col[1656];
  Double_t clus_x[1656];
  Double_t clus_y[1656];

  Int_t ndata_vtp;
  //Double_t ndata_clus_arr[1656];
  //Double_t clus_ch_id[1656];
  Double_t vtp_e[1656];
  // Double_t clus_row[1656];
  //Double_t clus_col[1656];
  Double_t vtp_x[1656];
  Double_t vtp_y[1656];
  Double_t vtp_detid[1656];









  double Parsamp[3]={0.};
  double ParErrsamp[3]={0.};
  double Parstime[3]={0.};
  double ParErrstime[3]={0.};
  double Parsint[3]={0.};
  double ParErrsint[3]={0.};
  double Parsped[3]={0.};
  double ParErrsped[3]={0.};
  double hamp_min = 0.; 
  double hamp_max = 45.;






  


  TH2F* time_peak_vs_ch = new TH2F("time_peak_vs_ch", "", 1656, 0, 1656, 200, 0, 350);

  TH1F* adc_amp_histo[cells];
  TH1F* adc_time_histo[cells];

  //Int_t ndata_adc;
  Double_t adc_amp[cells];
  Double_t adc_time[cells];

  for(int i = 0; i < cells; i++){
    adc_amp_histo[i] = new TH1F("","", 90, 0, 210);
    adc_time_histo[i] = new TH1F("","", 400, 0, 400);
  }

  TChain chain("T");
  //chain.Add(Form("/adaqfs/home/a-onl/sbs/Rootfiles/ecal_%i_-1.root",runis));
  chain.Add(Form("/adaqfs/home/a-onl/sbs/Rootfiles/gep5_replayed_%i_full.root",runis));
  

  
  chain.SetBranchStatus("*",0);
  chain.SetBranchStatus("*ecal*",1);
  //chain.SetBranchStatus("*MC*",1);

  //primary clusters are clusters with *ideally* highest energy

  chain.SetBranchAddress("Ndata.earm.ecal.adcelemID", &ndata_adc);//ndata array defining number of adc events in an entry
  chain.SetBranchAddress("earm.ecal.a_amp_p", adc_amp_mv);
  chain.SetBranchAddress("earm.ecal.adcelemID", adc_ch_id);
  chain.SetBranchAddress("earm.ecal.a_p", adc_int_pC);
  chain.SetBranchAddress("earm.ecal.adcrow", adc_row);
  chain.SetBranchAddress("earm.ecal.adccol", adc_col);
  chain.SetBranchAddress("earm.ecal.a_time", adc_time);

  chain.SetBranchAddress("Ndata.earm.ecal.clus_blk.id", &ndata_clus_blk);//ndata array defining number of events for blocks found in primary clusters
  chain.SetBranchAddress("earm.ecal.clus_blk.id", clus_blk_ch_id);
  chain.SetBranchAddress("earm.ecal.clus_blk.e", clus_blk_e);
  chain.SetBranchAddress("earm.ecal.clus_blk.row", clus_blk_row);
  chain.SetBranchAddress("earm.ecal.clus_blk.col", clus_blk_col);
  chain.SetBranchAddress("earm.ecal.clus_blk.x", clus_blk_x);
  chain.SetBranchAddress("earm.ecal.clus_blk.y", clus_blk_y);

  chain.SetBranchAddress("Ndata.earm.ecal.clus.id", &ndata_clus);//ndata array defining number of events for all blocks before any clustering
  //chain.SetBranchAddress("Ndata.earm.ecal.clus.e", ndata_clus_arr);
  chain.SetBranchAddress("earm.ecal.clus.id", clus_ch_id);
  chain.SetBranchAddress("earm.ecal.clus.e", clus_e);
  chain.SetBranchAddress("earm.ecal.clus.row", clus_row);
  chain.SetBranchAddress("earm.ecal.clus.col", clus_col);
  chain.SetBranchAddress("earm.ecal.clus.x", clus_x);
  chain.SetBranchAddress("earm.ecal.clus.y", clus_y);

  chain.SetBranchAddress("Ndata.earm.ecal.vtp.clus.size", &ndata_vtp);//ndata array defining number of events for blocks found in primary clusters
  //chain.SetBranchAddress("earm.ecal.clus_blk.id", clus_blk_ch_id);
  chain.SetBranchAddress("earm.ecal.vtp.clus.e", vtp_e);
  //chain.SetBranchAddress("earm.ecal.clus_blk.row", clus_blk_row);
  //chain.SetBranchAddress("earm.ecal.clus_blk.col", clus_blk_col);
  chain.SetBranchAddress("earm.ecal.vtp.clus.x", vtp_x);
  chain.SetBranchAddress("earm.ecal.vtp.clus.y", vtp_y);
  chain.SetBranchAddress("earm.ecal.vtp.detid", vtp_detid);


  int time_mean[1656] = {0};
  int time_bad_list[1656] = {0};
  int time_entries_list[1656] = {0};
  int time_outlier_list[1656] = {0};

  Int_t nentries = chain.GetEntries();

  TH1D *raw_occu = new TH1D("raw_occu", "", 100, 0.5, 100.5);

  

  for(Int_t i=0; i<nentries; i++){

    chain.GetEntry(i);


    int k = 0;
    for(Int_t j=0; j<ndata_adc; j++){

      //cout << adc_time[j] << " " << adc_amp[j] << endl; 

      Int_t pmtindx1 = adc_ch_id[j];

      adc_amp_histo[pmtindx1]->Fill(adc_amp[j]);
      adc_time_histo[pmtindx1]->Fill(adc_time[j]);

      if(k+1 == ndata_adc){
	//raw occupancy: how many hits per event
	raw_occu->Fill(double (k));
      }
      k++;

    }
  }










  //fitting

  TF1 *fgaus2 = new TF1("fgaus","gaus");

  for(Int_t n = 0; n < 1656; n++){
    int lowerBinCamp;
    int upperBinCamp; 
    int lowerBinCtime; 
    int upperBinCtime; 
    int lowerBinCint; 
    int upperBinCint; 
    int lowerBinCped; 
    int upperBinCped;

    int maxBintime = adc_time_histo[n]->GetMaximumBin();
    double maxBinCentertime = adc_time_histo[n]->GetXaxis()->GetBinCenter( maxBintime );
    double maxCounttime = adc_time_histo[n]->GetMaximum();
    double binWidthtime = adc_time_histo[n]->GetBinWidth(maxBintime);
    double stdDevtime = adc_time_histo[n]->GetStdDev();

         
    if(adc_time_histo[n]->GetEntries()>=0){ 
      
      // Create fit functions for each module
      fgaus2->SetLineColor(2);
      fgaus2->SetNpx(1000);
      
      // first fit
      
      if(stdDevtime > 15){
	lowerBinCtime = (maxBintime)*binWidthtime - 15.0;
	upperBinCtime = (maxBintime)*binWidthtime + 15.0;
      }else{
	lowerBinCtime = (maxBintime)*binWidthtime - 1.0*stdDevtime;
	upperBinCtime = (maxBintime)*binWidthtime + 1.0*stdDevtime;
	
      }
      
      fgaus2->SetParameters( maxCounttime,maxBinCentertime,stdDevtime );
      fgaus2->SetRange( lowerBinCtime, upperBinCtime );
      //h_pmt_is[n]->Fit(fgaus2,"NO+RQ");
      fgaus2->GetParameters(Parstime);
      
      
      adc_time_histo[n]->Fit( fgaus2,"+RQ" );//kip
      time_mean[n] = fgaus2->GetParameter(1);
      //cout << n << " " << fgaus2->GetParameter(0) << " " << fgaus2->GetParameter(1) << endl;
      fgaus2->GetParameters(Parstime);
      
      time_outlier_list[n] = 215 - time_mean[n];
      //cout time_mean[n];
      
      //cout << n << " " << time_mean[n] << endl;
      
      //if(adc_time_histo[n]->GetEntries() < 100 && abs(time_outlier_list[n] > 20)){
      
      
      if(abs(time_outlier_list[n]) > 11){
	time_bad_list[n] = 1;
      }
      
      time_entries_list[n] = adc_time_histo[n]->GetEntries();
      
      
      
    }
    
    
  }

  TString time_outlier_txt;
  
  time_outlier_txt = Form("textfile/time_outliers_%i.txt", runis);
  
  ofstream time_bad_file(time_outlier_txt);
  
  time_bad_file << "Ch, nEntries, offset from avg, time peak(ns)" << endl;
  
  for(int i = 0; i < 1656; i++){
    if(time_bad_list[i] == 1){
      time_bad_file << i << " " << time_entries_list[i] << " " << time_outlier_list[i] << " " << time_mean[i] << endl;
      cout  << "Channel index: " << i << ", num of entries for that channel: " << time_entries_list[i] << ", difference in ns from desired time:  " << time_outlier_list[i] << ", timing peak: " << time_mean[i] << endl;
    }
  }
  cout << "this list of channels with bad timing peaks is also a text file found at /textfiles/time_outliers_" << runis << ".txt, now you should use the root script generate_fadc_timing_delays.C in order to get the new delays in order to fix these bad timing peaks." << endl;
  
  





  ofstream time_mean_list;
  time_mean_list.open(Form("textfile/time_peak_list_run%i.txt",runis));

  
  
  for(int n=0; n < 1656; n++){
    time_mean_list << n << " " << time_mean[n] << endl;
    time_peak_vs_ch->Fill(n, time_mean[n]);
  }

  time_mean_list.close();

  TFile f1(Form("hist/time_outlier_fits%i.root",runis), "RECREATE");
  for(int i = 0; i < 1656; i++){
    if(time_bad_list[i] == 1){
      adc_time_histo[i]->Write(Form("time_outlier_ch%i",i));
    }
  }

  //TFile f(Form("timing_distribution_%i.root",runis), "RECREATE");
  time_peak_vs_ch->Write("time_peak_vs_ch");
  raw_occu->Write("raw_occu");
  
  gSystem->Exit(0);


  

}
