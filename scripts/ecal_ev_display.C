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

void ecal_ev_display(const char* inelastic_sim_name){

  bool fit_arg = 0;

  //cout << "Type 0 for Gaus fit, type 1 for expo fit." << endl;
  //cin >> fit_arg;
  
  int arr_index = 0;
  int cells = 1656;

  double ex[cells], mean_ECAL[cells], sigma_ECAL[cells], dmean_ECAL[cells], dsigma_ECAL[cells];

  double cells_ECAL[cells];
  
  
  int row_arr[cells];
  int col_arr[cells];
  double xpos_arr[cells];
  double ypos_arr[cells];

  double mc_ep[2000];
  int ecal_adcID[1656];
  
  ifstream ecalfile("maps/ECAL_r_c_x_y_cpr.csv");

  TString currentline;

  int check = 0;

  while( currentline.ReadLine(ecalfile) ){

    if( check != 0 ){ //new cell information:
      TObjArray *tokens = ( (TObjArray*) currentline.Tokenize(",") );
      int ntokens = tokens->GetEntries();

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

      }
    }
    check = 1;
  }

  TH1F *htemp;

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
  
  TChain chain("T");
  chain.Add(inelastic_sim_name);

  chain.SetBranchStatus("*",0);
  chain.SetBranchStatus("*ecal*",1);
  chain.SetBranchStatus("*MC*",1);

  //primary clusters are clusters with *ideally* highest energy

  chain.SetBranchAddress("Ndata.earm.ecal.a_amp_p", &ndata_adc);//ndata array defining number of adc events in an entry
  chain.SetBranchAddress("earm.ecal.a_amp_p", adc_amp_mv);
  chain.SetBranchAddress("earm.ecal.adcelemID", adc_ch_id);
  chain.SetBranchAddress("earm.ecal.a_p", adc_int_pC);
  chain.SetBranchAddress("earm.ecal.adcrow", adc_row);
  chain.SetBranchAddress("earm.ecal.adccol", adc_col);

  chain.SetBranchAddress("Ndata.earm.ecal.clus_blk.id", &ndata_clus_blk);//ndata array defining number of events for blocks found in primary clusters
  chain.SetBranchAddress("earm.ecal.clus_blk.id", clus_blk_ch_id);
  chain.SetBranchAddress("earm.ecal.clus_blk.e", clus_blk_e);
  chain.SetBranchAddress("earm.ecal.clus_blk.row", clus_blk_row);
  chain.SetBranchAddress("earm.ecal.clus_blk.col", clus_blk_col);
  chain.SetBranchAddress("earm.ecal.clus_blk.x", clus_blk_x);
  chain.SetBranchAddress("earm.ecal.clus_blk.y", clus_blk_y);

  chain.SetBranchAddress("Ndata.earm.ecal.clus.id", &ndata_clus);//ndata array defining number of events for all blocks before any clustering
  chain.SetBranchAddress("Ndata.earm.ecal.clus.e", ndata_clus_arr);
  chain.SetBranchAddress("earm.ecal.clus.id", clus_ch_id);
  chain.SetBranchAddress("earm.ecal.clus.e", clus_e);
  chain.SetBranchAddress("earm.ecal.clus.row", clus_row);
  chain.SetBranchAddress("earm.ecal.clus.col", clus_col);
  chain.SetBranchAddress("earm.ecal.clus.x", clus_x);
  chain.SetBranchAddress("earm.ecal.clus.y", clus_y);

  bool plot_type = 0;

  cout << "Type 0 to plot energy colz for primary clusters, Type 1 to plot groups of clusters before primary cut" << endl;
  cin >> plot_type;

  //cs->Clear();
  //cs->Divide(4,4);
  /*
  TH1F *test = new TH1F("test", "", 100, 1.5, 3.0);
  cs->cd(1);
  
  int new_print = 0;
  */

  Int_t nentries = chain.GetEntries();

  //for(Int_t ent=0; ent<nentries; ent++){
  for(Int_t ent=1; ent<10; ent++){

    chain.GetEntry(ent);

    /*
    for(Int_t j=0; j<ndata_earm_ecal_a; j++){//get the info from the adc arrays for this entry

    }
    */
    TH2D *e_vs_rowcol = new TH2D("e_vs_rowcol", "", 27, 1, 27, 69, 1, 69);
    TH2D *e_vs_xy = new TH2D("e_vs_xy", "", 27, -0.61, 0.65, 69, -1.5, 1.5);
    e_vs_xy->SetTitle(Form("energy of blocks in primary cluster, ev: %i", ent));

    TH2D *clus_vs_xy = new TH2D("clus_vs_xy", "", 27, -0.61, 0.65, 69, -1.5, 1.5);
    clus_vs_xy->SetTitle(Form("group blocks by cluster before cutting out non-primary clusters, ev: %i", ent));

    if(plot_type == 0 && ndata_clus_blk >= 1){
    
      for(Int_t j=0; j<ndata_clus_blk; j++){//get info from the blk's in primary clusters for this entry
	//want 3d histo row vs col, with energy deposition as z color axis, 2d histo with energy weight   
	e_vs_xy->Fill(clus_blk_y[j], clus_blk_x[j], clus_blk_e[j]);
      }
      auto cs = new TCanvas("cs", "",  700, 900);
      //cs->Divide(2);
      cs->Clear();
      cs->cd(0);
      cs->SetCanvasSize(1000, 1000);
      cs->SetWindowSize(1000, 1000);
      e_vs_xy->Draw("Colz");
      
      cs->SetLogz();
      gPad->Modified();//needed to make drawn canvas show up in every iteration of loop
      gPad->Update();

      gPad->WaitPrimitive();
      
      //cout << " Continue? " << endl;
      //cin.ignore();	     			     
			       
    }

    if(plot_type == 1 && ndata_clus >= 1){
    //if(plot_type == 1){
      cout << "HERE " << endl << endl << endl << endl;
    
      for(Int_t j=0; j<ndata_clus; j++){//get info from the blk's in primary clusters for this entry
	//want 3d histo row vs col, with energy deposition as z color axis, 2d histo with energy weight   
	//clus_vs_xy->Fill(clus_y[j], clus_x[j], ndata_clus_arr[j]);
	clus_vs_xy->Fill(clus_y[j], clus_x[j], ndata_clus_arr[j]);
	cout << ndata_clus_arr[j] << endl;
      }
      auto cs = new TCanvas("cs", "",  700, 900);
      //cs->Divide(2);
      cs->Clear();
      cs->cd(0);
      cs->SetCanvasSize(1000, 1000);
      cs->SetWindowSize(1000, 1000);
      clus_vs_xy->Draw("Colz");
      
      //cs->SetLogz();
      gPad->Modified();//needed to make drawn canvas show up in every iteration of loop
      gPad->Update();
      
      gPad->WaitPrimitive();
      
      //cout << " Continue? " << endl;
      //cin.ignore();	     			     
			       
    }

  }

}



















