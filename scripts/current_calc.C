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

void current_calc(int run){

  TChain *C = new TChain("TSsbs");
  //TFile *C = new TFile ("TSsbs");
  //C->Add(rootfilename);
  C->Add(Form("/adaqfs/home/a-onl/sbs/Rootfiles/gep5_replayed_%i_full.root",run));

  //chain.SetBranchStatus("*",0);
  C->SetBranchStatus("sbs*", 1);

  //the 1 3 and 10 refers to gain

  //beam current is given in microamps, so multiplying the given current by time in seconds gives microcoulombs

  Double_t bcm_u1_curr;
  C->SetBranchAddress("sbs.bcm.u1.current", &bcm_u1_curr);
  Double_t bcm_u3_curr;
  C->SetBranchAddress("sbs.bcm.u3.current", &bcm_u3_curr);
  Double_t bcm_u10_curr;
  C->SetBranchAddress("sbs.bcm.u10.current", &bcm_u10_curr);
  Double_t bcm_d1_curr;
  C->SetBranchAddress("sbs.bcm.d1.current", &bcm_d1_curr);
  Double_t bcm_d3_curr;
  C->SetBranchAddress("sbs.bcm.d3.current", &bcm_d3_curr);
  Double_t bcm_d10_curr;
  C->SetBranchAddress("sbs.bcm.d10.current", &bcm_d10_curr);

  Double_t d104k_clk_count;
  C->SetBranchAddress("sbs.104kHz_CLK.cnt", &d104k_clk_count);
  Double_t d104k_clk_rate;
  C->SetBranchAddress("sbs.104kHz_CLK.rate", &d104k_clk_rate);


  Int_t nentries = C->GetEntries();
  //cout << nentries << endl;

  Double_t prev_time = 0;
  Double_t charge_accumulated_u1 = 0;
  Double_t charge_accumulated_u3 = 0;
  Double_t charge_accumulated_u10 = 0;
  Double_t charge_accumulated_d1 = 0;
  Double_t charge_accumulated_d3 = 0;
  Double_t charge_accumulated_d10 = 0;

  Double_t tot_time_sec = 0;
  
  for(Int_t i=0; i<nentries; i++){

    C->GetEntry(i);
    
    Double_t time_sec = d104k_clk_count/d104k_clk_rate;
    Double_t time_min = time_sec/60;

    Double_t time_interval = time_sec - prev_time;

    if(time_interval >=0.01 && time_interval <= 5){
      tot_time_sec += time_interval;
    
      Double_t charge_accu_this_int_u1 = bcm_u1_curr*time_interval;
      charge_accumulated_u1 += charge_accu_this_int_u1;
      Double_t charge_accu_this_int_u3 = bcm_u3_curr*time_interval;
      charge_accumulated_u3 += charge_accu_this_int_u3;
      Double_t charge_accu_this_int_u10 = bcm_u10_curr*time_interval;
      charge_accumulated_u10 += charge_accu_this_int_u10;
      Double_t charge_accu_this_int_d1 = bcm_d1_curr*time_interval;
      charge_accumulated_d1 += charge_accu_this_int_d1;
      Double_t charge_accu_this_int_d3 = bcm_d3_curr*time_interval;
      charge_accumulated_d3 += charge_accu_this_int_d3;
      Double_t charge_accu_this_int_d10 = bcm_d10_curr*time_interval;
      charge_accumulated_d10 += charge_accu_this_int_d10;
      
      cout << time_sec << " " << bcm_u1_curr << endl;
    }
    prev_time = time_sec;
    
    
  }
  printf("Time of run  %7.2f min \n",tot_time_sec/60);
  printf("Charge Monitors  (Micro Coulombs)\n");
  printf("Upstream BCM   gain x1 %8.2f     x3 %8.2f     x10 %8.2f\n",
	 charge_accumulated_u1, charge_accumulated_u3, charge_accumulated_u10);
  printf("Downstream BCM   gain x1 %8.2f     x3 %8.2f     x10 %8.2f\n",
	 charge_accumulated_d1, charge_accumulated_d3, charge_accumulated_d10);
}
