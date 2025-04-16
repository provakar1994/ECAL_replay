#include <TSystem.h>
#include <TChain.h>
#include <TString.h>
#include "TFile.h"
#include "TTree.h"
#include <TNtuple.h>
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TH1F.h"
#include <TH2.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPaveLabel.h>
#include <TProfile.h>
#include <TPolyLine.h>
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;

void make_hist(const char* replayed_file_name){
  gStyle->SetPalette(1,0);
  gStyle->SetOptStat(1000011);
  gStyle->SetOptFit(11);
  gStyle->SetTitleOffset(1.,"Y");
  gStyle->SetTitleOffset(.7,"X");
  gStyle->SetLabelSize(0.04,"XY");
  gStyle->SetTitleSize(0.06,"XY");
  gStyle->SetPadLeftMargin(0.12);
   //

  
  TString basename="ecal";
  TString Chainroot;
  //Chainroot = Form("~/sbs/Rootfiles/%s_%d_%d.root",basename.Data(),nrun,nev);
  TChain *fchain = new TChain("T");
  //fchain->Add(Chainroot);
  //cout << "check file " << Chainroot << endl;
  //Int_t npart = 1;
  //Chainroot = Form("~/sbs/Rootfiles/%s_%d_%d_%d.root",basename.Data(),nrun,nev,npart);
  //
  //while (gSystem->FindFile(".",Chainroot)) { 
  //cout << " add file " << Chainroot << endl;
  fchain->Add(Chainroot);
  //npart++;
  //Chainroot = Form("~/sbs/Rootfiles/%s_%d_%d_%d.root",basename.Data(),nrun,nev,npart);
  //}
  
 fchain->Add(replayed_file_name);
 
  //
 //TString fullname = Form("%s_%d_%d",basename.Data(),nrun,nev);
 //TString outputhist;
 //outputhist= "hist/"+fullname+"_hist.root";
  TObjArray HList(0);
  static const Int_t MaxHit=10000;
  fchain->SetBranchStatus("*",0);
  Int_t vtp_nclus;
   fchain->SetBranchAddress("Ndata.earm.ecal.vtp.clus.e",&vtp_nclus) ;
   fchain->SetBranchStatus("Ndata.earm.ecal.vtp.clus.e",1) ; 
  Double_t vtp_clus_detid;
   fchain->SetBranchAddress("earm.ecal.vtp.detid",&vtp_clus_detid) ;
   fchain->SetBranchStatus("earm.ecal.vtp.detid",1) ; 
  Double_t vtp_clus_e[1000];
   fchain->SetBranchAddress("earm.ecal.vtp.clus.e",&vtp_clus_e) ;
   fchain->SetBranchStatus("earm.ecal.vtp.clus.e",1) ; 
  Double_t vtp_clus_x[1000];
   fchain->SetBranchAddress("earm.ecal.vtp.clus.x",&vtp_clus_x) ;
   fchain->SetBranchStatus("earm.ecal.vtp.clus.x",1) ; 
  Double_t vtp_clus_y[1000];
   fchain->SetBranchAddress("earm.ecal.vtp.clus.y",&vtp_clus_y) ;
   fchain->SetBranchStatus("earm.ecal.vtp.clus.y",1) ; 
  Double_t vtp_clus_size[1000];
   fchain->SetBranchAddress("earm.ecal.vtp.clus.size",&vtp_clus_size) ;
   fchain->SetBranchStatus("earm.ecal.vtp.clus.size",1) ; 
  Double_t vtp_clus_time[1000];
   fchain->SetBranchAddress("earm.ecal.vtp.clus.time",&vtp_clus_time) ;
   fchain->SetBranchStatus("earm.ecal.vtp.clus.time",1) ;
   //
 Int_t hvtp_nclus;
   fchain->SetBranchAddress("Ndata.sbs.hcal.vtp.clus.e",&hvtp_nclus) ;
   fchain->SetBranchStatus("Ndata.sbs.hcal.vtp.clus.e",1) ; 
  Double_t hvtp_clus_detid;
   fchain->SetBranchAddress("sbs.hcal.vtp.detid",&hvtp_clus_detid) ;
   fchain->SetBranchStatus("sbs.hcal.vtp.detid",1) ; 
  Double_t hvtp_clus_e[1000];
   fchain->SetBranchAddress("sbs.hcal.vtp.clus.e",&hvtp_clus_e) ;
   fchain->SetBranchStatus("sbs.hcal.vtp.clus.e",1) ; 
  Double_t hvtp_clus_x[1000];
   fchain->SetBranchAddress("sbs.hcal.vtp.clus.x",&hvtp_clus_x) ;
   fchain->SetBranchStatus("sbs.hcal.vtp.clus.x",1) ; 
  Double_t hvtp_clus_y[1000];
   fchain->SetBranchAddress("sbs.hcal.vtp.clus.y",&hvtp_clus_y) ;
   fchain->SetBranchStatus("sbs.hcal.vtp.clus.y",1) ; 
  Double_t hvtp_clus_size[1000];
   fchain->SetBranchAddress("sbs.hcal.vtp.clus.size",&hvtp_clus_size) ;
   fchain->SetBranchStatus("sbs.hcal.vtp.clus.size",1) ; 
  Double_t hvtp_clus_time[1000];
   fchain->SetBranchAddress("sbs.hcal.vtp.clus.time",&hvtp_clus_time) ;
   fchain->SetBranchStatus("sbs.hcal.vtp.clus.time",1) ;
   //
 Int_t ecal_nclus_hits;
   fchain->SetBranchAddress("Ndata.earm.ecal.clus.e",&ecal_nclus_hits) ;
   fchain->SetBranchStatus("Ndata.earm.ecal.clus.e",1) ;
 Double_t ecal_clus_e[100];
   fchain->SetBranchAddress("earm.ecal.clus.e",&ecal_clus_e) ;
   fchain->SetBranchStatus("earm.ecal.clus.e",1) ;
 Double_t ecal_clus_atime[100];
   fchain->SetBranchAddress("earm.ecal.clus.atimeblk",&ecal_clus_atime) ;
   fchain->SetBranchStatus("earm.ecal.clus.atimeblk",1) ;
 //
      Double_t ecal_nclus; // raw adc
     fchain->SetBranchAddress("earm.ecal.nclus",&ecal_nclus) ;
     fchain->SetBranchStatus("earm.ecal.nclus",1) ;
      Double_t hcal_nclus; // raw adc
     fchain->SetBranchAddress("sbs.hcal.nclus",&hcal_nclus) ;
     fchain->SetBranchStatus("sbs.hcal.nclus",1) ;
      Double_t hcal_eblk; // raw adc
     fchain->SetBranchAddress("sbs.hcal.eblk",&hcal_eblk) ;
     fchain->SetBranchStatus("sbs.hcal.eblk",1) ;
      Double_t hcal_e; // raw adc
     fchain->SetBranchAddress("sbs.hcal.e",&hcal_e) ;
     fchain->SetBranchStatus("sbs.hcal.e",1) ;
    Double_t ecal_e; // raw adc
     fchain->SetBranchAddress("earm.ecal.e",&ecal_e) ;
     fchain->SetBranchStatus("earm.ecal.e",1) ;
    Double_t ecal_atimeblk; // raw adc
     fchain->SetBranchAddress("earm.ecal.atimeblk",&ecal_atimeblk) ;
     fchain->SetBranchStatus("earm.ecal.atimeblk",1) ;
   Double_t hcal_atimeblk; // raw adc
     fchain->SetBranchAddress("sbs.hcal.atimeblk",&hcal_atimeblk) ;
     fchain->SetBranchStatus("sbs.hcal.atimeblk",1) ;
    Double_t ecal_eblk; // raw adc
     fchain->SetBranchAddress("earm.ecal.eblk",&ecal_eblk) ;
     fchain->SetBranchStatus("earm.ecal.eblk",1) ;
    Double_t ecal_colblk; // raw adc
     fchain->SetBranchAddress("earm.ecal.colblk",&ecal_colblk) ;
     fchain->SetBranchStatus("earm.ecal.colblk",1) ;
    Double_t ecal_rowblk; // raw adc
     fchain->SetBranchAddress("earm.ecal.rowblk",&ecal_rowblk) ;
     fchain->SetBranchStatus("earm.ecal.rowblk",1) ;
    Double_t ecal_idblk; // raw adc
     fchain->SetBranchAddress("earm.ecal.idblk",&ecal_idblk) ;
     fchain->SetBranchStatus("earm.ecal.idblk",1) ;
    Double_t ecal_nblk; // raw adc
     fchain->SetBranchAddress("earm.ecal.nblk",&ecal_nblk) ;
     fchain->SetBranchStatus("earm.ecal.nblk",1) ;
     // ecal hit array
   Int_t ecalHits;
   fchain->SetBranchAddress("Ndata.earm.ecal.a_amp_p",&ecalHits) ;
   fchain->SetBranchStatus("Ndata.earm.ecal.a_amp_p",1) ;
    Double_t ecal_ped[MaxHit]; // raw adc
     fchain->SetBranchAddress("earm.ecal.ped",&ecal_ped) ;
     fchain->SetBranchStatus("earm.ecal.ped",1) ;
   Double_t ecal_time[MaxHit]; // raw adc
     fchain->SetBranchAddress("earm.ecal.a_time",&ecal_time) ;
     fchain->SetBranchStatus("earm.ecal.a_time",1) ;
   Double_t ecal_amp[MaxHit]; // raw adc
     fchain->SetBranchAddress("earm.ecal.a_amp_p",&ecal_amp) ;
     fchain->SetBranchStatus("earm.ecal.a_amp_p",1) ;
    Double_t ecal_int[MaxHit]; // raw adc
     fchain->SetBranchAddress("earm.ecal.a_p",&ecal_int) ;
     fchain->SetBranchStatus("earm.ecal.a_p",1) ;
   Double_t ecal_id[MaxHit]; // raw adc
     fchain->SetBranchAddress("earm.ecal.adcelemID",&ecal_id) ;
     fchain->SetBranchStatus("earm.ecal.adcelemID",1) ;
   Double_t ecal_row[MaxHit]; // raw adc
     fchain->SetBranchAddress("earm.ecal.adcrow",&ecal_row) ;
     fchain->SetBranchStatus("earm.ecal.adcrow",1) ;
    Double_t ecal_col[MaxHit]; // raw adc
     fchain->SetBranchAddress("earm.ecal.adccol",&ecal_col) ;
     fchain->SetBranchStatus("earm.ecal.adccol",1) ;
   Double_t ecal_x[MaxHit]; // raw adc
     fchain->SetBranchAddress("earm.ecal.adcxpos",&ecal_x) ;
     fchain->SetBranchStatus("earm.ecal.adcxpos",1) ;
   Double_t ecal_y[MaxHit]; // raw adc
     fchain->SetBranchAddress("earm.ecal.adcypos",&ecal_y) ;
     fchain->SetBranchStatus("earm.ecal.adcypos",1) ;
//
       Double_t sbs_hcal_nclus = 0.;
  Double_t hcal_reftdc[7];
  Double_t hcal_reftdc_elemID[7];

  Int_t Ndata_hcal_reftdc=0;
  
  TH1D *h1_ECAL = new TH1D("h1_ECAL","ECAL Trig time; ns", 500, 0.1, 500.1);
       HList.Add(h1_ECAL);
  TH1D *h1_HCAL = new TH1D("h1_HCAL","HCAL Trig time; ns", 500,0.1, 500.1);
        HList.Add(h1_HCAL);
  TH2D *h1_HCAL_ECAL = new TH2D("h1_HCAL_ECAL",";HCAL Trig time ns;ECAL Trig time ns", 500,0.1, 500.1, 500,0.1, 500.1);
        HList.Add(h1_HCAL_ECAL);
 TH1D *h1_COIN = new TH1D("h1_COIN","COIN Trig time; ns", 500, 0.1, 500.1);
        HList.Add(h1_COIN);
  TH1D *h1_calDiff = new TH1D("h1_calDiff","HCal/ECAL Difference; ns", 500, -250, 250);
         HList.Add(h1_calDiff);
 TH1D *h1_RF_time = new TH1D("h1_RF_time","Accelerator RF Time; ns;",351,-0.5,350.5);
         HList.Add(h1_RF_time);

  double trigbits;
  
  // Declare branches
  fchain->SetBranchStatus("sbs.hcal.Ref.tdc",1);
  fchain->SetBranchStatus("sbs.hcal.Ref.tdcelemID",1);
  fchain->SetBranchStatus("Ndata.sbs.hcal.Ref.tdcelemID",1);
  fchain->SetBranchStatus("g.trigbits",1);
  
  fchain->SetBranchAddress("sbs.hcal.nclus", &sbs_hcal_nclus);
  fchain->SetBranchAddress("sbs.hcal.Ref.tdc",hcal_reftdc);
  fchain->SetBranchAddress("sbs.hcal.Ref.tdcelemID",hcal_reftdc_elemID);
  fchain->SetBranchAddress("Ndata.sbs.hcal.Ref.tdcelemID",&Ndata_hcal_reftdc);
  fchain->SetBranchAddress("g.trigbits", &trigbits);

     //
     TH2D* h_vnclus_enclus;
     h_vnclus_enclus = new TH2D("h_vnclus_enclus"," VTP DETID =12; VTP nclus ; ECAL nclus",40,0,40,20,0,20);
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
      TH2D* h_hvtpEclus_hcalEclus;
     h_hvtpEclus_hcalEclus = new TH2D("h_hvtpEclus_hcalEclus"," ; hVTP E clus ; HCAL E clus",200,0,1000,200,0,1);
       HList.Add(h_hvtpEclus_hcalEclus);
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
     
     //
 int nbadRF_HCAL = 0;
    vector<Int_t> goodhits;
  Int_t nentries = fchain->GetEntries();
   cout << " Total Entry = " << nentries << endl;
  for (int ii = 0; ii < nentries; ii++) {
      fchain->GetEntry(ii);
     if (ii%50000==0) cout << " Entry = " << ii << endl;
     //
       Double_t coin_time=0., hcal_time=0, ecal_time=0, rf_time=0;

    bool anyHCALRF = false;
    
    for( int ihit=0; ihit<Ndata_hcal_reftdc; ihit++ ){
      if( hcal_reftdc_elemID[ihit]==3 ){
	rf_time = hcal_reftdc[ihit];
	anyHCALRF = true;
      }
      if( hcal_reftdc_elemID[ihit]==4 ){ //BigBite trigger time:
	coin_time = hcal_reftdc[ihit];
      }
      if( hcal_reftdc_elemID[ihit]==5 ){ //BigBite trigger time:
	ecal_time = hcal_reftdc[ihit];
      }
      if( hcal_reftdc_elemID[ihit]==6 ){ //BigBite trigger time:
	hcal_time = hcal_reftdc[ihit];
      }
    }
    
    if( !anyHCALRF ) nbadRF_HCAL++;
    
    //Take difference to eliminate 25ns jitter

    //if( trigbits == 4 ){
      h1_ECAL->Fill(ecal_time);
      h1_HCAL->Fill(hcal_time);
      h1_HCAL_ECAL->Fill(hcal_time,ecal_time);
      h1_calDiff->Fill(ecal_time-hcal_time);
      h1_COIN->Fill(coin_time);
      
      h1_RF_time->Fill(rf_time);
  
     //
    if (vtp_clus_detid>=0) {
      h_vnclus_enclus->Fill(float(vtp_nclus),float(ecal_nclus_hits));
	h_hvnclus_hnclus->Fill(float(hvtp_nclus),float(ecal_nclus_hits));
      Int_t max_vtp_index=0;
      for (Int_t nc=0;nc<vtp_nclus;nc++) {
	if (vtp_clus_e[nc]>=vtp_clus_e[max_vtp_index]) max_vtp_index=nc;
	//	cout << nc << " " << vtp_clus_e[nc] << " " << max_vtp_index <<endl;
        }
      Int_t max_hvtp_index=0;
      for (Int_t nc=0;nc<vtp_nclus;nc++) {
	if (hvtp_clus_e[nc]>=hvtp_clus_e[max_hvtp_index]) max_hvtp_index=nc;
	//	cout << nc << " " << hvtp_clus_e[nc] << " " << max_hvtp_index <<endl;
        }
      if (hvtp_nclus>0 && max_hvtp_index>=0 ) {
        h_hvtpEclus_hcalEclus->Fill(hvtp_clus_e[max_hvtp_index],hcal_e);
	h_hvtpTclus->Fill(hvtp_clus_time[max_hvtp_index]);
	h_hvtpEclus->Fill(hvtp_clus_e[max_hvtp_index]);
	h_hvtpSizeclus->Fill(hvtp_clus_size[max_hvtp_index]);
	h_hvtp_x_y->Fill(hvtp_clus_x[max_hvtp_index],hvtp_clus_y[max_hvtp_index]);
	h_hvtp_x_y_time->Fill(hvtp_clus_x[max_hvtp_index],hvtp_clus_y[max_hvtp_index],hvtp_clus_time[max_hvtp_index]);
	  }
      if (vtp_nclus>0 && max_vtp_index>=0 && ecal_nclus_hits>0) {
	h_vtpEclus_ecalEclus->Fill(vtp_clus_e[max_vtp_index],ecal_e);
	h_vtpEclus->Fill(vtp_clus_e[max_vtp_index]);
	h_ecalSize_ecalEclus->Fill(ecal_nblk,ecal_e);
	h_vtpSizeclus->Fill(vtp_clus_size[max_vtp_index]);
	h_hcalTime_ecaltime->Fill(hcal_atimeblk,ecal_atimeblk);
	h_hcalTime->Fill(hcal_atimeblk);
	h_hcalE->Fill(hcal_e);
	h_vtpTclus->Fill(vtp_clus_time[max_vtp_index]);
  	h_ecalEclus->Fill(ecal_e);
 	h_ecalSizeclus->Fill(ecal_nblk);
 	h_ecalTclus->Fill(ecal_atimeblk);
        h_ecal_row_col->Fill(ecal_colblk+1,ecal_rowblk+1); 
        h_vtp_row_col->Fill(vtp_clus_x[max_vtp_index],vtp_clus_y[max_vtp_index]); 
        h_ecal_id_e->Fill(ecal_idblk,ecal_eblk); 
        h_ecal_row_col_e->Fill(ecal_colblk+1,ecal_rowblk+1,ecal_eblk); 
        h_ecal_row_col_t->Fill(ecal_colblk+1,ecal_rowblk+1,ecal_atimeblk); 
      }	
    } else {
     h_vnclus_enclus_detidzero->Fill(float(vtp_nclus),float(ecal_nclus_hits));
       Int_t max_vtp_index=0;
      for (Int_t nc=0;nc<vtp_nclus;nc++) {
	if (vtp_clus_e[nc]>=vtp_clus_e[max_vtp_index]) max_vtp_index=nc;
	//	cout << nc << " " << vtp_clus_e[nc] << " " << max_vtp_index <<endl;
        }
      if (vtp_nclus>0 && max_vtp_index>=0 && ecal_nclus_hits>0) {
	h_vtpEclus_ecalEclus_detidzero->Fill(vtp_clus_e[max_vtp_index],ecal_e);
      }
   }
  }
    //
  h_ecal_row_col_e->Divide(h_ecal_row_col);
  h_ecal_row_col_t->Divide(h_ecal_row_col);
  h_hvtp_x_y_time->Divide(h_hvtp_x_y);
  //
  //TFile hsimc(outputhist,"recreate");
  TFile hsimc(Form("hist/clus_vtp.root"),"UPDATE");
   HList.Write();
  //
}
