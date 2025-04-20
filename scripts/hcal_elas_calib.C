/*
  This script has been prepared for the energy calibration of ECAL detector located in GEpV electron arm. 
  It does so by minimizing the chi2 of the difference between calorimeter cluster energy and the reconstructed 
  scattered electron momentum calculated based on the reconstructed proton track assuming elastic scattering. 
  It gets the old adc gain coefficients (GeV/pC) from tree and writes out the new adc gain coeffs. and 
  ratios (New/Old) to file. One needs a configfile to execute this script [see cfg/example.cfg]. To execute, do:
  ----
  [a-onl@aonl2 macros]$ pwd
  /adaqfs/home/a-onl/sbs/ECal_replay/scripts
  [a-onl@aonl2 scripts]$ root -l 
  root [0] .x elas_calib.C("cfg/example.cfg")
  ----
  P. Datta  <pdbforce@jlab.org>  CREATED  07 April 2025
*/

#include <iostream>
#include <fstream>
#include <vector>

#include "LookUpTableReader.h"

/*
  Convention:
  Pei + Ppi = Pef + Ppf
  where,
  Pei = (E_beam,0,0,E_beam), 
  Ppi = (Mp,0,0,0),
  Pef = (p_e, p_e*sin(th_e), 0, p_e*cos(th_e)),
  Ppf = (E_p, p_p*sin(th_p), 0, p_p*cos(th_p))
*/

// 
Double_t const Mp = 0.938272081;  // +/- 6E-9 GeV

// ECAL geometry
const int kNblks = 288;
const int kNrows = 24;
const int kNcols = 12; 

// Utility functions
string GetDate();
void CustmProfHisto(TH1*);
TString GetOutFileBase(TString);
std::vector<int> ReadNumFromTextFileToV(const std::string &filename);
void BuildMatrix(const double *A, TVectorD &B, TMatrixD &M, Double_t p_e);
std::vector<std::string> SplitString(char const delim, std::string const str);
std::unordered_set<int> GetEdgeBlkID(LookUpTableReader& ecalmap, bool isdedug);
std::unordered_map<int,int> GetNColPerRow(LookUpTableReader& ecalmap, bool isdebug);
std::vector<std::pair<double, double>> GetXYPosBADChan(LookUpTableReader& ecalmap, std::vector<int>& badch, bool isdebug);
bool ISNearBADChan(std::vector<std::pair<double, double>> const &xyBADchan, double xSEED, double ySEED, double r2_max, int debug);
  
// Main function
void hcal_elas_calib(char const *configfilename,
		std::string prefix="",  //prefix to the output file names
		int isdebug=0)          //0=>False, >0=>True
{
  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  // Define a clock to check macro processing time
  TStopwatch *sw = new TStopwatch();
  TStopwatch *sw2 = new TStopwatch();
  sw->Start(); sw2->Start();  

  TChain *C = new TChain("T");
  //creating base for outfile names
  TString cfgfilebase = GetOutFileBase(configfilename);

  // Defining important variables
  // misc
  TString macros_dir, badchan_file;
  Int_t Nmin=10;
  Double_t minMBratio=0.1;  
  // HCAL related
  Double_t ECAL_dist=8.0; //m
  Double_t ECAL_voff=0., ECAL_hoff=0., ECAL_zoff=0.; //m
  // histogram related
  Double_t h_W2_bin=200, h_W2_min=0., h_W2_max=5.;
  Double_t h_Q2_bin=200, h_Q2_min=0., h_Q2_max=5.;  
  Double_t h_eHCAL_bin=200, h_eHCAL_min=0., h_eHCAL_max=5.;
  Double_t h_EovEexp_bin=200, h_EovEexp_min=0., h_EovEexp_max=5., EovEexp_fit_width=1.5;
  Double_t h2_EovEexp_bin=200, h2_EovEexp_min=0., h2_EovEexp_max=5.;
  Double_t h_eSF_bin=200, h_eSF_min=0., h_eSF_max=0.45, eSF_fit_width=1.;  
  Double_t h2_p_bin = 200, h2_p_min = 0., h2_p_max = 5.;  
  Double_t h2_p_coarse_bin=25, h2_p_coarse_min=0., h2_p_coarse_max=5.;
  // analysis related
  bool cut_on_W2=0, cut_on_EovEexp=0, cut_on_eHCAL=0, cut_on_nearBADchan=0;
  bool indtADC=0, isHCALedge=0, isNearBadChan=0;
  Double_t r_max, r2_max;
  Double_t W2_mean, W2_sigma, W2_nsigma;
  Double_t dt_mean, dt_sigma, dt_nsigma;
  Double_t dt_mean_1, dt_sigma_1, dt_nsigma_1;
  Double_t dt_mean_2, dt_sigma_2, dt_nsigma_2;      
  Double_t EovEexp_cut_limit, eHCAL_cut_limit;
  Double_t hit_threshold=0., engFrac_cut=0., tmax_cut=1000.;    
  // kinematic related
  Double_t E_beam, th_bb;
  Double_t p_p, p_px, p_py, p_pz, th_p, ph_p, E_p, KE_p, KE_p_frac, p_p_offset; 
  Double_t th_e, ph_e;
  Double_t nu, Q2, W2, dx, dy;
  Double_t nu_4vec, Q2_4vec, W2_4vec, dx_4vec, dy_4vec;  
  // -- Defining variables for scattered e- momentum
  Double_t p_e;      // calculated using recontructed momnetum
  Double_t p_e_th;   // calculated using recontructed angles
  Double_t p_e_CURR; // this variable will be used for calibration
  // analysis related
  int recon_type=1;  // 1=> use p_e. Otherwise, p_e_th is used  
  // bool enable_4vec=false;
  Double_t HCAL_sf=0.07; // HCAL sampling fraction
  
  TMatrixD M(kNblks,kNblks), M_inv(kNblks,kNblks);
  TVectorD B(kNblks), CoeffR(kNblks);  
  Double_t A[kNblks];
  
  // Reading config file
  ifstream configfile(configfilename);
  char runlistfile[1000]; 
  TString currentline, readline;
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endRunlist") ){
    if( !currentline.BeginsWith("#") ){
      sprintf(runlistfile,"%s",currentline.Data());
      ifstream run_list(runlistfile);
      while( readline.ReadLine( run_list ) && !readline.BeginsWith("endlist") ){
  	if( !readline.BeginsWith("#") ){
	  std::cout << readline << "\n";
	  C->Add(readline);
  	}
      }   
    } 
  }  
  TCut globalcut = ""; TString gcutstr;
  while (currentline.ReadLine(configfile) && !currentline.BeginsWith("endcut")) {
    if (!currentline.BeginsWith("#")) {
      globalcut += currentline;
      gcutstr += currentline;
    }    
  }
  std::vector<std::string> gCutList = SplitString('&', gcutstr.Data());
  TTreeFormula *GlobalCut = new TTreeFormula("GlobalCut", globalcut, C);
  while (currentline.ReadLine(configfile)) {
    if (currentline.BeginsWith("#")) continue;
    TObjArray *tokens = currentline.Tokenize(" ");
    Int_t ntokens = tokens->GetEntries();
    if( ntokens>1 ){
      TString skey = ((TObjString*)(*tokens)[0])->GetString();
      if( skey == "macros_dir" ){
	macros_dir = ((TObjString*)(*tokens)[1])->GetString();
      }
      if( skey == "badchan_file" ){
	badchan_file = ((TObjString*)(*tokens)[1])->GetString();
      }      
      if( skey == "E_beam" ){
	E_beam = ((TObjString*)(*tokens)[1])->GetString().Atof();
      }
      if( skey == "BB_theta" ){
	th_bb = ((TObjString*)(*tokens)[1])->GetString().Atof();
	th_bb *= TMath::DegToRad(); 
      }
      if( skey == "ECAL_dist" ){
	ECAL_dist = ((TObjString*)(*tokens)[1])->GetString().Atof();
      }
      if( skey == "ECAL_posOFF" ){
	ECAL_voff = ((TObjString*)(*tokens)[1])->GetString().Atof();
	ECAL_hoff = ((TObjString*)(*tokens)[2])->GetString().Atof();
	ECAL_zoff = ((TObjString*)(*tokens)[3])->GetString().Atof();
      }
      if( skey == "nearBADchan_cut" ){
	cut_on_nearBADchan = ((TObjString*)(*tokens)[1])->GetString().Atoi();
	r_max = ((TObjString*)(*tokens)[2])->GetString().Atof();
	r2_max = pow(r_max,2);
      }            
      if( skey == "hit_threshold" ){
      	hit_threshold = ((TObjString*)(*tokens)[1])->GetString().Atof();
      }
      if( skey == "tmax_cut" ){
      	tmax_cut = ((TObjString*)(*tokens)[1])->GetString().Atof();
      }
      if( skey == "engFrac_cut" ){
      	engFrac_cut = ((TObjString*)(*tokens)[1])->GetString().Atof();
      }
      if( skey == "Min_Event_Per_Channel" ){
	Nmin = ((TObjString*)(*tokens)[1])->GetString().Atof();
      }
      if( skey == "Min_MB_Ratio" ){
	minMBratio = ((TObjString*)(*tokens)[1])->GetString().Atoi();
      }
      if( skey == "eHCAL_cut" ){
	cut_on_eHCAL = ((TObjString*)(*tokens)[1])->GetString().Atoi();
	eHCAL_cut_limit = ((TObjString*)(*tokens)[2])->GetString().Atof();
      }      
      if( skey == "EovEexp_cut" ){
	cut_on_EovEexp = ((TObjString*)(*tokens)[1])->GetString().Atoi();
	EovEexp_cut_limit = ((TObjString*)(*tokens)[2])->GetString().Atof();
      }
      if( skey == "W2_cut" ){
	cut_on_W2 = ((TObjString*)(*tokens)[1])->GetString().Atoi();
	W2_mean = ((TObjString*)(*tokens)[2])->GetString().Atof();
	W2_sigma = ((TObjString*)(*tokens)[3])->GetString().Atof();
	W2_nsigma = ((TObjString*)(*tokens)[4])->GetString().Atof();
      }
      if( skey == "dt_ADC_cut_1" ){
	dt_mean_1 = ((TObjString*)(*tokens)[1])->GetString().Atof();
	dt_sigma_1 = ((TObjString*)(*tokens)[2])->GetString().Atof();
	dt_nsigma_1 = ((TObjString*)(*tokens)[3])->GetString().Atof();
      }
      if( skey == "dt_ADC_cut_2" ){
	dt_mean_2 = ((TObjString*)(*tokens)[1])->GetString().Atof();
	dt_sigma_2 = ((TObjString*)(*tokens)[2])->GetString().Atof();
	dt_nsigma_2 = ((TObjString*)(*tokens)[3])->GetString().Atof();
      }            
      if( skey == "h_W2" ){
	h_W2_bin = ((TObjString*)(*tokens)[1])->GetString().Atoi();
	h_W2_min = ((TObjString*)(*tokens)[2])->GetString().Atof();
	h_W2_max = ((TObjString*)(*tokens)[3])->GetString().Atof();
      }
      if( skey == "h_Q2" ){
	h_Q2_bin = ((TObjString*)(*tokens)[1])->GetString().Atoi();
	h_Q2_min = ((TObjString*)(*tokens)[2])->GetString().Atof();
	h_Q2_max = ((TObjString*)(*tokens)[3])->GetString().Atof();
      }
      if( skey == "h_EovEexp" ){
	h_EovEexp_bin = ((TObjString*)(*tokens)[1])->GetString().Atoi();
	h_EovEexp_min = ((TObjString*)(*tokens)[2])->GetString().Atof();
	h_EovEexp_max = ((TObjString*)(*tokens)[3])->GetString().Atof();
      }
      if( skey == "EovEexp_fit_width" ){
	EovEexp_fit_width = ((TObjString*)(*tokens)[1])->GetString().Atof();
      }
      if( skey == "h_eHCAL" ){
	h_eHCAL_bin = ((TObjString*)(*tokens)[1])->GetString().Atoi();
	h_eHCAL_min = ((TObjString*)(*tokens)[2])->GetString().Atof();
	h_eHCAL_max = ((TObjString*)(*tokens)[3])->GetString().Atof();
      }
      if( skey == "h2_p" ){
	h2_p_bin = ((TObjString*)(*tokens)[1])->GetString().Atoi();
	h2_p_min = ((TObjString*)(*tokens)[2])->GetString().Atof();
	h2_p_max = ((TObjString*)(*tokens)[3])->GetString().Atof();
      }
      if( skey == "h2_p_coarse" ){
	h2_p_coarse_bin = ((TObjString*)(*tokens)[1])->GetString().Atoi();
	h2_p_coarse_min = ((TObjString*)(*tokens)[2])->GetString().Atof();
	h2_p_coarse_max = ((TObjString*)(*tokens)[3])->GetString().Atof();
      }
      if( skey == "h2_EovEexp" ){
	h2_EovEexp_bin = ((TObjString*)(*tokens)[1])->GetString().Atoi();
	h2_EovEexp_min = ((TObjString*)(*tokens)[2])->GetString().Atof();
	h2_EovEexp_max = ((TObjString*)(*tokens)[3])->GetString().Atof();
      }
      if( skey == "p_p_Offset" ){
	p_p_offset = ((TObjString*)(*tokens)[1])->GetString().Atof();
      }
      if( skey == "*****" ){
	break;
      }
    } 
    delete tokens;
  }

  // // ------------------------
  // // Reading Various Maps
  // // ------------------------
  // // Reading the master map of HCAL
  // LookUpTableReader HCALmap;
  // HCALmap.readCSV("maps/HCAL_r_c_x_y_cpr.csv");  
  // // Reading module IDs of HCAL blocks at the edge
  // std::unordered_set<int> edgeHCALblks = GetEdgeBlkID(HCALmap,isdebug);
  // // Reading no. of columns per row
  // std::unordered_map<int,int> ncolprow = GetNColPerRow(HCALmap,isdebug);  
  // // Reading bad HCAL module channels
  // std::vector<int> BADchanIDs = ReadNumFromTextFileToV(badchan_file.Data());
  // std::vector<std::pair<double, double>> xyBADchan = GetXYPosBADChan(HCALmap,BADchanIDs,isdebug); 
  // // ---  

  // Turning on relevant tree branches
  int maxNtr = 200; //max # of tracks expected per event
  C->SetBranchStatus("*", 0);
  // sbs.hcal branches
  C->SetBranchStatus("sbs.hcal.*", 1);
  Double_t nclusHCAL;              C->SetBranchAddress("sbs.hcal.nclus", &nclusHCAL);
  Double_t idblkHCAL;              C->SetBranchAddress("sbs.hcal.idblk", &idblkHCAL);
  Double_t rowblkHCAL;             C->SetBranchAddress("sbs.hcal.rowblk", &rowblkHCAL);
  Double_t colblkHCAL;             C->SetBranchAddress("sbs.hcal.colblk", &colblkHCAL);
  Double_t nblkHCAL;               C->SetBranchAddress("sbs.hcal.nblk", &nblkHCAL);
  Double_t atimeHCAL;              C->SetBranchAddress("sbs.hcal.atimeblk", &atimeHCAL);
  Double_t eHCAL;                  C->SetBranchAddress("sbs.hcal.e", &eHCAL);
  Double_t xHCAL;                  C->SetBranchAddress("sbs.hcal.x", &xHCAL);
  Double_t yHCAL;                  C->SetBranchAddress("sbs.hcal.y", &yHCAL);
  Double_t idclblkHCAL[maxNtr];    C->SetBranchAddress("sbs.hcal.clus_blk.id", &idclblkHCAL);
  Double_t eclblkHCAL[maxNtr];     C->SetBranchAddress("sbs.hcal.clus_blk.e", &eclblkHCAL);
  Double_t xclblkHCAL[maxNtr];     C->SetBranchAddress("sbs.hcal.clus_blk.x", &xclblkHCAL);
  Double_t yclblkHCAL[maxNtr];     C->SetBranchAddress("sbs.hcal.clus_blk.y", &yclblkHCAL);
  Double_t rowclblkHCAL[maxNtr];   C->SetBranchAddress("sbs.hcal.clus_blk.row", &rowclblkHCAL);
  Double_t colclblkHCAL[maxNtr];   C->SetBranchAddress("sbs.hcal.clus_blk.col", &colclblkHCAL);
  Double_t atimeclblkHCAL[maxNtr]; C->SetBranchAddress("sbs.hcal.clus_blk.atime", &atimeclblkHCAL);
  Double_t againblkHCAL;           C->SetBranchAddress("sbs.hcal.againblk", &againblkHCAL);  
  // earm.ecal branches
  Double_t eECAL;                C->SetBranchStatus("earm.ecal.e",1); C->SetBranchAddress("earm.ecal.e", &eECAL);
  Double_t xECAL;                C->SetBranchStatus("earm.ecal.x",1); C->SetBranchAddress("earm.ecal.x", &xECAL);
  Double_t yECAL;                C->SetBranchStatus("earm.ecal.y",1); C->SetBranchAddress("earm.ecal.y", &yECAL);
  Double_t idblkECAL;            C->SetBranchStatus("earm.ecal.idblk",1); C->SetBranchAddress("earm.ecal.idblk", &idblkECAL);  
  Double_t rowblkECAL;           C->SetBranchStatus("earm.ecal.rowblk",1); C->SetBranchAddress("earm.ecal.rowblk", &rowblkECAL);
  Double_t colblkECAL;           C->SetBranchStatus("earm.ecal.colblk",1); C->SetBranchAddress("earm.ecal.colblk", &colblkECAL);     
  Double_t atimeECAL;            C->SetBranchStatus("earm.ecal.atimeblk",1);   C->SetBranchAddress("earm.ecal.atimeblk", &atimeECAL); 
  // sbs.tr branches
  C->SetBranchStatus("sbs.tr.*", 1);
  Double_t trN;                  C->SetBranchAddress("sbs.tr.n", &trN);
  Double_t trP[maxNtr];          C->SetBranchAddress("sbs.tr.p", &trP);
  Double_t trPx[maxNtr];         C->SetBranchAddress("sbs.tr.px", &trPx);
  Double_t trPy[maxNtr];         C->SetBranchAddress("sbs.tr.py", &trPy);
  Double_t trPz[maxNtr];         C->SetBranchAddress("sbs.tr.pz", &trPz);
  Double_t trX[maxNtr];          C->SetBranchAddress("sbs.tr.x", &trX);
  Double_t trY[maxNtr];          C->SetBranchAddress("sbs.tr.y", &trY);
  Double_t trTh[maxNtr];         C->SetBranchAddress("sbs.tr.th", &trTh);
  Double_t trPh[maxNtr];         C->SetBranchAddress("sbs.tr.ph", &trPh);
  Double_t trVz[maxNtr];         C->SetBranchAddress("sbs.tr.vz", &trVz);
  Double_t trVy[maxNtr];         C->SetBranchAddress("sbs.tr.vy", &trVy);
  Double_t trTgth[maxNtr];       C->SetBranchAddress("sbs.tr.tg_th", &trTgth);
  Double_t trTgph[maxNtr];       C->SetBranchAddress("sbs.tr.tg_ph", &trTgph);
  Double_t trRx[maxNtr];         C->SetBranchAddress("sbs.tr.r_x", &trRx);
  Double_t trRy[maxNtr];         C->SetBranchAddress("sbs.tr.r_y", &trRy);
  Double_t trRth[maxNtr];        C->SetBranchAddress("sbs.tr.r_th", &trRth);
  Double_t trRph[maxNtr];        C->SetBranchAddress("sbs.tr.r_ph", &trRph);
  // heep variables
  C->SetBranchStatus("heep.*", 1);
  Double_t heep_pp_pth;          C->SetBranchAddress("heep.pp_pth", &heep_pp_pth);
  Double_t heep_pp_eth;          C->SetBranchAddress("heep.pp_eth", &heep_pp_eth);
  Double_t heep_eth_pth;         C->SetBranchAddress("heep.eth_pth", &heep_eth_pth);  
  Double_t heep_dpp;             C->SetBranchAddress("heep.dpp", &heep_dpp);
  Double_t heep_dt_ADC;          C->SetBranchAddress("heep.dt_ADC", &heep_dt_ADC);    
  Double_t heep_dxECAL;          C->SetBranchAddress("heep.dxECAL", &heep_dxECAL);
  Double_t heep_dyECAL;          C->SetBranchAddress("heep.dyECAL", &heep_dyECAL);  
  // Event info
  C->SetMakeClass(1);
  C->SetBranchStatus("fEvtHdr.*", 1);
  UInt_t rnum;                   C->SetBranchAddress("fEvtHdr.fRun", &rnum);
  UInt_t trigbits;               C->SetBranchAddress("fEvtHdr.fTrigBits", &trigbits);
  ULong64_t gevnum;              C->SetBranchAddress("fEvtHdr.fEvtNum", &gevnum);
  // turning on additional branches for the global cut
  C->SetBranchStatus("sbs.gemFT.*", 1);
  // C->SetBranchStatus("heep.*", 1);
  // C->SetBranchStatus("sbs.gem.track.ngoodhits", 1);
  // C->SetBranchStatus("sbs.gem.track.chi2ndf", 1);

  
  bool badCells[kNblks]; // Cells that have events less than Nmin
  Int_t nevents_per_blk[kNblks];  
  // Clear arrays
  memset(nevents_per_blk, 0, kNblks*sizeof(int));
  memset(badCells, 0, kNblks*sizeof(bool));

  Double_t oldADCgain[kNblks];
  for (int i=0; i<kNblks; i++) { oldADCgain[i] = -1000; }    

  gStyle->SetOptStat(0);
  TH2D *h2_eHCAL_vs_blk_raw = new TH2D("h2_eHCAL_vs_blk_raw","Raw E_clus per block",kNcols,0,kNcols,kNrows,0,kNrows);
  TH2D *h2_EovEexp_vs_blk_raw = new TH2D("h2_EovEexp_vs_blk_raw","Raw E/Eexp per block",kNcols,0,kNcols,kNrows,0,kNrows);
  TH2D *h2_EovEexp_vs_blk_raw_calib = new TH2D("h2_EovEexp_vs_blk_raw_calib","Raw E/Eexp per block (Calibrated)",kNcols,0,kNcols,kNrows,0,kNrows);
  TH2D *h2_count = new TH2D("h2_count","Count for E/Eexp per block",kNcols,0,kNcols,kNrows,0,kNrows);
  TH2D *h2_count_calib = new TH2D("h2_count_calib","Count for E/Eexp per block",kNcols,0,kNcols,kNrows,0,kNrows);
  
  // Creating output ROOT file to contain histograms
  TString outFile, outPlot, outGain, outGainR;
  std::string prefix_str = prefix.empty() ? "" : prefix + "_";
  outFile = Form("%s/hist/%s%s.root",macros_dir.Data(),prefix_str.c_str(),cfgfilebase.Data());
  outPlot = Form("%s/plots/%s%s.pdf",macros_dir.Data(),prefix_str.c_str(),cfgfilebase.Data());
  outGain = Form("%s/gain/%s%s_gainCoeff.txt",macros_dir.Data(),prefix_str.c_str(),cfgfilebase.Data());  
  outGainR = Form("%s/gain/%s%s_gainRatio.txt",macros_dir.Data(),prefix_str.c_str(),cfgfilebase.Data());  

  //std::unique_ptr<TFile> fout( TFile::Open(outFile, "RECREATE") );
  TFile *fout = new TFile(outFile, "RECREATE");
  fout->cd();
   
  TH1D *h_eHCAL = new TH1D("h_eHCAL","Best HCAL Clus. Eng. | Before",h_eHCAL_bin,h_eHCAL_min,h_eHCAL_max);
  TH1D *h_eHCAL_calib = new TH1D("h_eHCAL_calib","Best HCAL Clus. Eng. | After",h_eHCAL_bin,h_eHCAL_min,h_eHCAL_max);
  TH1D *h_EovEexp = new TH1D("h_EovEexp","E/E_exp (Before)",h_EovEexp_bin,h_EovEexp_min,h_EovEexp_max);
  TH1D *h_EovEexp_calib = new TH1D("h_EovEexp_calib","E/E_exp (After)",h_EovEexp_bin,h_EovEexp_min,h_EovEexp_max);
  TH1D *h_eSF = new TH1D("h_eSF","Energy Sampling Fraction (Before)",h_eSF_bin,h_eSF_min,h_eSF_max);
  TH1D *h_eSF_calib = new TH1D("h_eSF_calib","Energy Sampling Fraction (After)",h_eSF_bin,h_eSF_min,h_eSF_max);  
  TH1D *h_xHCAL_diff = new TH1D("h_xHCAL_diff","Vertical Pos. Diff. | Before; xHCAL - xHCAL_exp (m)",200,-0.5,0.5);
  TH1D *h_xHCAL_diff_calib = new TH1D("h_xHCAL_diff_calib","Vertical Pos. Diff. | After; xHCAL - xHCAL_exp (m)",200,-0.5,0.5);
  TH1D *h_yHCAL_diff = new TH1D("h_yHCAL_diff","Horizontal Pos. Diff. | Before; yHCAL - yHCAL_exp (m)",200,-0.5,0.5);
  TH1D *h_yHCAL_diff_calib = new TH1D("h_yHCAL_diff_calib","Horizontal Pos. Diff. | After; yHCAL - yHCAL_exp (m)",200,-0.5,0.5);
  
  TH2D *h2_EovEexp_vs_P = new TH2D("h2_EovEexp_vs_P","E/Eexp vs p | Before; p (GeV); E/Eexp",h2_p_coarse_bin,h2_p_coarse_min,h2_p_coarse_max,h2_EovEexp_bin,h2_EovEexp_min,h2_EovEexp_max);
  TProfile *h2_EovEexp_vs_P_prof = new TProfile("h2_EovEexp_vs_P_prof","E/Eexp vs P (Profile) Before",h2_p_coarse_bin,h2_p_coarse_min,h2_p_coarse_max,h_EovEexp_min,h_EovEexp_max,"S");
  CustmProfHisto(h2_EovEexp_vs_P_prof);
  TH2D *h2_EovEexp_vs_P_calib = new TH2D("h2_EovEexp_vs_P_calib","E/Eexp vs p | After; p (GeV); E/Eexp",h2_p_coarse_bin,h2_p_coarse_min,h2_p_coarse_max,h2_EovEexp_bin,h2_EovEexp_min,h2_EovEexp_max);
  TProfile *h2_EovEexp_vs_P_prof_calib = new TProfile("h2_EovEexp_vs_P_prof_calib","E/Eexp vs P (Profile) After",h2_p_coarse_bin,h2_p_coarse_min,h2_p_coarse_max,h_EovEexp_min,h_EovEexp_max,"S");
  CustmProfHisto(h2_EovEexp_vs_P_prof_calib);
  
  TH2D *h2_EovEexp_vs_blk = new TH2D("h2_EovEexp_vs_blk","E/Eexp per block | Before",kNcols,0,kNcols,kNrows,0,kNrows);
  TH2D *h2_EovEexp_vs_blk_calib = new TH2D("h2_EovEexp_vs_blk_calib","E/Eexp per block | After",kNcols,0,kNcols,kNrows,0,kNrows);  

  TH2D *h2_EovEexp_vs_xHCALexp = new TH2D("h2_EovEexp_vs_xHCALexp","E/Eexp vs xHCAL Expected | Before",200,-1.7,1.7,h_EovEexp_bin,h_EovEexp_min,h_EovEexp_max);
  TH2D *h2_EovEexp_vs_xHCALexp_calib = new TH2D("h2_EovEexp_vs_xHCALexp_calib","E/Eexp vs xHCAL Expected | After",200,-1.7,1.7,h_EovEexp_bin,h_EovEexp_min,h_EovEexp_max);
  TH2D *h2_EovEexp_vs_yHCALexp = new TH2D("h2_EovEexp_vs_yHCALexp","E/Eexp vs yHCAL Expected | Before",200,-0.8,0.8,h_EovEexp_bin,h_EovEexp_min,h_EovEexp_max);
  TH2D *h2_EovEexp_vs_yHCALexp_calib = new TH2D("h2_EovEexp_vs_yHCALexp_calib","E/Eexp vs yHCAL Expected | After",200,-0.8,0.8,h_EovEexp_bin,h_EovEexp_min,h_EovEexp_max);
  TH2D *h2_eHCAL_vs_xHCALexp = new TH2D("h2_eHCAL_vs_xHCALexp","HCAL energy vs xHCAL Expected | Before",200,-1.7,1.7,h_eHCAL_bin,h_eHCAL_min,h_eHCAL_max);
  TH2D *h2_eHCAL_vs_xHCALexp_calib = new TH2D("h2_eHCAL_vs_xHCALexp_calib","HCAL energy vs xHCAL Expected | After",200,-1.7,1.7,h_eHCAL_bin,h_eHCAL_min,h_eHCAL_max);
  TH2D *h2_eHCAL_vs_yHCALexp = new TH2D("h2_eHCAL_vs_yHCALexp","HCAL energy vs yHCAL Expected | Before",200,-0.8,0.8,h_eHCAL_bin,h_eHCAL_min,h_eHCAL_max);
  TH2D *h2_eHCAL_vs_yHCALexp_calib = new TH2D("h2_eHCAL_vs_yHCALexp_calib","HCAL energy vs yHCAL Expected | After",200,-0.8,0.8,h_eHCAL_bin,h_eHCAL_min,h_eHCAL_max);  

  TH2D *h2_nev_per_blk = new TH2D("h2_nev_per_blk","Number of Good Events per Block;HCAL cols;HCAL rows",kNcols,0,kNcols,kNrows,0,kNrows);
  // TH2D *h2_nev_per_blk_bot = new TH2D("h2_nev_per_blk_bot","# good ev per blk (Bottom Half);HCAL cols;HCAL rows",kNcols,0,kNcols,33,0,33);
  // TH2D *h2_nev_per_blk_top = new TH2D("h2_nev_per_blk_top","# good ev per blk (Top Half);HCAL cols;HCAL rows",kNcols,0,kNcols,35,34,69);  
  
  // defining output ROOT tree (Set max size to 4GB)
  TTree *Tout = new TTree("Tout", cfgfilebase.Data()); 
  Tout->SetMaxTreeSize(4000000000LL);
  //
  UInt_t    T_rnum;     Tout->Branch("rnum", &T_rnum, "rnum/i");  // The run number for each set of data. This is important because run numbers are not always continuous.
  ULong64_t T_gevnum;   Tout->Branch("gevnum", &T_gevnum, "gevnum/l");  // Global event number    
  //
  bool T_indtADC;        Tout->Branch("indtADC", &T_indtADC, "indtADC/O");  
  bool T_isHCALedge;     Tout->Branch("isHCALedge", &T_isHCALedge, "isHCALedge/O");
  bool T_isNearBadChan;  Tout->Branch("isNearBadChan", &T_isNearBadChan, "isNearBadChan/O");  
  //
  Double_t T_heep_dpp;    Tout->Branch("heep_dpp", &T_heep_dpp, "heep_dpp/D");
  Double_t T_heep_pp_pth; Tout->Branch("heep_pp_pth", &T_heep_pp_pth, "heep_pp_pth/D");
  Double_t T_heep_pp_eth; Tout->Branch("heep_pp_eth", &T_heep_pp_eth, "heep_pp_eth/D");
  Double_t T_heep_dt_ADC; Tout->Branch("heep_dt_ADC", &T_heep_dt_ADC, "heep_dt_ADC/D");  
  Double_t T_heep_dxECAL; Tout->Branch("heep_dxECAL", &T_heep_dxECAL, "heep_dxECAL/D");
  Double_t T_heep_dyECAL; Tout->Branch("heep_dyECAL", &T_heep_dyECAL, "heep_dyECAL/D");    
  Double_t T_heep_eth_pth;Tout->Branch("heep_eth_pth", &T_heep_eth_pth, "heep_eth_pth/D");
  //
  Double_t T_E_beam;     Tout->Branch("E_beam", &T_E_beam, "E_beam/D");
  Double_t T_p_p;        Tout->Branch("p_p", &T_p_p, "p_p/D");
  Double_t T_p_px;       Tout->Branch("p_px", &T_p_px, "p_px/D");
  Double_t T_p_py;       Tout->Branch("p_py", &T_p_py, "p_py/D");
  Double_t T_p_pz;       Tout->Branch("p_pz", &T_p_pz, "p_pz/D");  
  Double_t T_E_p;        Tout->Branch("E_p", &T_E_p, "E_p/D");
  Double_t T_KE_p;       Tout->Branch("KE_p", &T_KE_p, "KE_p/D");
  Double_t T_KE_p_frac;  Tout->Branch("KE_p_frac", &T_KE_p_frac, "KE_p_frac/D");    
  Double_t T_th_p;       Tout->Branch("th_p", &T_th_p, "th_p/D");
  Double_t T_ph_p;       Tout->Branch("ph_p", &T_ph_p, "ph_p/D");
  Double_t T_p_e;        Tout->Branch("p_e", &T_p_e, "p_e/D");
  Double_t T_th_e;       Tout->Branch("th_e", &T_th_e, "th_e/D");
  Double_t T_ph_e;       Tout->Branch("ph_e", &T_ph_e, "ph_e/D");
  Double_t T_p_e_th;     Tout->Branch("p_e_th", &T_p_e_th, "p_e_th/D");
  //
  Double_t T_nu;         Tout->Branch("nu", &T_nu, "nu/D"); 
  Double_t T_W2;         Tout->Branch("W2", &T_W2, "W2/D"); 
  Double_t T_Q2;         Tout->Branch("Q2", &T_Q2, "Q2/D");
  Double_t T_dxECAL;     Tout->Branch("dxECAL", &T_dxECAL, "dxECAL/D");
  Double_t T_dyECAL;     Tout->Branch("dyECAL", &T_dyECAL, "dyECAL/D");    
  Double_t T_xECAL_exp;  Tout->Branch("xECAL_exp", &T_xECAL_exp, "xECAL_exp/D"); 
  Double_t T_yECAL_exp;  Tout->Branch("yECAL_exp", &T_yECAL_exp, "yECAL_exp/D");
  //
  // Double_t T_nu_4vec;    if (enable_4vec) Tout->Branch("nu_4vec", &T_nu_4vec, "nu_4vec/D"); 
  // Double_t T_W2_4vec;    if (enable_4vec) Tout->Branch("W2_4vec", &T_W2_4vec, "W2_4vec/D"); 
  // Double_t T_Q2_4vec;    if (enable_4vec) Tout->Branch("Q2_4vec", &T_Q2_4vec, "Q2_4vec/D");
  // Double_t T_xECAL_exp_4vec; if (enable_4vec) Tout->Branch("xECAL_exp_4vec", &T_xECAL_exp_4vec, "xECAL_exp_4vec/D"); 
  // Double_t T_yECAL_exp_4vec; if (enable_4vec) Tout->Branch("yECAL_exp_4vec", &T_yECAL_exp_4vec, "yECAL_exp_4vec/D");
  // HCAL related
  Double_t T_eHCAL;      Tout->Branch("eHCAL", &T_eHCAL, "eHCAL/D");
  Double_t T_xHCAL;      Tout->Branch("xHCAL", &T_xHCAL, "xHCAL/D"); 
  Double_t T_yHCAL;      Tout->Branch("yHCAL", &T_yHCAL, "yHCAL/D"); 
  Double_t T_rowblkHCAL; Tout->Branch("rowblkHCAL", &T_rowblkHCAL, "rowblkHCAL/D"); 
  Double_t T_colblkHCAL; Tout->Branch("colblkHCAL", &T_colblkHCAL, "colblkHCAL/D"); 
  Double_t T_idblkHCAL;  Tout->Branch("idblkHCAL", &T_idblkHCAL, "idblkHCAL/D"); 
  Double_t T_atimeHCAL;  Tout->Branch("atimeHCAL", &T_atimeHCAL, "atimeHCAL/D"); 
  Double_t T_EovEexp;    Tout->Branch("EovEexp", &T_EovEexp, "EovEexp/D");
  // ECAL related
  Double_t T_eECAL;      Tout->Branch("eECAL", &T_eECAL, "eECAL/D"); 
  Double_t T_xECAL;      Tout->Branch("xECAL", &T_xECAL, "xECAL/D"); 
  Double_t T_yECAL;      Tout->Branch("yECAL", &T_yECAL, "yECAL/D"); 
  Double_t T_rowblkECAL; Tout->Branch("rowblkECAL", &T_rowblkECAL, "rowblkECAL/D"); 
  Double_t T_colblkECAL; Tout->Branch("colblkECAL", &T_colblkECAL, "colblkECAL/D"); 
  Double_t T_idblkECAL;  Tout->Branch("idblkECAL", &T_idblkECAL, "idblkECAL/D"); 
  Double_t T_atimeECAL;  Tout->Branch("atimeECAL", &T_atimeECAL, "atimeECAL/D");   
  
  // calculating ECAL co-ordinates
  TVector3 ECAL_zaxis(sin(th_bb),0,cos(th_bb)); // use BB angle to calculate the center of Ecal
  TVector3 ECAL_xaxis(0,-1,0); 
  TVector3 ECAL_yaxis = ECAL_zaxis.Cross(ECAL_xaxis).Unit();
  // Define the center of Ecal in 3D space
  TVector3 ECAL_origin = (ECAL_dist+ECAL_zoff)*ECAL_zaxis + ECAL_voff*ECAL_xaxis + ECAL_hoff*ECAL_yaxis; 
  
  std::cout << std::endl;
  Long64_t Ngoodevs=0, Nelasevs=0; 
  Long64_t Nevents = C->GetEntries(), nevent=0; UInt_t runnum=0; 
  Double_t timekeeper=0., timeremains=0.;
  Int_t treenum=0, currenttreenum=0, itrrun=0;
  std::vector<std::string> lrnum;    // list of run numbers  
  while(C->GetEntry(nevent++)) {
    // Calculating remaining time 
    sw2->Stop();
    timekeeper += sw2->RealTime();
    if (nevent % 5000 == 0 && nevent != 0) 
      timeremains = timekeeper * (double(Nevents) / double(nevent) - 1.); 
    sw2->Reset();
    sw2->Continue();

    if(nevent % 100 == 0) std::cout << nevent << "/" << Nevents  << ", " << int(timeremains/60.) << "m \r";;
    std::cout.flush();
    // ------

    // get old gain coefficients (do it before applying global cut)
    oldADCgain[int(idblkHCAL)-1] = againblkHCAL;
        
    // apply global cuts efficiently (AJRP method)
    currenttreenum = C->GetTreeNumber();
    if (nevent == 1 || currenttreenum != treenum) {
      treenum = currenttreenum;
      GlobalCut->UpdateFormulaLeaves();

      // track change of runnum
      if (nevent == 1 || rnum != runnum) {
	runnum = rnum; itrrun++;
	lrnum.push_back(to_string(rnum));
      }
    } 
    bool passedgCut = GlobalCut->EvalInstance(0) != 0;   
    if (passedgCut) {
      memset(A, 0, kNblks*sizeof(double));      

      // HCAL active area cut
      isHCALedge = rowblkHCAL == 0 || rowblkHCAL == 23 || colblkHCAL == 0 || colblkHCAL == 11; 
      if (isHCALedge) continue;

      // Applying dt cut
      if (rnum<2660) {
	dt_mean = dt_mean_1;
	dt_sigma = dt_sigma_1;
	dt_nsigma = dt_nsigma_1;
      } else {
	dt_mean = dt_mean_2;
	dt_sigma = dt_sigma_2;
	dt_nsigma = dt_nsigma_2;	
      }
      indtADC = abs(heep_dt_ADC-dt_mean)<dt_sigma*dt_nsigma;
      if (!indtADC) continue;
      
      // keeping count of events passing globcal cuts
      Ngoodevs++;
      
      TVector3 vertex(0,0,trVz[0]);
      
      // scattered proton kinematics (reconstructed)
      p_p = heep_pp_pth;  //trP[0]; // [04:16:25] Calculating based on angles until momentum is better calibrated.
      p_px = trPx[0];
      p_py = trPy[0];
      p_pz = trPz[0];    
      E_p = sqrt(p_p*p_p + Mp*Mp);
      KE_p = E_p - Mp;
      // Applying HCAL sampling fraction
      KE_p_frac = KE_p*HCAL_sf; 
	
      th_p = TMath::ACos(trPz[0]/trP[0]);
      ph_p = atan2(trPy[0],-trPx[0]); 
      
      // scattered electron kinematics (from p kinematics)
      p_e = E_beam + Mp - E_p;
      th_e = 2.0 * TMath::ATan((Mp/(Mp+E_beam))*(1./TMath::Tan(th_p)));
      ph_e = - ph_p;
      p_e_th = E_beam/(1. + (E_beam/Mp)*(1.0 - cos(th_e)));

      if (recon_type==1) p_e_CURR = p_e; //=> Use momentum recon method
      else p_e_CURR = p_e_th; //=> Use angle recon method (Preferred)

      // kinematic variables
      nu = E_beam - p_e_CURR; 
      Q2 = 2.*E_beam*p_e_CURR*(1.-cos(th_e)); 
      W2 = Mp*Mp + 2.*Mp*nu - Q2;
      TVector3 Pef_3vec(p_e_CURR*sin(th_e)*cos(ph_e), p_e_CURR*sin(th_e)*sin(ph_e), p_e_CURR*cos(th_e));    
      TVector3 Pefhat = Pef_3vec.Unit();
      // calculating expected hit positions on ECAL (angles only method)
      Double_t sintersect = (ECAL_origin - vertex).Dot(ECAL_zaxis) / (Pefhat.Dot(ECAL_zaxis));
      TVector3 ECAL_intersect = vertex + sintersect*Pefhat; 
      T_xECAL_exp = (ECAL_intersect - ECAL_origin).Dot(ECAL_xaxis);
      T_yECAL_exp = (ECAL_intersect - ECAL_origin).Dot(ECAL_yaxis);
      T_dxECAL = xECAL - T_xECAL_exp;
      T_dyECAL = yECAL - T_yECAL_exp;    

      T_rnum = rnum;
      T_gevnum = gevnum;      
      
      T_indtADC = indtADC;      
      T_isHCALedge = isHCALedge;
      T_isNearBadChan = isNearBadChan;

      T_E_beam = E_beam;
      
      T_heep_pp_pth = heep_pp_pth;
      T_heep_pp_eth = heep_pp_eth;
      T_heep_eth_pth = heep_eth_pth;
      T_heep_dt_ADC = heep_dt_ADC;            
      T_heep_dpp = heep_dpp;
      T_heep_dxECAL = heep_dxECAL;
      T_heep_dyECAL = heep_dyECAL;      
      
      T_p_p = p_p;
      T_p_px = p_px;
      T_p_py = p_py;
      T_p_pz = p_pz;    
      T_E_p = E_p;
      T_KE_p = KE_p;
      T_KE_p_frac = KE_p_frac;
      
      T_th_p = th_p;
      T_ph_p = ph_p; 

      // scattered electron kinematics (from p kinematics)
      T_p_e = p_e;
      T_th_e = th_e;
      T_ph_e = ph_e;
      T_p_e_th = p_e_th;

      // kinematic variables
      T_nu = nu;
      T_Q2 = Q2;
      T_W2 = W2;
      
      // Filling HCAL related branches
      T_eHCAL = eHCAL;
      T_xHCAL = xHCAL;
      T_yHCAL = yHCAL;
      T_rowblkHCAL = rowblkHCAL;
      T_colblkHCAL = colblkHCAL;
      T_idblkHCAL = idblkHCAL;
      T_atimeHCAL = atimeHCAL;

      // Filling ECAL related branches
      T_eECAL = eECAL;
      T_xECAL = xECAL;
      T_yECAL = yECAL;
      T_rowblkECAL = rowblkECAL;
      T_colblkECAL = colblkECAL;
      T_idblkECAL = idblkECAL;
      T_atimeECAL = atimeECAL;      
      
      // Loop over all the blocks in main cluster and fill in A's
      Double_t eHCALcl = 0.;
      for(Int_t blk=0; blk<nblkHCAL; blk++){
	Int_t blkID = int(idclblkHCAL[blk]) - 1;	
	// if (eclblkHCAL[blk]>sh_hit_threshold) {
	// Double_t shtdiff = shClBlkAtime[blk]-shClBlkAtime[0];
	// Double_t shengFrac = eclblkHCAL[blk]/eclblkHCAL[0];
	// if (fabs(shtdiff)<sh_tmax_cut && shengFrac>=sh_engFrac_cut) {
	Double_t eclblkHCAL_i = eclblkHCAL[blk];
	A[blkID] += eclblkHCAL_i;
	eHCALcl += eclblkHCAL_i;
	// filling cluster level histos
	// if (blk!=0) {
	//   h_HCALcltdiff->Fill(shtdiff);
	//   h2_HCALtdiff_vs_engFrac->Fill(shengFrac,shtdiff);
	// }
	// }
	// }
	h2_nev_per_blk->Fill(colclblkHCAL[blk],rowclblkHCAL[blk],1.);
	// if (blkID<=791) h2_nev_per_blk_bot->Fill(colclblkHCAL[blk],rowclblkHCAL[blk],1.);
	// else h2_nev_per_blk_top->Fill(colclblkHCAL[blk],rowclblkHCAL[blk],1.);
	nevents_per_blk[blkID]++; 
      }      
      Double_t EovEexp = eHCALcl/KE_p_frac;
      Double_t eSF = eHCALcl/KE_p;      
      
      T_EovEexp = EovEexp; // ********* need to move? After elastic cut formation?
      
      // filling diagnostic histos
      h_eSF->Fill(eSF);
      h_eHCAL->Fill(eHCALcl);
      h_EovEexp->Fill(EovEexp);

      h_xHCAL_diff->Fill(dx);
      h_yHCAL_diff->Fill(dy);      
     
      // E/Eexp vs. p
      h2_EovEexp_vs_P->Fill(KE_p_frac,EovEexp);
      h2_EovEexp_vs_P_prof->Fill(KE_p_frac, EovEexp, 1.);

      // E/Eexp vs. blk
      h2_EovEexp_vs_blk_raw->Fill(colblkHCAL,rowblkHCAL,EovEexp);
      h2_count->Fill(colblkHCAL,rowblkHCAL,1.);
      
      // // Correlation with expected track position
      // h2_eHCAL_vs_xHCALexp->Fill(T_xECAL_exp, eHCALcl);
      // h2_eHCAL_vs_yHCALexp->Fill(T_yHCAL_exp, eHCALcl);      
      // h2_EovEexp_vs_xHCALexp->Fill(T_xHCAL_exp, EovEexp);
      // h2_EovEexp_vs_yHCALexp->Fill(T_yHCAL_exp, EovEexp);      

      // E/Eexp vs. rnum

      // Clus variables vs. rnum
      
      // Lets form the matrix of linear equations
      BuildMatrix(A,B,M,KE_p_frac);

      Tout->Fill();
    } // global cut
  } // event loop
  std::cout << "\n\n";
  h2_EovEexp_vs_blk->Divide(h2_EovEexp_vs_blk_raw, h2_count);

  // Let's customize the histogram ranges
  h2_EovEexp_vs_blk->GetZaxis()->SetRangeUser(0.6,1.2);  

  // B.Print();  
  // M.Print();

  // ////////////////////////////////////////////////////
  // // Time to calculate and report gain coefficients //
  // ////////////////////////////////////////////////////

  TH1D *h_nevent_blk = new TH1D("h_nevent_blk", "No. of Good Events; HCAL Blocks", kNblks, 0, kNblks);
  TH1D *h_coeff_Ratio_blk = new TH1D("h_coeff_Ratio_blk", "Ratio of Gain Coefficients(new/old); HCAL Blocks", kNblks, 0, kNblks);
  TH1D *h_coeff_blk = new TH1D("h_coeff_blk", "ADC Gain Coefficients(GeV/pC); HCAL Blocks", kNblks, 0, kNblks);
  TH1D *h_coeff_blk_old = new TH1D("h_old_coeff_blk", "Old ADC Gain Coefficients(GeV/pC); HCAL Blocks", kNblks, 0, kNblks);
  TH2D *h2_coeff_detView_bot = new TH2D("h2_coeff_detView_bot", "New ADC Gain Coefficients (Bottom 23 Rows)", kNcols, 0, kNcols, 23, 0, 23);
  TH2D *h2_coeff_detView_old_bot = new TH2D("h2_coeff_detView_old_bot", "Old ADC Gain Coefficients (Bottom 23 Rows)", kNcols, 0, kNcols, 23, 0, 23);    
  TH2D *h2_coeff_detView_mid = new TH2D("h2_coeff_detView_mid", "New ADC Gain Coefficients (Middle 23 Rows)", kNcols, 0, kNcols, 23, 23, 46);
  TH2D *h2_coeff_detView_old_mid = new TH2D("h2_coeff_detView_old_mid", "Old ADC Gain Coefficients (Middle 23 Rows)", kNcols, 0, kNcols, 23, 23, 46);
  TH2D *h2_coeff_detView_top = new TH2D("h2_coeff_detView_top", "New ADC Gain Coefficients (Top 23 Rows)", kNcols, 0, kNcols, 23, 46, 69);
  TH2D *h2_coeff_detView_old_top = new TH2D("h2_coeff_detView_old_top", "Old ADC Gain Coefficients (Top 23 Rows)", kNcols, 0, kNcols, 23, 46, 69);
  
  // Leave the bad channels out of the calculation
  for(Int_t j = 0; j<kNblks; j++){
    badCells[j]=false;
    if (nevents_per_blk[j] < Nmin || M(j,j) < minMBratio*B(j)) {
      B(j) = 1.;
      M(j, j) = 1.;
      for(Int_t k = 0; k<kNblks; k++){
	if(k!=j){
	  M(j, k) = 0.;
	  M(k, j) = 0.;
	}
      }
      badCells[j]=true;
    }
  }  
  
  // Getting coefficients (rather ratios)
  M_inv = M.Invert();
  CoeffR = M_inv*B;

  // Filing histos and files with coefficients
  int iblk=0;
  ofstream outGain_data, outGainR_data;
  outGain_data.open(outGain);
  outGainR_data.open(outGainR);
  Double_t newADCgratio[kNblks];
  for (int i=0; i<kNblks; i++) { newADCgratio[i] = -1000; }  
  for (int irow=0; irow<kNrows; irow++) {
    for (int icol=0; icol<kNcols; icol++) {      
      Double_t oldCoeff = oldADCgain[iblk];
      Double_t gainRatio = CoeffR(iblk);
      if(!badCells[iblk]) {
	h_nevent_blk->Fill(iblk, nevents_per_blk[iblk]);
	h_coeff_Ratio_blk->Fill(iblk, gainRatio);
	h_coeff_blk->Fill(iblk, gainRatio*oldCoeff);
	if (irow<23) {
	  h2_coeff_detView_bot->Fill(icol, irow, gainRatio*oldCoeff);
	  h2_coeff_detView_old_bot->Fill(icol, irow, oldCoeff);
	} else if (irow<46) {
	  h2_coeff_detView_mid->Fill(icol, irow, gainRatio*oldCoeff);
	  h2_coeff_detView_old_mid->Fill(icol, irow, oldCoeff);
	} else {
	  h2_coeff_detView_top->Fill(icol, irow, gainRatio*oldCoeff);
	  h2_coeff_detView_old_top->Fill(icol, irow, oldCoeff);
	}

	// std::cout << gainRatio << "  ";
	outGain_data << gainRatio * oldCoeff << " ";
	outGainR_data << gainRatio << " ";	
	newADCgratio[iblk] = gainRatio;
      }else {
	gainRatio = 1.; // (Important)
	h_nevent_blk->Fill(iblk,nevents_per_blk[iblk]);
	h_coeff_Ratio_blk->Fill(iblk, gainRatio);
	h_coeff_blk->Fill(iblk, gainRatio*oldCoeff);
	if (irow<23) {
	  h2_coeff_detView_bot->Fill(icol, irow, gainRatio*oldCoeff);
	  h2_coeff_detView_old_bot->Fill(icol, irow, oldCoeff);
	} else if (irow<46) {
	  h2_coeff_detView_mid->Fill(icol, irow, gainRatio*oldCoeff);
	  h2_coeff_detView_old_mid->Fill(icol, irow, oldCoeff);
	} else {
	  h2_coeff_detView_top->Fill(icol, irow, gainRatio*oldCoeff);
	  h2_coeff_detView_old_top->Fill(icol, irow, oldCoeff);
	}

	// std::cout << gainRatio << "  ";
	outGain_data << gainRatio * oldCoeff << " ";
	outGainR_data << gainRatio << " ";	
	newADCgratio[iblk] = gainRatio;
      }

      iblk += 1;
    } //col
    // std::cout << std::endl;
    outGain_data << std::endl;
    outGainR_data << std::endl;    
  } //row
  std::cout << std::endl;

  // //////////////////////////////////////////////////////////////////////
  // // 2nd Loop over all events to check the performance of calibration //
  // //////////////////////////////////////////////////////////////////////

  // add branches to Tout to store values after calibration
  Double_t T_eHCAL_calib;  TBranch *T_eHCAL_c = Tout->Branch("eHCAL_calib", &T_eHCAL_calib, "eHCAL_calib/D");  
  Double_t T_xHCAL_calib;  TBranch *T_xHCAL_c = Tout->Branch("xHCAL_calib", &T_xHCAL_calib, "xHCAL_calib/D");
  Double_t T_yHCAL_calib;  TBranch *T_yHCAL_c = Tout->Branch("yHCAL_calib", &T_yHCAL_calib, "yHCAL_calib/D");
  Double_t T_EovEexp_calib;   TBranch *T_EovEexp_c = Tout->Branch("EovEexp_calib", &T_EovEexp_calib, "EovEexp_calib/D");  

  itrrun=0; runnum=0;     
  Nevents = C->GetEntries(), nevent=0;
  std::cout << "Looping over events again to check calibration.." << std::endl; 
  while(C->GetEntry(nevent++)) {
    // Calculating remaining time 
    sw2->Stop();
    timekeeper += sw2->RealTime();
    if (nevent % 25000 == 0 && nevent != 0) 
      timeremains = timekeeper * (double(Nevents) / double(nevent) - 1.); 
    sw2->Reset();
    sw2->Continue();

    if(nevent % 100 == 0) std::cout << nevent << "/" << Nevents  << ", " << int(timeremains/60.) << "m \r";;
    std::cout.flush();
    // ------

    // apply global cuts efficiently (AJRP method)
    currenttreenum = C->GetTreeNumber();
    if (nevent == 1 || currenttreenum != treenum) {
      treenum = currenttreenum;
      GlobalCut->UpdateFormulaLeaves();

      // track change of runnum
      if (nevent == 1 || rnum != runnum) {
	runnum = rnum; itrrun++;
      }
    } 
    bool passedgCut = GlobalCut->EvalInstance(0) != 0;   
    if (passedgCut) {

      // HCAL active area cut
      isHCALedge = rowblkHCAL == 0 || rowblkHCAL == 23 || colblkHCAL == 0 || colblkHCAL == 11; 
      if (isHCALedge) continue;
      
      // Applying dt cut
      if (rnum<2660) {
	dt_mean = dt_mean_1;
	dt_sigma = dt_sigma_1;
	dt_nsigma = dt_nsigma_1;
      } else {
	dt_mean = dt_mean_2;
	dt_sigma = dt_sigma_2;
	dt_nsigma = dt_nsigma_2;	
      }
      indtADC = abs(heep_dt_ADC-dt_mean)<dt_sigma*dt_nsigma;
      if (!indtADC) continue;      

      TVector3 vertex(0,0,trVz[0]);
      
      // scattered proton kinematics (reconstructed)
      p_p = trP[0];
      p_px = trPx[0];
      p_py = trPy[0];
      p_pz = trPz[0];
      E_p = sqrt(p_p*p_p + Mp*Mp);
      KE_p = E_p - Mp;
      // Applying HCAL sampling fraction
      KE_p_frac = KE_p*HCAL_sf;      
      
      th_p = TMath::ACos(trPz[0]/trP[0]);
      ph_p = atan2(trPy[0],-trPx[0]);

      // scattered electron kinematics (from p kinematics)
      p_e = E_beam + Mp - E_p;
      th_e = 2.0 * TMath::ATan((Mp/(Mp+E_beam))*(1./TMath::Tan(th_p)));
      ph_e = - ph_p;
      p_e_th = E_beam/(1. + (E_beam/Mp)*(1.0 - cos(th_e)));

      if (recon_type==1) p_e_CURR = p_e; //=> Use momentum recon method
      else p_e_CURR = p_e_th; //=> Use angle recon method (Preferred)      

      // kinematic variables
      nu = E_beam - p_e_CURR; 
      Q2 = 2.*E_beam*p_e_CURR*(1.-cos(th_e)); 
      W2 = Mp*Mp + 2.*Mp*nu - Q2;
      TVector3 Pef_3vec(p_e_CURR*sin(th_e)*cos(ph_e), p_e_CURR*sin(th_e)*sin(ph_e), p_e_CURR*cos(th_e));    
      TVector3 Pefhat = Pef_3vec.Unit();
      // calculating expected hit positions on HCAL (angles only method)
      Double_t sintersect = (ECAL_origin - vertex).Dot(ECAL_zaxis) / (Pefhat.Dot(ECAL_zaxis));
      TVector3 ECAL_intersect = vertex + sintersect*Pefhat; 
      Double_t xECAL_exp = (ECAL_intersect - ECAL_origin).Dot(ECAL_xaxis);
      Double_t yECAL_exp = (ECAL_intersect - ECAL_origin).Dot(ECAL_yaxis);

      // Calculating calibrated energy and cluster centroid
      Double_t eHCALcl_calib = 0., xHCALcl_calib = 0., yHCALcl_calib = 0., eclblkHCAL_calib_HE = 0.;   
      for(Int_t blk=0; blk<nblkHCAL; blk++){
	Int_t blkID = int(idclblkHCAL[blk] - 1);	
	// calculating the updated cluster centroid
	xHCALcl_calib = (xHCALcl_calib*eHCALcl_calib + xclblkHCAL[blk]*eclblkHCAL[blk]) / (eHCALcl_calib+eclblkHCAL[blk]);
	yHCALcl_calib = (yHCALcl_calib*eHCALcl_calib + yclblkHCAL[blk]*eclblkHCAL[blk]) / (eHCALcl_calib+eclblkHCAL[blk]);	

	if (blk==0) eclblkHCAL_calib_HE = eclblkHCAL[blk] * newADCgratio[blkID];
	Double_t eclblkHCAL_calib = eclblkHCAL[blk] * newADCgratio[blkID];
	
	// if (eclblkHCAL[blk]>sh_hit_threshold) {
	// Double_t shtdiff = shClBlkAtime[blk]-shClBlkAtime[0];
	// Double_t shengFrac = eclblkHCAL[blk]/eclblkHCAL[0];
	// if (fabs(shtdiff)<sh_tmax_cut && shengFrac>=sh_engFrac_cut) {
	eHCALcl_calib += eclblkHCAL_calib;
	// filling cluster level histos
	// if (blk!=0) {
	//   h_HCALcltdiff->Fill(shtdiff);
	//   h2_HCALtdiff_vs_engFrac->Fill(shengFrac,shtdiff);
	// }
	// }
	// }
      }
      Double_t EovEexp_calib = eHCALcl_calib/KE_p_frac;
      Double_t eSF_calib = eHCALcl_calib/KE_p;

      // filling output tree with calibrated entries
      T_eHCAL_calib = eHCALcl_calib; T_eHCAL_c->Fill();
      T_xHCAL_calib = xHCALcl_calib; T_xHCAL_c->Fill();
      T_yHCAL_calib = yHCALcl_calib; T_yHCAL_c->Fill();      
      T_EovEexp_calib = EovEexp_calib; T_EovEexp_c->Fill();      

      // filling diagnostic histos
      h_eSF_calib->Fill(eSF_calib);
      h_eHCAL_calib->Fill(eHCALcl_calib);
      h_EovEexp_calib->Fill(EovEexp_calib);

      // h_xHCAL_diff_calib->Fill(xHCALcl_calib - xHCAL_exp);
      // h_yHCAL_diff_calib->Fill(yHCALcl_calib - yHCAL_exp);      
     
      // E/Eexp vs. p
      h2_EovEexp_vs_P_calib->Fill(p_e_CURR,EovEexp_calib);
      h2_EovEexp_vs_P_prof_calib->Fill(p_e_CURR, EovEexp_calib, 1.);

      // E/Eexp vs. blk
      h2_EovEexp_vs_blk_raw_calib->Fill(colblkHCAL,rowblkHCAL,EovEexp_calib);
      h2_count_calib->Fill(colblkHCAL,rowblkHCAL,1.);
      
      // Correlation with expected track position
      // h2_eHCAL_vs_xHCALexp_calib->Fill(xHCAL_exp, eHCALcl_calib);
      // h2_eHCAL_vs_yHCALexp_calib->Fill(yHCAL_exp, eHCALcl_calib);      
      // h2_EovEexp_vs_xHCALexp_calib->Fill(xHCAL_exp, EovEexp_calib);
      // h2_EovEexp_vs_yHCALexp_calib->Fill(yHCAL_exp, EovEexp_calib);      

      // E/Eexp vs. rnum

      // Clus variables vs. rnum
          
    }
  }
  std::cout << "\n\n";
  h2_EovEexp_vs_blk_calib->Divide(h2_EovEexp_vs_blk_raw_calib, h2_count_calib);

  // Let's customize the histogram ranges
  h2_EovEexp_vs_blk_calib->GetZaxis()->SetRangeUser(0.6,1.2);    

  /////////////////////////////////
  // Generating diagnostic plots //
  /////////////////////////////////
  /**** Global settings ****/
  //gStyle->SetPalette(kRainBow);

  /**** Canvas 1 (E/Eexp) ****/
  TCanvas *c1 = new TCanvas("c1","E/Eexp",1500,1000);
  c1->Divide(2,2);
  c1->cd(1); //
  gPad->SetGridx();
  Double_t param[3], param_bc[3], sigerr, sigerr_bc;
  Int_t maxBin_bc = h_EovEexp->GetMaximumBin();
  Double_t binW_bc = h_EovEexp->GetBinWidth(maxBin_bc), norm_bc = h_EovEexp->GetMaximum();
  Double_t mean_bc = h_EovEexp->GetMean(), stdev_bc = h_EovEexp->GetStdDev();
  Double_t lower_lim_bc = h_EovEexp_min + maxBin_bc*binW_bc - EovEexp_fit_width*stdev_bc;
  Double_t upper_lim_bc = h_EovEexp_min + maxBin_bc*binW_bc + EovEexp_fit_width*stdev_bc; 
  TF1* fitg_bc = new TF1("fitg_bc","gaus",h_EovEexp_min,h_EovEexp_max);
  fitg_bc->SetRange(lower_lim_bc,upper_lim_bc);
  fitg_bc->SetParameters(norm_bc,mean_bc,stdev_bc);
  fitg_bc->SetLineWidth(2); fitg_bc->SetLineColor(1);
  h_EovEexp->Fit(fitg_bc,"QR"); fitg_bc->GetParameters(param_bc); sigerr_bc = fitg_bc->GetParError(2);
  h_EovEexp->SetLineWidth(2); h_EovEexp->SetLineColor(kRed); h_EovEexp->SetFillColorAlpha(kRed,0.4);
  Int_t maxBin = h_EovEexp_calib->GetMaximumBin();
  Double_t binW = h_EovEexp_calib->GetBinWidth(maxBin), norm = h_EovEexp_calib->GetMaximum();
  Double_t mean = h_EovEexp_calib->GetMean(), stdev = h_EovEexp_calib->GetStdDev();
  Double_t lower_lim = h_EovEexp_min + maxBin*binW - EovEexp_fit_width*stdev;
  Double_t upper_lim = h_EovEexp_min + maxBin*binW + EovEexp_fit_width*stdev; 
  TF1* fitg = new TF1("fitg","gaus",h_EovEexp_min,h_EovEexp_max);
  fitg->SetRange(lower_lim,upper_lim);
  fitg->SetParameters(norm,mean,stdev);
  fitg->SetLineWidth(2); fitg->SetLineColor(1);
  h_EovEexp_calib->Fit(fitg,"QR"); fitg->GetParameters(param); sigerr = fitg->GetParError(2);
  h_EovEexp_calib->SetLineWidth(2); h_EovEexp_calib->SetLineColor(kGreen+2); h_EovEexp_calib->SetFillColorAlpha(kGreen+2,0.4);
  // adjusting histogram height for the legend to fit properly
  h_EovEexp_calib->GetYaxis()->SetRangeUser(0.,max(norm,norm_bc)*1.2);
  h_EovEexp_calib->Draw(); h_EovEexp->Draw("same");
  // draw the legend
  TLegend *l = new TLegend(0.10,0.78,0.90,0.90);
  l->SetTextFont(42);
  l->AddEntry(h_EovEexp,Form("Before calib., #mu = %.2f, #sigma = (%.3f #pm %.3f) p",param_bc[1],param_bc[2]*100,sigerr_bc*100),"l");
  l->AddEntry(h_EovEexp_calib,Form("After calib., #mu = %.2f, #sigma = (%.3f #pm %.3f) p",param[1],param[2]*100,sigerr*100),"l");
  l->Draw();
  // c1->cd(2); //
  // TPad *pad1 = new TPad("pad1", "Top Pad", 0.0, 0.5, 1.0, 1.0);
  // pad1->Draw();
  // pad1->cd();
  // gPad->SetGridy();
  // gStyle->SetErrorX(0.0001);
  // h2_EovEexp_vs_P->SetStats(0);
  // h2_EovEexp_vs_P->Draw("colz");
  // h2_EovEexp_vs_P_prof->Draw("same");
  // c1->cd(2); //
  // TPad *pad2 = new TPad("pad2", "Bottom Pad", 0.0, 0.0, 1.0, 0.5);
  // pad2->Draw();
  // pad2->cd();
  // gPad->SetGridy();
  // gStyle->SetErrorX(0.0001);
  // h2_EovEexp_vs_P_calib->SetStats(0);
  // h2_EovEexp_vs_P_calib->Draw("colz");
  // h2_EovEexp_vs_P_prof_calib->Draw("same");
  // c1->cd(4); //
  // h2_EovEexp_vs_blk->SetStats(0);
  // h2_EovEexp_vs_blk->Draw("colz text");
  // c1->cd(5); //
  // h2_EovEexp_vs_blk_calib->SetStats(0);
  // h2_EovEexp_vs_blk_calib->Draw("colz text");
  // c1->cd(6); //
  // h2_EovEexp_vs_PSblk_calib->SetStats(0);
  // h2_EovEexp_vs_PSblk_calib->Draw("colz text");
  c1->cd(2);
  gPad->SetGridx();
  h_eSF->SetLineWidth(2); h_eSF->SetLineColor(kRed); h_eSF->SetFillColorAlpha(kRed,0.4);
  h_eSF_calib->SetLineWidth(2); h_eSF_calib->SetLineColor(kGreen+2); h_eSF_calib->SetFillColorAlpha(kGreen+2,0.4);  
  h_eSF_calib->Draw();
  h_eSF->Draw("same");
  c1->cd(3);
  gPad->SetGridx();
  h_eHCAL->SetLineWidth(2); h_eHCAL->SetLineColor(kRed); h_eHCAL->SetFillColorAlpha(kRed,0.4);
  h_eHCAL_calib->SetLineWidth(2); h_eHCAL_calib->SetLineColor(kGreen+2); h_eHCAL_calib->SetFillColorAlpha(kGreen+2,0.4);  
  h_eHCAL_calib->Draw();
  h_eHCAL->Draw("same");
  c1->SaveAs(Form("%s[",outPlot.Data())); c1->SaveAs(Form("%s",outPlot.Data())); c1->Write();
  //**** -- ***//

  // /**** Canvas 2 (Corr. with tr vars.) ****/
  TCanvas *c2 = new TCanvas("c2","E/Eexp vs blk",1500,1200);
  c2->Divide(2,1);
  c2->cd(1); //
  h2_EovEexp_vs_blk->SetStats(0);
  h2_EovEexp_vs_blk->Draw("colz");
  c2->cd(2); //
  h2_EovEexp_vs_blk_calib->SetStats(0);
  h2_EovEexp_vs_blk_calib->Draw("colz");
  c2->SaveAs(Form("%s",outPlot.Data())); c2->Write();
  //**** -- ***//
  
  // /**** Canvas 3 (Corr. with tr vars.) ****/
  // TCanvas *c3 = new TCanvas("c3","E/Eexp vs x, y exp",1500,1200);
  // c3->Divide(2,2);
  // c3->cd(1); //
  // gPad->SetGridy();
  // h2_EovEexp_vs_xHCALexp->SetStats(0);
  // h2_EovEexp_vs_xHCALexp->Draw("colz");
  // c3->cd(2); //
  // gPad->SetGridy();
  // h2_EovEexp_vs_yHCALexp->SetStats(0);
  // h2_EovEexp_vs_yHCALexp->Draw("colz");
  // c3->cd(3); //
  // gPad->SetGridy();
  // h2_EovEexp_vs_xHCALexp_calib->SetStats(0);
  // h2_EovEexp_vs_xHCALexp_calib->Draw("colz");
  // c3->cd(4); //
  // gPad->SetGridy();
  // h2_EovEexp_vs_yHCALexp_calib->SetStats(0);
  // h2_EovEexp_vs_yHCALexp_calib->Draw("colz");
  // c3->SaveAs(Form("%s",outPlot.Data())); c3->Write();
  // //**** -- ***//

  // /**** Canvas 4 (Corr. with tr vars.) ****/
  // TCanvas *c4 = new TCanvas("c4","eHCAL vs x, y exp",1500,1200);
  // c4->Divide(2,2);
  // c4->cd(1); //
  // gPad->SetGridy();
  // h2_eHCAL_vs_xHCALexp->SetStats(0);
  // h2_eHCAL_vs_xHCALexp->Draw("colz");
  // c4->cd(2); //
  // gPad->SetGridy();
  // h2_eHCAL_vs_yHCALexp->SetStats(0);
  // h2_eHCAL_vs_yHCALexp->Draw("colz");
  // c4->cd(3); //
  // gPad->SetGridy();
  // h2_eHCAL_vs_xHCALexp_calib->SetStats(0);
  // h2_eHCAL_vs_xHCALexp_calib->Draw("colz");
  // c4->cd(4); //
  // gPad->SetGridy();
  // h2_eHCAL_vs_yHCALexp_calib->SetStats(0);
  // h2_eHCAL_vs_yHCALexp_calib->Draw("colz");
  // c4->SaveAs(Form("%s",outPlot.Data())); c4->Write();
  // //**** -- ***//

  // /**** Canvas 5 (position resolution) ****/
  // TCanvas *c5 = new TCanvas("c5","pos. res.",1200,1000);
  // c5->Divide(2,2);  gStyle->SetOptFit(1111);
  // c5->cd(1); //
  // TF1* fit_c51 = new TF1("fit_c51","gaus",-0.5,0.5);
  // h_xHCAL_diff->Fit(fit_c51,"QR");
  // h_xHCAL_diff->SetStats(1);
  // c5->cd(2); //
  // TF1* fit_c52 = new TF1("fit_c52","gaus",-0.5,0.5);
  // h_yHCAL_diff->Fit(fit_c52,"QR");
  // h_yHCAL_diff->SetStats(1);
  // c5->cd(3); //
  // TF1* fit_c53 = new TF1("fit_c53","gaus",-0.5,0.5);
  // h_xHCAL_diff_calib->Fit(fit_c53,"QR");
  // h_xHCAL_diff_calib->SetStats(1);
  // c5->cd(4); //
  // TF1* fit_c54 = new TF1("fit_c54","gaus",-0.5,0.5);
  // h_yHCAL_diff_calib->Fit(fit_c54,"QR");
  // h_yHCAL_diff_calib->SetStats(1);
  // c5->SaveAs(Form("%s",outPlot.Data())); c5->Write();
  // //**** -- ***//

  // // /**** Canvas 5 (E/Eexp vs. run number) ****/
  // // TCanvas *c5 = new TCanvas("c5","E/Eexp vs rnum",1200,1000);
  // // c5->Divide(1,2);
  // // // // manipulating urnum vector
  // // // std::size_t nrun = lrnum.size();
  // // // if (nrun!=Nruns)
  // // //   std::cout << "*!*[WARNING] 'Nruns' value in run list doesn't match with total # runs analyzed!\n\n"; 
  // // c5->cd(1); //
  // // gPad->SetGridy();
  // // gStyle->SetErrorX(0.0001); 
  // // Custm2DRnumHisto(h2_EovEexp_vs_rnum,lrnum);
  // // h2_EovEexp_vs_rnum->Draw("colz");
  // // h2_EovEexp_vs_rnum_prof->Draw("same");
  // // c5->cd(2); //
  // // gPad->SetGridy();
  // // gStyle->SetErrorX(0.0001);
  // // Custm2DRnumHisto(h2_EovEexp_vs_rnum_calib,lrnum);
  // // h2_EovEexp_vs_rnum_calib->Draw("colz");
  // // h2_EovEexp_vs_rnum_calib_prof->Draw("same");
  // // c5->SaveAs(Form("%s",outPlot.Data())); c5->Write();
  // // //**** -- ***//

  // /**** Canvas 7 (gain coefficients TOP) ****/
  // TCanvas *c7 = new TCanvas("c7","gain Coeff top",1200,1000);
  // TStyle* oldStyle = (TStyle*)gStyle->Clone("oldStyle");
  // gStyle->SetPaintTextFormat("0.4f");  
  // c7->Divide(1,2);
  // c7->cd(1); Double_t h_max;  
  // h_max = h_coeff_blk_old->GetMaximum();
  // h2_coeff_detView_old_top->GetZaxis()->SetRangeUser(0.002,h_max); h2_coeff_detView_old_top->Draw("text col");
  // c7->cd(2); //
  // h_max = h_coeff_blk->GetMaximum();
  // h2_coeff_detView_top->GetZaxis()->SetRangeUser(0.002,h_max); h2_coeff_detView_top->Draw("text col");
  // //
  // c7->SaveAs(Form("%s",outPlot.Data())); c7->Write();
  // //**** -- ***//

  // /**** Canvas 8 (gain coefficients MIDDLE) ****/
  // TCanvas *c8 = new TCanvas("c8","gain Coeff mid",1200,1000);
  // c8->Divide(1,2);
  // c8->cd(1); 
  // h_max = h_coeff_blk_old->GetMaximum();
  // h2_coeff_detView_old_mid->GetZaxis()->SetRangeUser(0.002,h_max); h2_coeff_detView_old_mid->Draw("text col");
  // c8->cd(2); //
  // h_max = h_coeff_blk->GetMaximum();
  // h2_coeff_detView_mid->GetZaxis()->SetRangeUser(0.002,h_max); h2_coeff_detView_mid->Draw("text col");
  // //
  // c8->SaveAs(Form("%s",outPlot.Data())); c8->Write();
  // //**** -- ***//

  // /**** Canvas 9 (gain coefficients BOTTOM) ****/
  // TCanvas *c9 = new TCanvas("c9","gain Coeff bot",1200,1000);
  // c9->Divide(1,2);
  // c9->cd(1); 
  // h_max = h_coeff_blk_old->GetMaximum();
  // h2_coeff_detView_old_bot->GetZaxis()->SetRangeUser(0.002,h_max); h2_coeff_detView_old_bot->Draw("text col");
  // c9->cd(2); //
  // h_max = h_coeff_blk->GetMaximum();
  // h2_coeff_detView_bot->GetZaxis()->SetRangeUser(0.002,h_max); h2_coeff_detView_bot->Draw("text col");
  // //
  // c9->SaveAs(Form("%s",outPlot.Data())); c9->Write();
  // //
  // gStyle->SetPaintTextFormat(oldStyle->GetPaintTextFormat());
  // delete oldStyle;
  // //**** -- ***//

  /**** Canvas 10 (# events per block) ****/
  TCanvas *c10 = new TCanvas("c10","good ev per blk",1500,1200);
  c10->cd();
  gStyle->SetPaintTextFormat("0.0f");
  h2_nev_per_blk->SetMarkerSize(0.6);
  h2_nev_per_blk->Draw("colz text");
  // c10->Divide(1,2);
  // c10->cd(1); //
  // h2_nev_per_blk_bot->Draw("colz text");
  // c10->cd(2); //
  // h2_nev_per_blk_top->Draw("colz text");
  c10->SaveAs(Form("%s",outPlot.Data())); c10->Write();
  //**** -- ***//

  // // /**** Canvas 9 (cluster size) ****/
  // // TCanvas *c9 = new TCanvas("c9","cl. size vs rnum",1200,1000);
  // // c9->Divide(1,4);
  // // c9->cd(1); //
  // // Custm2DRnumHisto(h2_PSclsize_vs_rnum,lrnum);
  // // h2_PSclsize_vs_rnum->Draw("colz");
  // // h2_PSclsize_vs_rnum_prof->Draw("same");
  // // c9->cd(2); //
  // // Custm2DRnumHisto(h2_PSclmult_vs_rnum,lrnum);
  // // h2_PSclmult_vs_rnum->Draw("colz");
  // // h2_PSclmult_vs_rnum_prof->Draw("same");
  // // c9->cd(3); //
  // // Custm2DRnumHisto(h2clsize_vs_rnum,lrnum); 
  // // h2clsize_vs_rnum->Draw("colz");
  // // h2clsize_vs_rnum_prof->Draw("same");
  // // c9->cd(4); //
  // // Custm2DRnumHisto(h2clmult_vs_rnum,lrnum);
  // // h2clmult_vs_rnum->Draw("colz");
  // // h2clmult_vs_rnum_prof->Draw("same");
  // // c9->SaveAs(Form("%s",outPlot.Data())); c9->Write();
  // // //**** -- ***//

  /**** Summary Canvas ****/
  TCanvas *cSummary = new TCanvas("cSummary","Summary");
  cSummary->cd();
  TPaveText *pt = new TPaveText(.05,.1,.95,.8);
  pt->AddText(Form(" Date of creation: %s",GetDate().c_str()));
  pt->AddText(Form("Configfile: ECAL_replay/scripts/cfg/%s.cfg",cfgfilebase.Data()));
  pt->AddText(Form(" Total # events analyzed: %lld ",Nevents));
  pt->AddText(Form(" E/Eexp (before calib.) | #mu = %.2f, #sigma = (%.3f #pm %.3f) p",param_bc[1],param_bc[2]*100,sigerr_bc*100));
  pt->AddText(Form(" E/Eexp (after calib.)    | #mu = %.2f, #sigma = (%.3f #pm %.3f) p",param[1],param[2]*100,sigerr*100));
  pt->AddText(" Global cuts: ");
  std::string tmpstr = "";
  for (std::size_t i=0; i<gCutList.size(); i++) {
    if (i>0 && i%3==0) {pt->AddText(Form(" %s",tmpstr.c_str())); tmpstr="";}
    tmpstr += gCutList[i] + ", "; 
  }
  if (!tmpstr.empty()) pt->AddText(Form(" %s",tmpstr.c_str()));
  if (cut_on_eHCAL) pt->AddText(Form(" HCAL cluster energy < %.1f GeV",eHCAL_cut_limit));
  if (cut_on_EovEexp) pt->AddText(Form(" |E/Eexp - 1| < %.1f",EovEexp_cut_limit));
  pt->AddText(Form(" # events passed global cuts: %lld", Ngoodevs));
  pt->AddText(" Other cuts: ");
  pt->AddText(Form(" Minimum # events per block: %d | Cluster hit threshold: %.2f GeV",Nmin,hit_threshold));
  pt->AddText(Form(" Cluster tmax cut: %.1f ns | Cluster energy fraction cut: %.1f GeV",tmax_cut,engFrac_cut));
  pt->AddText(" ECAL position offsets: ");
  pt->AddText(Form(" v_offset: %.3f (m), h_offset: %.3f (m), z_offset: %.3f (m)",ECAL_voff,ECAL_hoff,ECAL_zoff));
  sw->Stop(); sw2->Stop();
  pt->AddText(Form("Macro processing time: CPU %.1fs | Real %.1fs",sw->CpuTime(),sw->RealTime()));
  TText *t1 = pt->GetLineWith("Configfile"); t1->SetTextColor(kRed+2);
  TText *t2 = pt->GetLineWith(" E/Eexp (be"); t2->SetTextColor(kRed);
  TText *t3 = pt->GetLineWith(" E/Eexp (af"); t3->SetTextColor(kGreen+2);
  TText *t4 = pt->GetLineWith(" Global"); t4->SetTextColor(kBlue);
  TText *t5 = pt->GetLineWith(" Other"); t5->SetTextColor(kBlue);
  // TText *t6 = pt->GetLineWith(" Various"); t6->SetTextColor(kBlue);
  TText *t7 = pt->GetLineWith("Macro"); t7->SetTextColor(kGreen+3);
  pt->Draw();
  cSummary->SaveAs(Form("%s",outPlot.Data())); cSummary->SaveAs(Form("%s]",outPlot.Data())); cSummary->Write();  
  //**** -- ***//

  std::cout << "List of output files:" << "\n";
  std::cout << " --------- " << "\n";
  std::cout << " 1. Summary plots : "        << outPlot << "\n";
  std::cout << " 2. Resulting histograms : " << outFile << "\n";
  std::cout << " 3. Gain ratios (new/old) : " << outGainR << "\n";
  std::cout << " 5. New ADC gain coeffs. (GeV/pC) : " << outGain << "\n";
  std::cout << " --------- " << "\n";
  
  std::cout << "CPU time = " << sw->CpuTime() << "s. Real time = " << sw->RealTime() << "s.\n\n";

  ///////////////////////////////////////////////////
  // Write individual memories to file explicitely //
  // to be able to read them using uproot          //
  ///////////////////////////////////////////////////  
  Tout->Write("", TObject::kOverwrite);
//   // diagnostic histos
  h_eSF->Write(); h_eSF_calib->Write();      
  h_eHCAL->Write(); h_eHCAL_calib->Write();    
  h_EovEexp->Write(); h_EovEexp_calib->Write();
//   h_xHCAL_diff->Write(); h_xHCAL_diff_calib->Write();
//   h_yHCAL_diff->Write(); h_yHCAL_diff_calib->Write();      
//   // correlations
//   h2_EovEexp_vs_P->Write(); h2_EovEexp_vs_P_calib->Write();
//   h2_EovEexp_vs_P_prof->Write(); h2_EovEexp_vs_P_prof_calib->Write();
//   h2_EovEexp_vs_blk->Write(); h2_EovEexp_vs_blk_calib->Write();
//   h2_EovEexp_vs_xHCALexp->Write(); h2_EovEexp_vs_xHCALexp_calib->Write();
//   h2_EovEexp_vs_yHCALexp->Write(); h2_EovEexp_vs_yHCALexp_calib->Write();
//   h2_eHCAL_vs_xHCALexp->Write(); h2_eHCAL_vs_xHCALexp_calib->Write();
//   h2_eHCAL_vs_yHCALexp->Write(); h2_eHCAL_vs_yHCALexp_calib->Write();
//   // gains
//   h_nevent_blk->Write(); h_coeff_Ratio_blk->Write();
//   h_coeff_blk->Write(); h_coeff_blk_old->Write();
//   h2_coeff_detView_bot->Write(); h2_coeff_detView_old_bot->Write();  
//   h2_coeff_detView_mid->Write(); h2_coeff_detView_old_mid->Write();  
//   h2_coeff_detView_top->Write(); h2_coeff_detView_old_top->Write();  
//   h2_nev_per_blk->Write();
}

// **** ========== Useful functions ========== ****  
// returns today's date
std::string GetDate(){
  time_t now = time(0);
  tm ltm = *localtime(&now);
  std::string yyyy = to_string(1900 + ltm.tm_year);
  std::string mm = to_string(1 + ltm.tm_mon);
  std::string dd = to_string(ltm.tm_mday);
  std::string date = mm + '/' + dd + '/' + yyyy;
  return date;
}

// splits a string by a delimiter (doesn't include empty sub-strings)
std::vector<std::string> SplitString(char const delim, std::string const myStr) {
  std::stringstream ss(myStr);
  std::vector<std::string> out;
  while (ss.good()) {
    std::string substr;
    std::getline(ss, substr, delim);
    if (!substr.empty()) out.push_back(substr);
  }
  if (out.empty()) std::cerr << "WARNING! No substrings found!\n";
  return out;
}

// returns output file base from configfilename
TString GetOutFileBase(TString configfilename) {
  std::vector<std::string> result;
  result = SplitString('/',configfilename.Data());
  TString temp = result[result.size() - 1];
  return temp.ReplaceAll(".cfg", "");
}

// customizes profile histograms
void CustmProfHisto(TH1* hprof) {
  hprof->SetStats(0);
  hprof->SetMarkerStyle(20);
  hprof->SetMarkerColor(2);
}

// Forms the matrices
void BuildMatrix(const double *A, TVectorD &B, TMatrixD &M, Double_t p_e) {
#pragma omp parallel for schedule(static)
  for (Int_t icol = 0; icol < kNblks; icol++) {
    Double_t a_col = A[icol];  // Cache A[icol]
    B(icol) += a_col;          // Update B using TVectorD
    for (Int_t irow = 0; irow < kNblks; irow++) {
      M(icol, irow) += a_col * A[irow] / p_e;
    }
  }
}

// Reads from a single column text file (no header)
std::vector<int> ReadNumFromTextFileToV(const std::string &filename) {
  std::vector<int> numVec;
  std::ifstream file(filename);
  int num;
  // Check if the file was successfully opened
  if (!file) {
    throw std::runtime_error("Error: Unable to open file " + filename);
  }
  while (file >> num) {
    numVec.push_back(num);
  }
  return numVec;
}

// Reads from a text file
std::unordered_set<int> GetEdgeBlkID(LookUpTableReader& ecalmap, bool isdebug) {
  std::unordered_set<int> edgeblkid;
  for (int iblk=0; iblk<kNblks; iblk++) {
    bool isonedge = int(ecalmap.GetValueByKey(iblk,5));
    if (isonedge) {
      edgeblkid.insert(iblk);
      if (isdebug) {
	std::cout << Form("Block ID on Edge: %d\n",iblk);
      }  
    }
  }
  return edgeblkid;
}

// Extracts # of Columns per Row from the HCAL map
std::unordered_map<int,int> GetNColPerRow(LookUpTableReader& ecalmap, bool isdebug) {
  std::unordered_map<int,int> ncolprow;
  for (int iblk=0; iblk<kNblks; iblk++) {
    int row = int(ecalmap.GetValueByKey(iblk,0));
    int ncol = int(ecalmap.GetValueByKey(iblk,4));    
    ncolprow[row] = ncol;
  }
  if (isdebug) {
    std::vector<int> keys;
    for (const auto& pair : ncolprow) {
        keys.push_back(pair.first);
    }
    std::sort(keys.begin(), keys.end());
    for (int key : keys) {
      std::cout << Form("Row %d has %d Columns\n",key,ncolprow[key]);
    }    
  }  
  return ncolprow;  
}

// Read the x,y positions of the known BAD channels
std::vector<std::pair<double, double>> GetXYPosBADChan(LookUpTableReader& ecalmap,
						       std::vector<int>& badchanlist, bool isdebug) {
  std::vector<std::pair<double, double>> positions;
  for (int chan : badchanlist) {
    double xpos = ecalmap.GetValueByKey(chan,2);
    double ypos = ecalmap.GetValueByKey(chan,3);
    positions.emplace_back(xpos, ypos);

    if (isdebug) {
      std::cout << Form("BAD Chan ID %d, xpos %f, ypos %f\n",chan,xpos,ypos);
    }    
  }
  return positions;
}

// Checks if the seed is near a known bad channel
bool ISNearBADChan(std::vector<std::pair<double, double>> const &xyBADchan, double xSEED, double ySEED, double r2_max, int debug) {
  bool nearBADchan = false;
  double dx, dy;
  for (const auto& [refX, refY] : xyBADchan) {
    dx = xSEED - refX;
    dy = ySEED - refY;
    if ((dx * dx + dy * dy) < r2_max) {
      nearBADchan = true;
      break;
    }
  }
  if (nearBADchan && debug>1) {
    std::cout << Form("xSEED - xBADch=%.6f, ySEED - yBADch=%.6f\n",dx,dy);
  }
  return nearBADchan;
}
