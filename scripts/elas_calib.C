#include <iostream>
#include <fstream>
#include <vector>

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
const int kNblks = 1656;
const int kNrows = 69;
const int kNcols = 27; // varies by row. 27 is the maxc.

void elas_calib()
{
  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  // Define a clock to check macro processing time
  TStopwatch *sw = new TStopwatch();
  TStopwatch *sw2 = new TStopwatch();
  sw->Start(); sw2->Start();  

  TChain *C = new TChain("T");
  //creating base for outfile names
  TString cfgfilebase = "elas_calib"; //= GetOutFileBase(configfilename);

  C->Add("/lustre24/expphy/volatile/halla/sbs/pdbforce/gep/mc/elas/replayed_GEP1_job*.root");

  // Double_t E_beam=0.;
  // Double_t p_p=0., E_p=0., th_p=0., ph_p=0.; 
  // Double_t p_e=0., th_e=0., ph_e=0.;
  // Double_t p_e_th; // e- momentum calculated using angles

  TCut globalcut = "sbs.tr.p[0]<6"; TString gcutstr;
  // while (currentline.ReadLine(configfile) && !currentline.BeginsWith("endcut")) {
  //   if (!currentline.BeginsWith("#")) {
  //     globalcut += currentline;
  //     gcutstr += currentline;
  //   }    
  // }
  // std::vector<std::string> gCutList = SplitString('&', gcutstr.Data());
  TTreeFormula *GlobalCut = new TTreeFormula("GlobalCut", globalcut, C);  

  // Turning on relevant tree branches
  int maxNtr = 200; //max # of tracks expected per event
  C->SetBranchStatus("*", 0);
  // bb.ps branches
  C->SetBranchStatus("earm.ecal.*", 1);
  Double_t nclusECAL;              C->SetBranchAddress("earm.ecal.nclus", &nclusECAL);
  Double_t idblkECAL;              C->SetBranchAddress("earm.ecal.idblk", &idblkECAL);
  Double_t rowblkECAL;             C->SetBranchAddress("earm.ecal.rowblk", &rowblkECAL);
  Double_t colblkECAL;             C->SetBranchAddress("earm.ecal.colblk", &colblkECAL);
  Double_t nblkECAL;               C->SetBranchAddress("earm.ecal.nblk", &nblkECAL);
  Double_t atimeECAL;              C->SetBranchAddress("earm.ecal.atimeblk", &atimeECAL);
  Double_t eECAL;                  C->SetBranchAddress("earm.ecal.e", &eECAL);
  Double_t xECAL;                  C->SetBranchAddress("earm.ecal.x", &xECAL);
  Double_t yECAL;                  C->SetBranchAddress("earm.ecal.y", &yECAL);
  Double_t idclblkECAL[maxNtr];    C->SetBranchAddress("earm.ecal.clus_blk.id", &idclblkECAL);
  Double_t eclblkECAL[maxNtr];     C->SetBranchAddress("earm.ecal.clus_blk.e", &eclblkECAL);
  Double_t xclblkECAL[maxNtr];     C->SetBranchAddress("earm.ecal.clus_blk.x", &xclblkECAL);
  Double_t yclblkECAL[maxNtr];     C->SetBranchAddress("earm.ecal.clus_blk.y", &yclblkECAL);
  Double_t rowclblkECAL[maxNtr];   C->SetBranchAddress("earm.ecal.clus_blk.row", &rowclblkECAL);
  Double_t colclblkECAL[maxNtr];   C->SetBranchAddress("earm.ecal.clus_blk.col", &colclblkECAL);
  Double_t atimeclblkECAL[maxNtr]; C->SetBranchAddress("earm.ecal.clus_blk.atime", &atimeclblkECAL);
  Double_t againblkECAL;           C->SetBranchAddress("earm.ecal.againblk", &againblkECAL);  
  // sbs.hcal branches
  Double_t eHCAL;                C->SetBranchStatus("sbs.hcal.e",1); C->SetBranchAddress("sbs.hcal.e", &eHCAL);
  Double_t xHCAL;                C->SetBranchStatus("sbs.hcal.x",1); C->SetBranchAddress("sbs.hcal.x", &xHCAL);
  Double_t yHCAL;                C->SetBranchStatus("sbs.hcal.y",1); C->SetBranchAddress("sbs.hcal.y", &yHCAL); 
  Double_t atimeHCAL;            C->SetBranchStatus("sbs.hcal.atimeblk",1);   C->SetBranchAddress("sbs.hcal.atimeblk", &atimeHCAL); 
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
  // Event info
  C->SetMakeClass(1);
  // C->SetBranchStatus("fEvtHdr.*", 1);
  // UInt_t rnum;                   C->SetBranchAddress("fEvtHdr.fRun", &rnum);
  // UInt_t trigbits;               C->SetBranchAddress("fEvtHdr.fTrigBits", &trigbits);
  // ULong64_t gevnum;              C->SetBranchAddress("fEvtHdr.fEvtNum", &gevnum);
  // turning on additional branches for the global cut
  // C->SetBranchStatus("sbs.gem.track.nhits", 1);
  // C->SetBranchStatus("sbs.gem.track.ngoodhits", 1);
  // C->SetBranchStatus("sbs.gem.track.chi2ndf", 1);

  // Creating output ROOT file to contain histograms
  TString outFile, outPlot;
  outFile = "test_elas_calib.root";
  // char const * debug = isdebug ? "_test" : "";
  // char const * elcut = elastic_cut ? "_elcut" : "";
  // outFile = Form("%s/hist/%s_prepass%d_bbcal_eng_calib%s%s.root",macros_dir.Data(),cfgfilebase.Data(),ppass,elcut,debug);
  // outPlot = Form("%s/plots/%s_prepass%d_bbcal_eng_calib%s%s.pdf",macros_dir.Data(),cfgfilebase.Data(),ppass,elcut,debug);

  //std::unique_ptr<TFile> fout( TFile::Open(outFile, "RECREATE") );
  TFile *fout = new TFile(outFile, "RECREATE");
  fout->cd();

  
  // // Clear arrays
  // memset(nevents_per_cell, 0, ncell*sizeof(int));
  // memset(badCells, 0, ncell*sizeof(bool));

  // defining output ROOT tree (Set max size to 4GB)
  TTree *Tout = new TTree("Tout", cfgfilebase.Data()); 
  Tout->SetMaxTreeSize(4000000000LL);
  //
  Double_t E_beam;     Tout->Branch("E_beam", &E_beam, "E_beam/D");
  Double_t p_p;        Tout->Branch("p_p", &p_p, "p_p/D");
  Double_t p_px;       Tout->Branch("p_px", &p_px, "p_px/D");
  Double_t p_py;       Tout->Branch("p_py", &p_py, "p_py/D");
  Double_t p_pz;       Tout->Branch("p_pz", &p_pz, "p_pz/D");  
  Double_t E_p;        Tout->Branch("E_p", &E_p, "E_p/D");
  Double_t th_p;       Tout->Branch("th_p", &th_p, "th_p/D");
  Double_t ph_p;       Tout->Branch("ph_p", &ph_p, "ph_p/D");
  Double_t p_e;        Tout->Branch("p_e", &p_e, "p_e/D");
  Double_t th_e;       Tout->Branch("th_e", &th_e, "th_e/D");
  Double_t ph_e;       Tout->Branch("ph_e", &ph_e, "ph_e/D");
  Double_t p_e_th;     Tout->Branch("p_e_th", &p_e_th, "p_e_th/D");
  //
  Double_t nu;         Tout->Branch("nu", &nu, "nu/D"); 
  Double_t W2;         Tout->Branch("W2", &W2, "W2/D"); 
  Double_t Q2;         Tout->Branch("Q2", &Q2, "Q2/D");
  Double_t dx;         Tout->Branch("dx", &dx, "dx/D"); 
  Double_t dy;         Tout->Branch("dy", &dy, "dy/D");
  //
  Double_t nu_4vec;    Tout->Branch("nu_4vec", &nu_4vec, "nu_4vec/D"); 
  Double_t W2_4vec;    Tout->Branch("W2_4vec", &W2_4vec, "W2_4vec/D"); 
  Double_t Q2_4vec;    Tout->Branch("Q2_4vec", &Q2_4vec, "Q2_4vec/D");
  Double_t dx_4vec;    Tout->Branch("dx_4vec", &dx_4vec, "dx_4vec/D"); 
  Double_t dy_4vec;    Tout->Branch("dy_4vec", &dy_4vec, "dy_4vec/D");  
  

  Double_t th_bb = 29.46*TMath::DegToRad();
  Double_t ECAL_dist = 8.0;
  Double_t ECAL_zoff = 0., ECAL_voff = 0., ECAL_hoff = 0.;
  
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

      // // track change of runnum
      // if (nevent == 1 || rnum != runnum) {
      // 	runnum = rnum; itrrun++;
      // 	lrnum.push_back(to_string(rnum));
      // }
    } 
    bool passedgCut = GlobalCut->EvalInstance(0) != 0;   
    if (passedgCut) {        

      // CHANGE THIS
      E_beam = 6.476;

      TVector3 vertex(0,0,trVz[0]);
      
      // scattered proton kinematics (reconstructed)
      p_p = trP[0];
      p_px = trPx[0];
      p_py = trPy[0];
      p_pz = trPz[0];    
      E_p = sqrt(p_p*p_p + Mp*Mp);
      th_p = TMath::ACos(trPz[0]/trP[0]);
      ph_p = atan2(trPy[0],-trPx[0]); 

      // scattered electron kinematics (from p kinematics)
      p_e = E_beam + Mp - E_p;
      th_e = 2.0 * TMath::ATan((Mp/(Mp+E_beam))*(1./TMath::Tan(th_p)));
      ph_e = - ph_p;
      p_e_th = E_beam/(1. + (E_beam/Mp)*(1.0 - cos(th_e)));

      // kinematic variables
      nu = E_beam - p_e_th; 
      Q2 = 2.*E_beam*p_e_th*(1.-cos(th_e)); 
      W2 = Mp*Mp + 2.*Mp*nu - Q2;
      TVector3 Pef_3vec(p_e_th*sin(th_e)*cos(ph_e), p_e_th*sin(th_e)*sin(ph_e), p_e_th*cos(th_e));    
      TVector3 Pefhat = Pef_3vec.Unit();
      // calculating expected hit positions on ECAL (angles only method)
      Double_t sintersect = (ECAL_origin - vertex).Dot(ECAL_zaxis) / (Pefhat.Dot(ECAL_zaxis));
      TVector3 ECAL_intersect = vertex + sintersect*Pefhat; 
      Double_t xECAL_exp = (ECAL_intersect - ECAL_origin).Dot(ECAL_xaxis);
      Double_t yECAL_exp = (ECAL_intersect - ECAL_origin).Dot(ECAL_yaxis);
      dx = xECAL - xECAL_exp;
      dy = yECAL - yECAL_exp;    

      // elastic calculations (Using 4-vector method)
      // Relevant 4-vectors
      /* Reaction    : e + p -> e' + p'
	 Conservation: Pei + Ppi = Pef + Ppf */
      TLorentzVector Pei(0,0,E_beam,E_beam);  // incoming e- 4-vector
      TLorentzVector Ppi(0,0,0,Mp);           // target nucleon 4-vector
      TLorentzVector Ppf(p_px,                // Struck nucleon 4-vector
			 p_py,
			 p_pz,
			 E_p);
      TLorentzVector Pef = Pei + Ppi - Ppf;   // scattered e- 4-vector
      //
      TLorentzVector q = Ppf - Ppi;                // 4-momentum of virtual photon
      nu_4vec = q.E();
      Q2_4vec = -q.M2();
      W2_4vec = Ppf.M2();
      TVector3 Pefhat_4vec = Pef.Vect().Unit();
      // calculating expected hit positions on ECAL (4-vect method)
      Double_t sintersect_4vec = (ECAL_origin - vertex).Dot(ECAL_zaxis) / (Pefhat_4vec.Dot(ECAL_zaxis));
      TVector3 ECAL_intersect_4vec = vertex + sintersect_4vec*Pefhat_4vec; 
      Double_t xECAL_exp_4vec = (ECAL_intersect_4vec - ECAL_origin).Dot(ECAL_xaxis);
      Double_t yECAL_exp_4vec = (ECAL_intersect_4vec - ECAL_origin).Dot(ECAL_yaxis);
      dx_4vec = xECAL - xECAL_exp_4vec;
      dy_4vec = yECAL - yECAL_exp_4vec;    

    
      Tout->Fill();
    }
  }

  Tout->Write("", TObject::kOverwrite);
}
