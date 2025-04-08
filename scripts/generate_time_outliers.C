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

//update for new cuts 1-17-25, 69, 39, 1
void generate_time_outliers(int runis, int mid_multiplicity = 1, int col_start = 1, int chain_amount = 1){

  int row_start = 68;
  int row_end = 68;

  //string runis;
  int entry_pass_count = 0;

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

  double new_col_xpos_list[27];
  double block_width = 4.292;//cm
  double block_width_SM_boundary = 4.455;//cm //this is the distance between two block who are in separate SM's according to JT file measurement extraction

  int block_mult = 0;
  int block_mult_SM_boundary = 0;
  
  for(int i = 0; i < 28; i++){
    
    new_col_xpos_list[i] = xpos_arr[0] + (block_mult*block_width) + (block_mult_SM_boundary*block_width_SM_boundary);

    if((i+1)%3 != 0){
      block_mult++;
    }else{
      block_mult_SM_boundary++;
    } 
    
  }
  
  int new_col[cells];

  int col_index = 0;

  for(int i = 0; i < cells; i++){
    
    for(int j = 0; j < 28; j++){
      if(xpos_arr[i] > new_col_xpos_list[j] - (block_width/2.0) && xpos_arr[i] < new_col_xpos_list[j] + (block_width/2.0)){
	new_col[i] = j; //vertically alligned column indices, middle section of ecal which has first col starting roughly 1 block width after the first col of the very bottom row, now has col index starting at 1 instead of 0, since they line of up with col index 1 of bottom most row instead, kind of a proxy for xpos, but allows for checking repeated hits along this rough via multiplicity check of unordered map for new_col to repeated hits at that new_col index
      }
    }
    
  }

  //for(int er = 0; er < cells; er++){
  //cout << " New col " << new_col[er] << " pmt index is: " << er << endl;
  //}

  bool column_cut = true;
  bool adj_check = false;;
  
  //cout << " adj_check? 1 or 0" <<endl;
  //cin >> adj_check;
  Int_t ndata_is; 
  Double_t amp_mv[2000];

  Int_t ndata_is_time, ndata_is_int, ndata_is_row, ndata_is_col, ndata_is_ped; 
  Double_t time[2000] = {0.};
  Double_t intg[2000] = {0.};
  Double_t row[1000] = {0.};
  Double_t col[1000] = {0.};
  Double_t ped[2000] = {0.};
  Double_t fa_rc[27][69] = {};
  Double_t ft_rc[27][69] = {};
  Double_t fi_rc[27][69] = {};
  Double_t fp_rc[27][69] = {};
  Double_t maxa[2000] = {0.};
  Double_t maxt[2000] = {0.};
  Double_t maxi[2000] = {0.};
  Double_t maxp[2000] = {0.};
  Double_t fita[2000] = {0.};
  Double_t fitt[2000] = {0.};
  Double_t fiti[2000] = {0.};
  Double_t fitp[2000] = {0.};
  Double_t stda[2000] = {0.};
  Double_t stdt[2000] = {0.};
  Double_t stdi[2000] = {0.};
  Double_t stdp[2000] = {0.};
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

  Int_t ndata_pmt; 
  Double_t pmt[2000] = {0.};

  Double_t x_pos[2000] = {1.34};
  Int_t ndata_x_pos;
  
  TH1F* h_pmt_is[2000];
  TH1F* time_pmt_is[2000];
  TH1F* int_pmt_is[2000];
  TH1F* ped_pmt_is[2000];

  

  TH1F* all_pmt_intg = new TH1F("all_pmt_intg", "", 1656, 0, 1655);
  TH1F* all_pmt_amp_mv = new TH1F("all_pmt_amp_mv", "", 1656, 0, 1655);

  TH2F* col_vs_x = new TH2F("col_vs_x", "", 28, 0, 28, 200, -150.0, 150.0);

  TH2F* time_peak_vs_ch = new TH2F("time_peak_vs_ch_mult5", "", 1656, 0, 1656, 200, 0, 350);

  TH2F* int_vs_int = new TH2F("int_vs_int", "", 100, 0, 150, 100, 0, 150);//simona wants integral vs integral for two hits that occur in the same entry, what about for multiplicity higher than two?

  //simona wants a histogram showing the number of hits per entry, one histo for each pmt, this  would just be the frequency distribution of different numbers of hits per entry. So x axis is number of hits per entry, and y axis is the frequency over all entries. Want to say, for this pmt and the entries in which it was hit, how many many hits were involved
  TH1F* hit_mult_per_entry[1656];

  TH1F* xpos_histo[1656];
  TH1F* col_histo[1656];

  //another plot is the integral as a function of number of hits per entry, this would also be a plot for every pmt
  TH2F* avg_intg_vs_hit_mult_per_entry[1656];

  for(Int_t k=0; k < 1656; k++){

    h_pmt_is[k] = new TH1F("","",90,0,210);
    time_pmt_is[k] = new TH1F("","",200,0,400);
    int_pmt_is[k] = new TH1F("","",150,0,150);
    ped_pmt_is[k]= new TH1F("","",150,0,150);

    hit_mult_per_entry[k] = new TH1F("","",20,0,20);
    avg_intg_vs_hit_mult_per_entry[k] = new TH2F("","",20,0,20, 50,0,100);

    xpos_histo[k] = new TH1F("","", 62, -65.0, 65.0);
    col_histo[k] = new TH1F("","", 28, -1, 27);
  }

  TH1F* h_rc[27][69];
  TH1F* time_rc[27][69];
  TH1F* int_rc[27][69];
  TH1F* ped_rc[27][69];

  for(Int_t r=0; r < 69; r++){
    for(Int_t c=0; c<27; c++){
    h_rc[c][r] = new TH1F("","",500,0,500);
    time_rc[c][r] = new TH1F("","",360,0,360);
    int_rc[c][r] = new TH1F("","",600,0,600);
    ped_rc[c][r]= new TH1F("","",150,0,150);
    }
  }

  TChain chain("T");
  









  /*
  for(int i = runis; i< (runis+chain_amount); i++){
    chain.Add(Form("/adaqfs/home/a-onl/sbs/Rootfiles/ecal_%i_-1.root", i));
  }
  */

  chain.Add(Form("/adaqfs/home/a-onl/sbs/Rootfiles/ecal_%i_-1.root", runis));
  //chain.Add(Form("/adaqfs/home/a-onl/sbs/Rootfiles/ecal_%i_-1_1.root", runis));
  //chain.Add(Form("/adaqfs/home/a-onl/sbs/Rootfiles/ecal_%i_-1_2.root", runis));
  
  
  vector<int> pmt_list;

  int time_mean[1656];
  int time_bad_list[1656] = {0};
  int time_entries_list[1656] = {0};
  int time_outlier_list[1656] = {0};

 chain.SetBranchStatus("*",0);
 chain.SetBranchStatus("*earm.ecal.*",1);
 chain.SetBranchAddress("Ndata.earm.ecal.a_amp_p",&ndata_is);
 chain.SetBranchAddress("earm.ecal.a_amp_p",amp_mv);
 
 chain.SetBranchAddress("Ndata.earm.ecal.a_time",&ndata_is_time);
 chain.SetBranchAddress("earm.ecal.a_time",time);

 chain.SetBranchAddress("Ndata.earm.ecal.a_p",&ndata_is_int);// was p
 chain.SetBranchAddress("earm.ecal.a_p",intg);

 chain.SetBranchAddress("Ndata.earm.ecal.ped",&ndata_is_ped);
 chain.SetBranchAddress("earm.ecal.ped",ped);

 chain.SetBranchAddress("Ndata.earm.ecal.adcrow",&ndata_is_row);
 chain.SetBranchAddress("earm.ecal.adcrow",row);

 chain.SetBranchAddress("Ndata.earm.ecal.adccol",&ndata_is_col);
 chain.SetBranchAddress("earm.ecal.adccol",col);
 
 chain.SetBranchAddress("Ndata.earm.ecal.adcelemID",&ndata_pmt);
 chain.SetBranchAddress("earm.ecal.adcelemID",pmt);//0-1655

 chain.SetBranchAddress("Ndata.earm.ecal.y",&ndata_x_pos);
 //chain.SetBranchAddress("earm.ecal.y",x_pos);//0-1655

 int block_width_check = 1;
 double width_check = 4.3*block_width_check;

 vector<vector<double>> pmt_intg_list(1656);
 vector<vector<double>> pmt_amp_mv_list(1656);
 
  Int_t nentries = chain.GetEntries();


  int entry_passed_list[10000] = {0};
  int mult_count = 0;
  
  for(Int_t i=0; i<nentries; i++){
  //for(Int_t i=0; i<200; i++){

    chain.GetEntry(i);
    int mid = 0;
    int above = 0;
    int below = 0;

    bool ndata_check = false;













    //in the entry loop

    if(column_cut == true){//do adjacent hits cut, where we record an entry only if there were adjacent cols fired(or basically pmts fired within a desired block width
      vector<int> all_hit_new_cols;//just list of all cols that were hit using the new_col index
      vector<double> all_hit_x;//just list of all cols that were hit using the new_col index

      all_hit_new_cols.clear();
      all_hit_x.clear();
      
      
      for(Int_t j=0; j<ndata_pmt; j++){

	if(pmt[j] == 440){
	  continue;
	}

	if(row[j] >= 0 && row[j] <= 69){//look for hits in any row of ecal, 57-69 is top section of ecal
	  int pmtindx = pmt[j];
	  all_hit_x.push_back(xpos_arr[pmtindx]);
	}
      }
      
      unordered_map<int, int> same_col_hit_new_col_frequency;//map of cols and the number of hits minus 1, in that col
      unordered_map<double, int> same_x_hit_x_frequency;//map of cols and the number of hits minus 1, in that col
      
      for(int j = 0; j < all_hit_x.size(); j++){
	
	bool j_is_same_x_as_any_prev_hit_x = false;

	for(int l = 0; l < j; l++){
	  //if(all_hit_x[j] == all_hit_x[l]){
	  if(all_hit_x[j] >= all_hit_x[l] - 2.0 && all_hit_x[j] <= all_hit_x[l] + 2.0){	    
	    j_is_same_x_as_any_prev_hit_x = true; //dont want to double count the number of hits in a col, so when looping over every combination of two hits, every hit up to j has been checked already for repeated hits in that column, so if another hit later in the entry has the same col number as any of the already checked hits, then we dont want to check its combinations again, thus we count a proper multiplicity minus 1
	    //cout << "REPEAT PREVIOUS HIT X, current pmt " << all_hit_x[j] << " and previous one that is repeated " << all_hit_x[l] << " j and l are: " << j << " " << l << endl;
	  }
	}

	if(j_is_same_x_as_any_prev_hit_x == false){
	  for(int k = (j+1); k < all_hit_x.size(); k++){//want to compare all combinations of two hits, making k = j+1 ensures that no combination will occur twice, and dont want repeats so want to make sure we skip j if it is in the same col as any previous repeated col
	    if(all_hit_x[k] >= all_hit_x[j] - 2.0 && all_hit_x[k] <= all_hit_x[j] + 2.0){	    
	      same_x_hit_x_frequency[all_hit_x[j]]++;

	      //cout << "Same_x_hit_x_freq = " << same_x_hit_x_frequency[all_hit_x[j]] << " all_hit_x " << all_hit_x[j] << " j is : " << j << " and k(check all other pmts) is " << k << endl; 
	     
	    }
	  }
	}
	
      }
      //back in the entry loop, no longer in j loop
      
      vector<int> same_col_repeated_new_col;

      
      
      for(auto& entry : same_col_hit_new_col_frequency){
	if(entry.second >= (mid_multiplicity - 1)){//the frequency counter only iterates if there after finding two hits, so if there are only two hits in a col, then the entry.second will be 1
	  same_col_repeated_new_col.push_back(entry.first);	  
	}
      }

      vector<int> same_x_repeated_x;
      vector<int> same_x_repeated_x_freq;

      same_x_repeated_x.clear();
      same_x_repeated_x_freq.clear();
      
      for(auto& entry : same_x_hit_x_frequency){
	if(entry.second >= (mid_multiplicity - 1)){//the frequency counter only iterates if there after finding two hits, so if there are only two hits in a col, then the entry.second will be 1
	  same_x_repeated_x.push_back(entry.first);
	  same_x_repeated_x_freq.push_back(entry.second + 1);
	  //cout << " entry.first " << entry.first << " entry.second " << entry.second << endl;
	}
      }

      int ndata_passed_count = 0;
      
      if(same_x_repeated_x.size() > 0){

	vector<int> col_list;
	vector<int> row_list;
	vector<double> xpos_list;

	for(int col_itr = 0; col_itr < same_x_repeated_x.size(); col_itr++){

	  ndata_passed_count = 0;

	
	  for(Int_t j=0; j<ndata_pmt; j++){
	  
	    Int_t pmtindx1 = pmt[j];
	 	    
	    //if(x_pos[j] == same_x_repeated_x[col_itr]){
	    if(xpos_arr[pmtindx1] >= same_x_repeated_x[col_itr] - 2.0 && xpos_arr[pmtindx1] <= same_x_repeated_x[col_itr] + 2.0){//??????????????
	    if(pmtindx1 != 440){
	      ndata_passed_count++;
	    }
	      //if(j >= 12){
	      //cout << " entry: " << i << " ndatapmt: " << j << " ndata"  << " Pmt: " << pmtindx << " Col: " << col_arr[pmtindx] << " intg: " << intg[j] << " time: " << time[j] << endl;
		//}
	    //cout << " COLUMN CUT entry: " << i << " ndatapmt: " << ndata_passed_count << " " << same_x_repeated_x_freq[col_itr] << " ndata " << j << " Pmt: " << pmtindx1 << " Col: " << col_arr[pmtindx1] << " Row: " << row_arr[pmtindx1] << " xpos: " << xpos_arr[pmtindx1] << " intg: " << intg[j] << " time: " << time[j] << endl;


	      if(i < 10000){
	      if(entry_passed_list[i] == 0 && i < 10000){
		entry_passed_list[i] = same_x_repeated_x_freq[col_itr];
	      }
	      }
	      if(time[j] > 0 && time[j] < 512){//290 325
		h_pmt_is[pmtindx1]->Fill(amp_mv[j]);
		int_pmt_is[pmtindx1]->Fill(intg[j]);
		time_pmt_is[pmtindx1]->Fill(time[j]);

		//cout << same_x_repeated_x_freq[col_itr] << endl;
	
		ped_pmt_is[pmtindx1]->Fill(ped[j]);

		if(ndata_passed_count == same_x_repeated_x_freq[col_itr]){
		
		  hit_mult_per_entry[pmtindx1]->Fill(same_x_repeated_x_freq[col_itr]);
		  avg_intg_vs_hit_mult_per_entry[pmtindx1]->Fill(same_x_repeated_x_freq[col_itr], intg[j]);

		  col_histo[pmtindx1]->Fill(col_arr[pmtindx1]);
		  xpos_histo[pmtindx1]->Fill(xpos_arr[pmtindx1]);

		  if(same_x_repeated_x_freq[col_itr] == 5){
		    mult_count = mult_count + same_x_repeated_x_freq[col_itr];
		  }

		  
		  
		}
	      }
	      
	      
	    }
	    
	  }
	  
	}
	/*
	for(int list = 0; list < col_list.size(); list++){
	  col_histo[
	}
	*/
      }
    }
  
    //cout << mult_count << "  LOOK " << endl;








    














  }
    
  double avg_intg[1656] = {0.};
  
  for(int g = 0; g < 1656; g++){

    double tot = 0;

    if(pmt_intg_list[g].size() > 0){
      for(int h = 0; h < pmt_intg_list[g].size(); h++){
	tot += pmt_intg_list[g][h];
      }
      
      avg_intg[g] = (tot/pmt_intg_list[g].size());
      all_pmt_intg->Fill(g, avg_intg[g]);
      
    }

  }
  
  double avg_amp_mv[1656] = {0.};
  
  for(int g = 0; g < 1656; g++){

    double tot = 0;

    if(pmt_amp_mv_list[g].size() > 0){
      for(int h = 0; h < pmt_amp_mv_list[g].size(); h++){
	tot += pmt_amp_mv_list[g][h];
      }
      
      avg_amp_mv[g] = (tot/pmt_amp_mv_list[g].size());
      all_pmt_amp_mv->Fill(g, avg_amp_mv[g]);
    }

  }
  
  //cout << "Fraction of entries passing cuts: " << entry_pass_count << " / " << nentries << endl;


    gStyle->SetOptStat(0);
// Let's fit the histograms with Gauss (twice)
  TF1 *fgaus = new TF1("fgaus","gaus");
  TF1 *fgaus2 = new TF1("fgaus2","gaus");
  TF1 *fgaus3 = new TF1("fgaus3","gaus");
  TF1 *fgaus4 = new TF1("fgaus4","gaus");

      for(Int_t n=0; n<1656; n++){
      int lowerBinCamp;
      int upperBinCamp; 
      int lowerBinCtime; 
      int upperBinCtime; 
      int lowerBinCint; 
      int upperBinCint; 
      int lowerBinCped; 
      int upperBinCped;
      int maxBinamp = h_pmt_is[n]->GetMaximumBin();
      double maxBinCenteramp = h_pmt_is[n]->GetXaxis()->GetBinCenter( maxBinamp );
      double maxCountamp = h_pmt_is[n]->GetMaximum();
      double binWidthamp = h_pmt_is[n]->GetBinWidth(maxBinamp);
      double stdDevamp = h_pmt_is[n]->GetStdDev();
      int maxBintime = time_pmt_is[n]->GetMaximumBin();
      double maxBinCentertime = time_pmt_is[n]->GetXaxis()->GetBinCenter( maxBintime );
      double maxCounttime = time_pmt_is[n]->GetMaximum();
      double binWidthtime = time_pmt_is[n]->GetBinWidth(maxBintime);
      double stdDevtime = time_pmt_is[n]->GetStdDev();
      int maxBinint = int_pmt_is[n]->GetMaximumBin();
      double maxBinCenterint = int_pmt_is[n]->GetXaxis()->GetBinCenter( maxBinint );
      double maxCountint = int_pmt_is[n]->GetMaximum();
      double binWidthint = int_pmt_is[n]->GetBinWidth(maxBinint);
      double stdDevint = int_pmt_is[n]->GetStdDev();
      int maxBinped = ped_pmt_is[n]->GetMaximumBin();
      double maxBinCenterped = ped_pmt_is[n]->GetXaxis()->GetBinCenter( maxBinped );
      double maxCountped = ped_pmt_is[n]->GetMaximum();
      double binWidthped = ped_pmt_is[n]->GetBinWidth(maxBinped);
      double stdDevped = ped_pmt_is[n]->GetStdDev();

	//Amp fit
      if(h_pmt_is[n]->GetEntries()>1){ 

	// Create fit functions for each module
	fgaus->SetLineColor(2);
	fgaus->SetNpx(1000);

	// first fit
	lowerBinCamp = (maxBinamp)*binWidthamp - 3*stdDevamp;
	upperBinCamp = (maxBinamp)*binWidthamp + 3*stdDevamp;
	fgaus->SetParameters( maxCountamp,maxBinCenteramp,(stdDevamp*0.5) );
	fgaus->SetRange( lowerBinCamp, upperBinCamp );
	//h_pmt_is[n]->Fit( fgaus,"NO+RQ");
	//fgaus->GetParameters(Parsamp);

  
	// Second fit with tailored range
	//lowerBinCamp = hamp_min + (maxBinamp)*binWidthamp - Parsamp[2];
	//upperBinCamp = hamp_min + (maxBinamp)*binWidthamp + Parsamp[2];
	//fgaus->SetParameters( Parsamp[0],Parsamp[1],Parsamp[2] );
	//fgaus->SetRange( lowerBinCamp, upperBinCamp );

	//h_pmt_is[n]->Fit( fgaus,"+RQ" ); //kip
	fgaus->GetParameters(Parsamp);
        maxa[n] = Parsamp[0];
        fita[n] = Parsamp[1];
        stda[n] = Parsamp[2];
	for ( int i=0; i<3; i++ ) ParErrsamp[i] = fgaus->GetParError(i); 

	
	//file->cd();

}
	//Time fit
        if(time_pmt_is[n]->GetEntries()>1){ 

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

  
	// Second fit with tailored range
	//lowerBinCtime = hamp_min + (maxBintime)*binWidthtime - Parstime[2];
	//upperBinCtime = hamp_min + (maxBintime)*binWidthtime + Parstime[2];
	//fgaus2->SetParameters( Parstime[0],Parstime[1],Parstime[2] );
	//fgaus2->SetRange( lowerBinCtime, upperBinCtime );


	time_pmt_is[n]->Fit( fgaus2,"+RQ" );//kip
	time_mean[n] = fgaus2->GetParameter(1);
	fgaus2->GetParameters(Parstime);

	time_outlier_list[n] = 195 - time_mean[n];
	
	//if(time_pmt_is[n]->GetEntries() < 100 && abs(time_outlier_list[n] > 20)){
	if(abs(time_outlier_list[n] > 20)){
	  time_bad_list[n] = 1;
	}

	time_entries_list[n] = time_pmt_is[n]->GetEntries();
	
	
      maxt[n] = Parstime[0];
      fitt[n] = Parstime[1];
      stdt[n] = Parstime[2];
	for ( int i=0; i<3; i++ ) ParErrstime[i] = fgaus2->GetParError(i); 
	//file->cd();
             
}
        // Integral fit
        if(int_pmt_is[n]->GetEntries()>1){ 

	// Create fit functions for each module
	fgaus3->SetLineColor(2);
	fgaus3->SetNpx(1000);

	// first fit
	lowerBinCint = (maxBinint)*binWidthint - 3*stdDevint;
	upperBinCint = (maxBinint)*binWidthint + 3*stdDevint;
	fgaus3->SetParameters( maxCountint,maxBinCenterint,(stdDevint*0.5) );
	fgaus3->SetRange( lowerBinCint, upperBinCint );
	//int_pmt_is[n]->Fit(fgaus3,"NO+RQ");
	//fgaus3->GetParameters(Parsint);

  
	// Second fit with tailored range
	//lowerBinCint = hamp_min + (maxBinint)*binWidthint - Parsint[2];
	//upperBinCint = hamp_min + (maxBinint)*binWidthint + Parsint[2];
	//fgaus3->SetParameters( Parsint[0],Parsint[1],Parsint[2] );
	//fgaus3->SetRange( lowerBinCint, upperBinCint );

	//int_pmt_is[n]->Fit( fgaus3,"+RQ" ); //kip
	fgaus3->GetParameters(Parsint);
	maxi[n] = Parsint[0];
	fiti[n] = Parsint[1];
	stdi[n] = Parsint[2];
	for ( int i=0; i<3; i++ ) ParErrsint[i] = fgaus3->GetParError(i); 
	//file->cd();

}
        //Pedestal fit
       if(ped_pmt_is[n]->GetEntries()>1){
	// Create fit functions for each module
	fgaus4->SetLineColor(2);
	fgaus4->SetNpx(1000);

	
	// first fit
	lowerBinCped = (maxBinped)*binWidthped - 3*stdDevped;
	upperBinCped = (maxBinped)*binWidthped + 3*stdDevped;
	fgaus4->SetParameters( maxCountped,maxBinCenterped,stdDevped );
	fgaus4->SetRange( lowerBinCped, upperBinCped );
	//ped_pmt_is[n]->Fit(fgaus4,"NO+RQ");
	//fgaus4->GetParameters(Parsped);

  
	// Second fit with tailored range
	//lowerBinCped = hamp_min + (maxBinped)*binWidthped - Parsped[2];
	//upperBinCped = hamp_min + (maxBinped)*binWidthped + Parsped[2];
	//fgaus4->SetParameters( Parsped[0],Parsped[1],Parsped[2] );
	//fgaus4->SetRange( lowerBinCped, upperBinCped );

	//ped_pmt_is[n]->Fit( fgaus4,"+RQ" );//kip
	fgaus4->GetParameters(Parsped);
	maxp[n] = Parsped[0];
	fitp[n] = Parsped[1];
	stdp[n] = Parsped[2];
	for ( int i=0; i<3; i++ ) ParErrsped[i] = fgaus4->GetParError(i); 

	
	//file->cd();
     
}
}

   
    for(Int_t r=0; r < 69; r++){
     if(r<3){
    for(Int_t c=0; c<12; c++){
    h_rc[c][r] = h_pmt_is[c+r*12];
    time_rc[c][r] = time_pmt_is[c+r*12];
    int_rc[c][r] = int_pmt_is[c+r*12];
    ped_rc[c][r]= ped_pmt_is[c+r*12];}
    }
    else if(r<6){
    for(Int_t c=0; c<18; c++){
    h_rc[c][r] = h_pmt_is[36+c+(r-3)*18];
    time_rc[c][r] = time_pmt_is[36+c+(r-3)*18];
    int_rc[c][r] = int_pmt_is[36+c+(r-3)*18];
    ped_rc[c][r]= ped_pmt_is[36+c+(r-3)*18];}
    }
    else if(r<12){
    for(Int_t c=0; c<21; c++){
    h_rc[c][r] = h_pmt_is[90+c+(r-6)*21];
    time_rc[c][r] = time_pmt_is[90+c+(r-6)*21];
    int_rc[c][r] = int_pmt_is[90+c+(r-6)*21];
    ped_rc[c][r]= ped_pmt_is[90+c+(r-6)*21];}
    }
    else if(r<18){
    for(Int_t c=0; c<24; c++){
    h_rc[c][r] = h_pmt_is[216+c+(r-12)*24];
    time_rc[c][r] = time_pmt_is[36+c+(r-3)*18];
    int_rc[c][r] = int_pmt_is[36+c+(r-3)*18];
    ped_rc[c][r]= ped_pmt_is[36+c+(r-3)*18];}
    }
    else if(r<54){
    for(Int_t c=0; c<27; c++){
    h_rc[c][r] = h_pmt_is[360+c+(r-18)*27];
    time_rc[c][r] = time_pmt_is[360+c+(r-18)*27];
    int_rc[c][r] = int_pmt_is[360+c+(r-18)*27];
    ped_rc[c][r]= ped_pmt_is[360+c+(r-18)*27];}
    }
    else if(r<60){
    for(Int_t c=0; c<24; c++){
    h_rc[c][r] = h_pmt_is[1332+c+(r-54)*24];
    time_rc[c][r] = time_pmt_is[1332+c+(r-54)*24];
    int_rc[c][r] = int_pmt_is[1332+c+(r-54)*24];
    ped_rc[c][r]= ped_pmt_is[1332+c+(r-54)*24];}
    }
    else if(r<66){
    for(Int_t c=0; c<21; c++){
    h_rc[c][r] = h_pmt_is[1476+c+(r-60)*21];
    time_rc[c][r] = time_pmt_is[1476+c+(r-60)*21];
    int_rc[c][r] = int_pmt_is[1476+c+(r-60)*21];
    ped_rc[c][r]= ped_pmt_is[1476+c+(r-60)*21];}
    }
    else {
    for(Int_t c=0; c<18; c++){
    h_rc[c][r] = h_pmt_is[1602+c+(r-66)*18];
    time_rc[c][r] = time_pmt_is[1602+c+(r-66)*18];
    int_rc[c][r] = int_pmt_is[1602+c+(r-66)*18];
    ped_rc[c][r]= ped_pmt_is[1602+c+(r-66)*18];}
    }
  }
     for(Int_t r=0; r < 69; r++){
     if(r<3){
    for(Int_t c=0; c<12; c++){
    fa_rc[c][r] = fita[c+r*12];
    ft_rc[c][r] = fitt[c+r*12];
    fi_rc[c][r] = fiti[c+r*12];
    fp_rc[c][r] = fitp[c+r*12];}
    }
    else if(r<6){
    for(Int_t c=0; c<18; c++){
    fa_rc[c][r] = fita[36+c+(r-3)*18];
    ft_rc[c][r] = fitt[36+c+(r-3)*18];
    fi_rc[c][r] = fiti[36+c+(r-3)*18];
    fp_rc[c][r] = fitp[36+c+(r-3)*18];}
    }
    else if(r<12){
    for(Int_t c=0; c<21; c++){
    fa_rc[c][r] = fita[90+c+(r-6)*21];
    ft_rc[c][r] = fitt[90+c+(r-6)*21];
    fi_rc[c][r] = fiti[90+c+(r-6)*21];
    fp_rc[c][r] = fitp[90+c+(r-6)*21];}
    }
    else if(r<18){
    for(Int_t c=0; c<24; c++){
    fa_rc[c][r] = fita[216+c+(r-12)*24];
    ft_rc[c][r] = fitt[216+c+(r-12)*24];
    fi_rc[c][r] = fiti[216+c+(r-12)*24];
    fp_rc[c][r] = fitp[216+c+(r-12)*24];}
    }
    else if(r<54){
    for(Int_t c=0; c<27; c++){
    fa_rc[c][r] = fita[360+c+(r-18)*27];
    ft_rc[c][r] = fitt[360+c+(r-18)*27];
    fi_rc[c][r] = fiti[360+c+(r-18)*27];
    fp_rc[c][r] = fitp[360+c+(r-18)*27];}
    }
    else if(r<60){
    for(Int_t c=0; c<24; c++){
    fa_rc[c][r] = fita[1332+c+(r-54)*24];
    ft_rc[c][r] = fitt[1332+c+(r-54)*24];
    fi_rc[c][r] = fiti[1332+c+(r-54)*24];
    fp_rc[c][r] = fitp[1332+c+(r-54)*24];}
    }
    else if(r<66){
    for(Int_t c=0; c<21; c++){
    fa_rc[c][r] = fita[1476+c+(r-60)*21];
    ft_rc[c][r] = fitt[1476+c+(r-60)*21];
    fi_rc[c][r] = fiti[1476+c+(r-60)*21];
    fp_rc[c][r] = fitp[1476+c+(r-60)*21];}
    }
    else {
    for(Int_t c=0; c<18; c++){
    fa_rc[c][r] = fita[1602+c+(r-66)*18];
    ft_rc[c][r] = fitt[1602+c+(r-66)*18];
    fi_rc[c][r] = fiti[1602+c+(r-66)*18];
    fp_rc[c][r] = fitp[1602+c+(r-66)*18];}
    }
  }

 /*TCanvas* c1f[1656];
    for(Int_t n=0; n<1656; n++){
    c1f[n] = new TCanvas("","",900,900);
  }*/
     
     //TFile f(Form("runTEST/output.root"), "RECREATE");
     //col_vs_x->Write("col_vs_x");

     
     int cell_from_row_start;
     int cell_from_row_end;

     for(int cell_count = 0; cell_count < cells; cell_count++){
       
       if(row_arr[cell_count] == (row_start-1) && row_arr[cell_count-1] != (row_start-1) ){//row_start and row_end will be given as 1-69, so need to subtract 1 to get to 0-68 index
	 cell_from_row_start = cell_count;
       }

       if(row_end < 69 ){
	 
	 if( (row_arr[cell_count] == (row_end-1) && row_arr[cell_count+1] != (row_end-1)) ){
	   cell_from_row_end = cell_count+1;
	 }
	 
       }else{
	 cell_from_row_end = 1656;
       }
       
     }

     TString time_outlier_txt;
     
     time_outlier_txt = Form("time_outliers%i.txt", runis);
     
     ofstream time_bad_file(time_outlier_txt);

     time_bad_file << "Ch, nEntries, offset from avg, time peak(ns)" << endl;
     
     for(int i = 0; i < 1656; i++){
       if(time_bad_list[i] == 1){
	 time_bad_file << i << " " << time_entries_list[i] << " " << time_outlier_list[i] << " " << time_mean[i] << endl;
       }
     }
     
     for(int m = 0; m < 1655; m++){//1620-1625,

       //if(new_col[m]==col_start && row_arr[m] > 0){
       if(1 == 0){
       
     //for(int m =0; m < 1; m++){
    auto c1f=new TCanvas("c1f","",900,900);
    //for (Int_t m=0; m<1656; m++){
    //for(auto &m : pmt_list){
    //int y = specific_pmt;
    c1f->Clear();
    c1f->cd();
    c1f->Divide(2,4);
    gPad->SetBottomMargin(0.13);

    c1f->cd(1);
    h_pmt_is[m]->SetLabelSize(0.045, "x");
    h_pmt_is[m]->SetLabelSize(0.045, "y");
    h_pmt_is[m]->GetXaxis()->SetTitleSize(0.055);
    h_pmt_is[m]->GetYaxis()->SetTitleSize(0.055);
    h_pmt_is[m]->GetXaxis()->CenterTitle(1);
    h_pmt_is[m]->GetYaxis()->CenterTitle(1);
    h_pmt_is[m]->GetXaxis()->SetTitleOffset(1.1);
    h_pmt_is[m]->GetYaxis()->SetTitleOffset(1.1);
    h_pmt_is[m]->GetXaxis()->SetRangeUser(0,110);
    //h_pmt_is[m]->GetYaxis()->SetRangeUser(0,10);
    h_pmt_is[m]->SetTitle(Form("Amp PMT %i, %i (mV)",(row_arr[m]), (col_arr[m])));
    h_pmt_is[m]->SetLineColor(kPink+3);
    h_pmt_is[m]->SetLineWidth(2);
    h_pmt_is[m]->SetFillStyle(3144);
    h_pmt_is[m]->SetFillColor(kAzure+7);
    h_pmt_is[m]->Draw();
    //TLatex *t1 = new TLatex(20,1000,Form("fit = %.2f ",fita[m]));
    //t1->Draw();

    c1f->cd(2);
    time_pmt_is[m]->SetLabelSize(0.045, "x");
    time_pmt_is[m]->SetLabelSize(0.045, "y");
    time_pmt_is[m]->GetXaxis()->SetTitleSize(0.055);
    time_pmt_is[m]->GetYaxis()->SetTitleSize(0.055);
    time_pmt_is[m]->GetXaxis()->CenterTitle(1);
    time_pmt_is[m]->GetYaxis()->CenterTitle(1);
    time_pmt_is[m]->GetXaxis()->SetTitleOffset(1.1);
    time_pmt_is[m]->GetYaxis()->SetTitleOffset(1.1);
    time_pmt_is[m]->GetXaxis()->SetRangeUser(0,400);
    //time_pmt_is[m]->GetYaxis()->SetRangeUser(0,1.1*maxt[m]);
   // gPad->SetLogy(1);
    time_pmt_is[m]->SetTitle(Form("Time PMT %i %i (ns)",(row_arr[m]), (col_arr[m])));
    time_pmt_is[m]->SetLineColor(kPink+3);
    time_pmt_is[m]->SetLineWidth(2);
    time_pmt_is[m]->SetFillStyle(3144);
    time_pmt_is[m]->SetFillColor(kViolet+7);
    time_pmt_is[m]->Draw();
    TLatex *t2 = new TLatex(250,300,Form("fit = %.2f",fitt[m]));
    t2->Draw();
    
    c1f->cd(3);
    int_pmt_is[m]->SetLabelSize(0.045, "x");
    int_pmt_is[m]->SetLabelSize(0.045, "y");
    int_pmt_is[m]->GetXaxis()->SetTitleSize(0.055);
    int_pmt_is[m]->GetYaxis()->SetTitleSize(0.055);
    int_pmt_is[m]->GetXaxis()->CenterTitle(1);
    int_pmt_is[m]->GetYaxis()->CenterTitle(1);
    int_pmt_is[m]->GetXaxis()->SetTitleOffset(1.1);
    int_pmt_is[m]->GetYaxis()->SetTitleOffset(1.1);
    int_pmt_is[m]->GetXaxis()->SetRangeUser(0,100);
    
    //int_pmt_is[m]->GetYaxis()->SetRangeUser(0);
    int_pmt_is[m]->SetTitle(Form("Int PMT %i, %i (pC)",(row_arr[m]), (col_arr[m])));
    int_pmt_is[m]->SetLineColor(kPink+3);
    int_pmt_is[m]->SetLineWidth(2);
    int_pmt_is[m]->SetFillStyle(3144);
    int_pmt_is[m]->SetFillColor(kOrange+6);
    gStyle->SetOptStat(11111111);
    int_pmt_is[m]->Draw();




    string intcheck = Form("Pedestal PMT %i (mV)",m);
    //int_pmt_is[m]->Write(Form("int_pmt_is%i",m));
   // TLatex *t3 = new TLatex(80,300,Form("fit = %.2f",fiti[m]));
   // t3->Draw();
    
    c1f->cd(4);
    ped_pmt_is[m]->SetLabelSize(0.045, "x");
    ped_pmt_is[m]->SetLabelSize(0.045, "y");
    ped_pmt_is[m]->GetXaxis()->SetTitleSize(0.055);
    ped_pmt_is[m]->GetYaxis()->SetTitleSize(0.055);
    ped_pmt_is[m]->GetXaxis()->CenterTitle(1);
    ped_pmt_is[m]->GetYaxis()->CenterTitle(1);
    ped_pmt_is[m]->GetXaxis()->SetTitleOffset(1.1);
    ped_pmt_is[m]->GetYaxis()->SetTitleOffset(1.1);
    ped_pmt_is[m]->GetXaxis()->SetRangeUser(fitp[m]-30*stdp[m],fitp[m]+30*stdp[m]);
    //ped_pmt_is[m]->GetYaxis()->SetRangeUser(0,1.1*maxp[m]);
    ped_pmt_is[m]->SetTitle(Form("Pedestal PMT %i, %i (mV)",(row_arr[m]), (col_arr[m])));
    ped_pmt_is[m]->SetLineColor(kPink+3);
    ped_pmt_is[m]->SetLineWidth(2);
    ped_pmt_is[m]->SetFillStyle(3144);
    ped_pmt_is[m]->SetFillColor(kGreen+1);
    ped_pmt_is[m]->Draw();
    //TLatex *t4 = new TLatex(40,6000,Form("fit = %.2f ",fitp[m]));
    //t4->Draw();

    c1f->cd(5);
    hit_mult_per_entry[m]->SetLabelSize(0.045, "x");
    hit_mult_per_entry[m]->SetLabelSize(0.045, "y");
    hit_mult_per_entry[m]->GetXaxis()->SetTitleSize(0.055);
    hit_mult_per_entry[m]->GetYaxis()->SetTitleSize(0.055);
    hit_mult_per_entry[m]->GetXaxis()->CenterTitle(1);
    hit_mult_per_entry[m]->GetYaxis()->CenterTitle(1);
    hit_mult_per_entry[m]->GetXaxis()->SetTitleOffset(1.1);
    hit_mult_per_entry[m]->GetYaxis()->SetTitleOffset(1.1);
    hit_mult_per_entry[m]->GetXaxis()->SetRangeUser(fitp[m]-30*stdp[m],fitp[m]+30*stdp[m]);
    //ped_pmt_is[m]->GetYaxis()->SetRangeUser(0,1.1*maxp[m]);
    hit_mult_per_entry[m]->SetTitle(Form("Pulse Recon. Mult. per Entry, PMT %i, %i (mV)",(row_arr[m]), (col_arr[m])));
    hit_mult_per_entry[m]->SetLineColor(kPink+4);
    hit_mult_per_entry[m]->SetLineWidth(2);
    hit_mult_per_entry[m]->SetFillStyle(3144);
    hit_mult_per_entry[m]->SetFillColor(kGreen+2);
    hit_mult_per_entry[m]->Draw();

    c1f->cd(6);
    avg_intg_vs_hit_mult_per_entry[m]->SetLabelSize(0.045, "x");
    avg_intg_vs_hit_mult_per_entry[m]->SetLabelSize(0.045, "y");
    avg_intg_vs_hit_mult_per_entry[m]->GetXaxis()->SetTitleSize(0.055);
    avg_intg_vs_hit_mult_per_entry[m]->GetYaxis()->SetTitleSize(0.055);
    avg_intg_vs_hit_mult_per_entry[m]->GetXaxis()->CenterTitle(1);
    avg_intg_vs_hit_mult_per_entry[m]->GetYaxis()->CenterTitle(1);
    avg_intg_vs_hit_mult_per_entry[m]->GetXaxis()->SetTitleOffset(1.1);
    avg_intg_vs_hit_mult_per_entry[m]->GetYaxis()->SetTitleOffset(1.1);
    avg_intg_vs_hit_mult_per_entry[m]->GetXaxis()->SetRangeUser(fitp[m]-30*stdp[m],fitp[m]+30*stdp[m]);
    //ped_pmt_is[m]->GetYaxis()->SetRangeUser(0,1.1*maxp[m]);
    avg_intg_vs_hit_mult_per_entry[m]->SetTitle(Form("Intg vs Pulse Recon. Mult. per Entry, PMT %i, %i (mV)",(row_arr[m]+1), (col_arr[m]+1)));
    avg_intg_vs_hit_mult_per_entry[m]->SetLineColor(kPink+5);
    avg_intg_vs_hit_mult_per_entry[m]->SetLineWidth(2);
    avg_intg_vs_hit_mult_per_entry[m]->SetFillStyle(3144);
    avg_intg_vs_hit_mult_per_entry[m]->SetFillColor(kGreen+3);
    avg_intg_vs_hit_mult_per_entry[m]->Draw();

    c1f->cd(7);
    col_histo[m]->SetLabelSize(0.045, "x");
    col_histo[m]->SetLabelSize(0.045, "y");
    col_histo[m]->GetXaxis()->SetTitleSize(0.055);
    col_histo[m]->GetYaxis()->SetTitleSize(0.055);
    col_histo[m]->GetXaxis()->CenterTitle(1);
    col_histo[m]->GetYaxis()->CenterTitle(1);
    col_histo[m]->GetXaxis()->SetTitleOffset(1.1);
    col_histo[m]->GetYaxis()->SetTitleOffset(1.1);
    col_histo[m]->GetXaxis()->SetRangeUser(fitp[m]-30*stdp[m],fitp[m]+30*stdp[m]);
    //ped_pmt_is[m]->GetYaxis()->SetRangeUser(0,1.1*maxp[m]);
    col_histo[m]->SetTitle(Form("Column index, PMT %i, %i (mV)",(row_arr[m]), (col_arr[m])));
    col_histo[m]->SetLineColor(kPink+7);
    col_histo[m]->SetLineWidth(2);
    col_histo[m]->SetFillStyle(3144);
    col_histo[m]->SetFillColor(kGreen+0);
    col_histo[m]->Draw();

    c1f->cd(8);
    xpos_histo[m]->SetLabelSize(0.045, "x");
    xpos_histo[m]->SetLabelSize(0.045, "y");
    xpos_histo[m]->GetXaxis()->SetTitleSize(0.055);
    xpos_histo[m]->GetYaxis()->SetTitleSize(0.055);
    xpos_histo[m]->GetXaxis()->CenterTitle(1);
    xpos_histo[m]->GetYaxis()->CenterTitle(1);
    xpos_histo[m]->GetXaxis()->SetTitleOffset(1.1);
    xpos_histo[m]->GetYaxis()->SetTitleOffset(1.1);
    xpos_histo[m]->GetXaxis()->SetRangeUser(fitp[m]-30*stdp[m],fitp[m]+30*stdp[m]);
    //ped_pmt_is[m]->GetYaxis()->SetRangeUser(0,1.1*maxp[m]);
    xpos_histo[m]->SetTitle(Form("X Position, PMT %i, %i (mV)",(row_arr[m]), (col_arr[m])));
    xpos_histo[m]->SetLineColor(kPink+1);
    xpos_histo[m]->SetLineWidth(2);
    xpos_histo[m]->SetFillStyle(3144);
    xpos_histo[m]->SetFillColor(kGreen+9);
    xpos_histo[m]->Draw();
    
    //c1f->Print(Form("run%i/Run_%i_PMT_%i.png",runis,runis,m));
    //c1f->Print(Form("run%i/Run_TEST_PMT_%i.png",runis, m));
     }
     
/*
double x[1656];
     Int_t nPMT = 1656;
     for (int i=0; i<nPMT; i++){
      x[i] = i;
}
     
    auto *Amp_pmt = new TGraph{1656,x,fita};
    auto *time_pmt = new TGraph{1656,x,fitt};
    auto *intg_pmt = new TGraph{1656,x,fiti};
    auto *ped_pmt = new TGraph{1656,x,fitp};
    auto c1a=new TCanvas("c1a","",900,900);

    c1a->Clear();
    c1a->Divide(2,2);
    gPad->SetBottomMargin(0.13);

    c1a->cd(1);
    Amp_pmt->GetXaxis()->SetTitleSize(0.055);
    Amp_pmt->GetYaxis()->SetTitleSize(0.055);
    Amp_pmt->GetXaxis()->CenterTitle(1);
    Amp_pmt->GetYaxis()->CenterTitle(1);
    Amp_pmt->GetXaxis()->SetTitleOffset(1.1);
    Amp_pmt->GetYaxis()->SetTitleOffset(1.1);
    Amp_pmt->SetTitle("Amplitude - PMT (mV)");
    Amp_pmt->SetMarkerStyle(20);
    Amp_pmt->GetYaxis()->SetRangeUser(0,500);
    Amp_pmt->SetMarkerColor(kAzure+7);
    Amp_pmt->Draw("AP");
    
    c1a->cd(2);
    time_pmt->GetXaxis()->SetTitleSize(0.055);
    time_pmt->GetYaxis()->SetTitleSize(0.055);
    time_pmt->GetXaxis()->CenterTitle(1);
    time_pmt->GetYaxis()->CenterTitle(1);
    time_pmt->GetXaxis()->SetTitleOffset(1.1);
    time_pmt->GetYaxis()->SetTitleOffset(1.1);
    time_pmt->SetTitle("Time - PMT (ns)");
    time_pmt->SetMarkerStyle(20);
    time_pmt->SetMarkerColor(kViolet+7);
    time_pmt->GetYaxis()->SetRangeUser(0,360);
    time_pmt->Draw("AP");  

    c1a->cd(3);
    intg_pmt->GetXaxis()->SetTitleSize(0.055);
    intg_pmt->GetYaxis()->SetTitleSize(0.055);
    intg_pmt->GetXaxis()->CenterTitle(1);
    intg_pmt->GetYaxis()->CenterTitle(1);
    intg_pmt->GetXaxis()->SetTitleOffset(1.1);
    intg_pmt->GetYaxis()->SetTitleOffset(1.1);
    intg_pmt->SetTitle("Integral - PMT (pC)");
    intg_pmt->SetMarkerStyle(20);
    intg_pmt->GetYaxis()->SetRangeUser(0,600);
    intg_pmt->SetMarkerColor(kOrange+7);
    intg_pmt->Draw("AP");



    c1a->cd(4);
    ped_pmt->GetXaxis()->SetTitleSize(0.055);
    ped_pmt->GetYaxis()->SetTitleSize(0.055);
    ped_pmt->GetXaxis()->CenterTitle(1);
    ped_pmt->GetYaxis()->CenterTitle(1);
    ped_pmt->GetXaxis()->SetTitleOffset(1.1);
    ped_pmt->GetYaxis()->SetTitleOffset(1.1);
    ped_pmt->SetTitle("Pedestal - PMT (mV)");
    ped_pmt->SetMarkerStyle(20);
    ped_pmt->GetYaxis()->SetRangeUser(0,150);
    ped_pmt->SetMarkerColor(kGreen+1);
    ped_pmt->Draw("AP");
    //c1a->Print(Form("run%i/Run_%i_fits.png",runis,runis));

    double xx[69];
     for (int i=0; i<69; i++){
      xx[i] = i;
}
    TGraph* Amp_c[27];
    TGraph* time_c[27];
    TGraph* intg_c[27];
    TGraph* ped_c[27];
    TCanvas* c1a[27];
    
     for(Int_t nc=0; nc<27; nc++){
       Amp_c[nc] = new TGraph{69,xx,fa_rc[nc]};
       time_c[nc] = new TGraph{69,xx,ft_rc[nc]};
       intg_c[nc] = new TGraph{69,xx,fi_rc[nc]};
       ped_c[nc] = new TGraph{69,xx,fp_rc[nc]};
       c1a[nc]; = new TCanvas("","",900,900);

    c1a[nc]->Clear();
    c1a[nc]->Divide(2,2);
    gPad->SetBottomMargin(0.13);

    c1a[nc]->cd(1);
    Amp_c[nc]->GetXaxis()->SetTitleSize(0.055);
    Amp_c[nc]->GetYaxis()->SetTitleSize(0.055);
    Amp_c[nc]->GetXaxis()->CenterTitle(1);
    Amp_c[nc]->GetYaxis()->CenterTitle(1);
    Amp_c[nc]->GetXaxis()->SetTitleOffset(1.1);
    Amp_c[nc]->GetYaxis()->SetTitleOffset(1.1);
    Amp_c[nc]->SetTitle(Form("Amplitude - col %d (mV)",2));
    Amp_c[nc]->SetMarkerStyle(20);
    Amp_c[nc]->GetYaxis()->SetRangeUser(0,500);
    Amp_c[nc]->SetMarkerColor(kAzure+7);
    Amp_c[nc]->Draw("AP");
    
    c1a[nc]->cd(2);
    time_c[nc]->GetXaxis()->SetTitleSize(0.055);
    time_c[nc]->GetYaxis()->SetTitleSize(0.055);
    time_c[nc]->GetXaxis()->CenterTitle(1);
    time_c[nc]->GetYaxis()->CenterTitle(1);
    time_c[nc]->GetXaxis()->SetTitleOffset(1.1);
    time_c[nc]->GetYaxis()->SetTitleOffset(1.1);
    time_c[nc]->SetTitle("Time - PMT (ns)");
    time_c[nc]->SetMarkerStyle(20);
    time_c[nc]->SetMarkerColor(kViolet+7);
    time_c[nc]->GetYaxis()->SetRangeUser(0,360);
    time_c[nc]->Draw("AP");  

    c1a[nc]->cd(3);
    intg_c[nc]->GetXaxis()->SetTitleSize(0.055);
    intg_c[nc]->GetYaxis()->SetTitleSize(0.055);
    intg_c[nc]->GetXaxis()->CenterTitle(1);
    intg_c[nc]->GetYaxis()->CenterTitle(1);
    intg_c[nc]->GetXaxis()->SetTitleOffset(1.1);
    intg_c[nc]->GetYaxis()->SetTitleOffset(1.1);
    intg_c[nc]->SetTitle("Integral - PMT (pC)");
    intg_c[nc]->SetMarkerStyle(20);
    intg_c[nc]->GetYaxis()->SetRangeUser(0,600);
    intg_c[nc]->SetMarkerColor(kOrange+7);
    intg_c[nc]->Draw("AP");

    c1a[nc]->cd(4);
    ped_pmt[nc]->GetXaxis()->SetTitleSize(0.055);
    ped_pmt[nc]->GetYaxis()->SetTitleSize(0.055);
    ped_pmt[nc]->GetXaxis()->CenterTitle(1);
    ped_pmt[nc]->GetYaxis()->CenterTitle(1);
    ped_pmt[nc]->GetXaxis()->SetTitleOffset(1.1);
    ped_pmt[nc]->GetYaxis()->SetTitleOffset(1.1);
    ped_pmt[nc]->SetTitle("Pedestal - PMT (mV)");
    ped_pmt[nc]->SetMarkerStyle(20);
    ped_pmt[nc]->GetYaxis()->SetRangeUser(0,150);
    ped_pmt[nc]->SetMarkerColor(kGreen+1);
    ped_pmt[nc]->Draw("AP");
    c1a[nc]->Print(Form("run%i/Run_%i_col_%i_fits.png",runis,runis,nc));
     }*/
/*
     TString rundir = Form("run
     

  gSystem->FreeDirectory(dir);
*/
     }
     /*
     TString title;
     title.Form("passed_col_cut_list_%d.txt", mid_multiplicity);
     ofstream col_mult_list(title);
     for(int l = 0; l < 10000; l++){
       col_mult_list << l << "   " << entry_passed_list[l] << endl;
     }
     */

     

  TString conv_comm,rm_comm;
  conv_comm = Form("convert -density 100x100 -quality 60 -compress jpeg run%d/*.png run%d/combined_run_%d_Col_%d_Mult%d.pdf",runis,runis,runis,col_start,mid_multiplicity);
  rm_comm = Form("rm run%d/Run_%d_PMT*.png",runis,runis);
  //system(conv_comm);
  //system(rm_comm);

  ofstream time_mean_list;
  time_mean_list.open(Form("time_peak_list_run%i.txt",runis));

  
  
  for(int n=0; n < 1656; n++){
    time_mean_list << n << " " << time_mean[n] << endl;
    time_peak_vs_ch->Fill(n, time_mean[n]);
  }

  time_mean_list.close();

  TFile f(Form("timing_distribution_%i.root",runis), "RECREATE");
  time_peak_vs_ch->Write("time_peak_vs_ch_mult5");
  
  gSystem->Exit(0);
}

