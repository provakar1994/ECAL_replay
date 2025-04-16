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

#include <math.h>

void make_gain_db_list(){








  ifstream ecalfile_a("PMTalphas1300_1350.txt");

  TString currentline_a;

  int arr_index_a = 0;

  double alpha_arr[1656];
  
  while( currentline_a.ReadLine(ecalfile_a) ){

    //if(){ //new cell information:
      TObjArray *tokens = ( (TObjArray*) currentline_a.Tokenize(" ") );
      int ntokens = tokens->GetEntries();

      if( ntokens == 2 ){
	int cell = ( (TObjString*) (*tokens)[0] )->GetString().Atoi();
	double alpha = ( (TObjString*) (*tokens)[1] )->GetString().Atof();

	//cout << alpha << endl;

	alpha_arr[arr_index_a] = alpha;
	arr_index_a++;

      }
  }












  ifstream ecalfile_rand("textfile/ecal_mc.dat");

  TString currentline_rand;

  int arr_index_rand = 0;

  double rand_arr[1656];
  
  while( currentline_rand.ReadLine(ecalfile_rand) ){

    //if(){ //new cell information:
      TObjArray *tokens = ( (TObjArray*) currentline_rand.Tokenize(" ") );
      int ntokens = tokens->GetEntries();

      if( ntokens == 1656 ){
	for(int i = 0; i < 1656; i++){
	  rand_arr[i] = ( (TObjString*) (*tokens)[i] )->GetString().Atof();
	  cout << rand_arr[i] << endl;
	}
	//rand_arr[arr_index_rand] = rand;

      }
  }










  


  //ifstream ecalfile1("textfile/ECal_run_2501_energy_ratio.txt");

  //ifstream ecalfile1("textfile/ECal_run_2501_same_gain_energy_ratio.txt");

  ifstream ecalfile1("textfile/ECal_run_2501_right_range_gain_no_fix_factor.txt");

  TString currentline1;

  int check1 = 0;

  int arr_index1 = 0;

  double energy_arr[1656];
  
  while( currentline1.ReadLine(ecalfile1) ){

    //if(){ //new cell information:
      TObjArray *tokens = ( (TObjArray*) currentline1.Tokenize(" ") );
      int ntokens = tokens->GetEntries();

      if( ntokens == 3 ){
	int cell = ( (TObjString*) (*tokens)[0] )->GetString().Atoi();
	double gain = ( (TObjString*) (*tokens)[2] )->GetString().Atof();

	//energy_arr[arr_index1] = gain*pow((0.007)/(rand_arr[arr_index1]), (1.0/alpha_arr[arr_index1]));
	//energy_arr[arr_index1] = gain*pow((rand_arr[arr_index1]*(1.0/0.007)), (1.0/alpha_arr[arr_index1]));
	//energy_arr[arr_index1] = gain*pow((0.007)/(rand_arr[arr_index1]), (1.0));
	//energy_arr[arr_index1] = gain*(rand_arr[arr_index1]*(1.0/0.007));
	//energy_arr[arr_index1] = gain*(rand_arr[arr_index1]*10000);
	energy_arr[arr_index1] = gain;
	//cout << gain << endl;

	arr_index1++;

      }
  }
  check1 = 1;






  

  ifstream ecalfile("textfile/ECal_run_2501_gain_factor.txt");
  //ifstream ecalfile("textfile/ECal_run_2501_random_gain_fit_factor.txt");

  TString currentline;

  int check = 0;

  int arr_index = 0;

  double gain_arr[1656];
  double gain_arr1[1656];
  double gain_arr2[1656];
  double gain_arr1a[1656];
  double gain_arr2a[1656];
  double gain_arr3[1656];
  double gain_arr4[1656];
  
  while( currentline.ReadLine(ecalfile) ){

    //if(){ //new cell information:
      TObjArray *tokens = ( (TObjArray*) currentline.Tokenize(" ") );
      int ntokens = tokens->GetEntries();

      if( ntokens == 2 ){
	int cell = ( (TObjString*) (*tokens)[0] )->GetString().Atoi();
	double gain = ( (TObjString*) (*tokens)[1] )->GetString().Atof();

	//cout << gain << endl;

	gain_arr1[arr_index] = gain*pow((0.007)/(rand_arr[arr_index]), (1.0/alpha_arr[arr_index]));
	gain_arr2[arr_index] = gain*pow((rand_arr[arr_index]*(1.0/0.007)), (1.0/alpha_arr[arr_index]));
	
	gain_arr1a[arr_index] = gain*pow((0.007)/(rand_arr[arr_index]), (alpha_arr[arr_index]));
	gain_arr2a[arr_index] = gain*pow((rand_arr[arr_index]*(1.0/0.007)), (alpha_arr[arr_index]));
	
	gain_arr3[arr_index] = gain*(0.007)/(rand_arr[arr_index]);
	gain_arr4[arr_index] = gain*(rand_arr[arr_index]*(1.0/0.007));

	//cout << rand_arr[arr_index]/0.007 << endl;
	
	gain_arr[arr_index] = gain;

	arr_index++;

      }
    }
    check = 1;

    int energy_count = 0;
    // }

    TH1F* gain_coeff = new TH1F("gain_coeff","", 1000, 0.0, 0.1);
    TH1F* gain_coeff1 = new TH1F("gain_coeff1","", 1000, 0.0, 0.1);
    TH1F* gain_coeff2 = new TH1F("gain_coeff2","", 1000, 0.0, 0.1);

    TH1F* gain_coeff1a = new TH1F("gain_coeff1a","", 1000, 0.0, 0.1);
    TH1F* gain_coeff2a = new TH1F("gain_coeff2a","", 1000, 0.0, 0.1);
    
    TH1F* gain_coeff3 = new TH1F("gain_coeff3","", 1000, 0.0, 0.1);
    TH1F* gain_coeff4 = new TH1F("gain_coeff4","", 1000, 0.0, 0.1);
    TH1F* energy_coeff = new TH1F("energy_coeff","", 1000, 0.0, 2.0);
    TH1F* energy_coeff_ind = new TH1F("energy_coeff_ind","", 1656, 0.5, 1656.5);
    TH1F* rand_coeff_ind = new TH1F("rand_coeff_ind","", 1656, 0.5, 1656.5);

    for(int i = 0; i < 1656; i++){

      gain_coeff->Fill(gain_arr[i]);
      gain_coeff1->Fill(gain_arr1[i]);
      gain_coeff2->Fill(gain_arr2[i]);

      gain_coeff1a->Fill(gain_arr1a[i]);
      gain_coeff2a->Fill(gain_arr2a[i]);
      
      gain_coeff3->Fill(gain_arr3[i]);
      gain_coeff4->Fill(gain_arr4[i]);
      
      if(gain_arr[i] >= 0.019){
	gain_arr[i] = 0.01;
      }
      //cout << gain_arr[i] << " ";
      energy_coeff->Fill(energy_arr[i]);
      energy_coeff_ind->Fill(i, energy_arr[i]);
      rand_coeff_ind->Fill(i, rand_arr[i]);

      cout << rand_arr[i] << endl;

      if(energy_arr[i] <= 1.25 && energy_arr[i] >= 0.75){
	energy_count++;
      }
    }

    TFile f1(Form("hist/gain_coeff_dist.root"), "RECREATE");
    gain_coeff->Write("gain_coeff");
    gain_coeff1->Write("gain_coeff1");
    gain_coeff2->Write("gain_coeff2");
    gain_coeff1a->Write("gain_coeff1a");
    gain_coeff2a->Write("gain_coeff2a");
    gain_coeff3->Write("gain_coeff3");
    gain_coeff4->Write("gain_coeff4");
    gain_coeff->Write("gain_coeff");
    energy_coeff->Write("energy_coeff");
    energy_coeff_ind->Write("energy_coeff_ind");
    rand_coeff_ind->Write("rand_coeff_ind");

    //cout << endl << "good energy" << energy_count << endl;



    
    
}
