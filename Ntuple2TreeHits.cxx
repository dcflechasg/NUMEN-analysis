#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TVector3.h"
#include <TFile.h>
#include <TNtuple.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <vector>

void Ntuple2Tree(const char* Ntuplefile, const char* Treefile, Int_t Nevs=0, Int_t MaxNhits=1000) {
	
	TFile *myFile1 = TFile::Open(Ntuplefile);
	if(myFile1->IsZombie()){cout<<"Problem opening Ntuple file"<<endl; return;}
	TFile outputFile (Treefile,"RECREATE");
	if(outputFile.IsZombie()){cout<<"Problem opening output file"<<endl; return;}
	TTree *tree= new TTree("tree","tree");
	
	Int_t NEvents; // event counter
	Int_t MaxNhitObs=0;

/*
	Int_t Fold; // number of detectors fired in event
	Int_t *DetectorID = new Int_t [Ndets]; // ID of detector
	Double_t *Energy = new Double_t [Ndets]; // total energy in each det
	Double_t *Time = new Double_t [Ndets]; // shortest time in det
*/

	Int_t Nhit;
	Int_t *DetIDhit=new Int_t [MaxNhits];
	Double_t *Ehit=new Double_t [MaxNhits];
	Double_t *Thit=new Double_t [MaxNhits];
	Int_t *PDGchit=new Int_t [MaxNhits];
	Double_t ExDCE=0.;
	Int_t NRhit=0;
	Int_t kR=0;

/*
 	vector <Int_t> IDvec;
	vector <Double_t> Evec;
	vector <Double_t> Tvec;
	vector <Int_t>::iterator it;
	
	tree->Branch("Fold",&Fold,"Fold/I");
	tree->Branch("DetectorID",DetectorID,"DetectorID[Fold]/I");
	tree->Branch("Energy",Energy,"Energy[Fold]/D");
	tree->Branch("Time",Time,"Time[Fold]/D");
*/
	
	tree->Branch("Nhit",&Nhit,"Nhit/I");
	tree->Branch("DetIDhit",DetIDhit,"DetIDhit[Nhit]/I");
	tree->Branch("Ehit",Ehit,"Ehit[Nhit]/D");
	tree->Branch("Thit",Thit,"Thit[Nhit]/D");
	tree->Branch("PDGchit",PDGchit,"PDGchit[Nhit]/I");
	tree->Branch("ExDCE",&ExDCE,"ExDCE/D");
	tree->Branch("NRhit",&NRhit,"NRhit/I");
	tree->Branch("kR",&kR,"kR/I");
	
//    TNtuple *Ntuple = new TNtuple("Ntuple","Combined Ntuple","Event:PDGcode:DetectorID:Energy:Time");

   TTreeReader myReader1("Ntuple", myFile1);
   TTreeReaderValue<Int_t> myEvent1(myReader1, "Event");
   TTreeReaderValue<Int_t> myPDGcode1(myReader1, "PDGcode");
   TTreeReaderValue<Int_t> myDetectorID1(myReader1, "DetectorID");
   TTreeReaderValue<Double_t> myEnergy1(myReader1, "Energy");
   TTreeReaderValue<Double_t> myTime1(myReader1, "Time");

   myReader1.Restart();

   // Loop over all entries of the Ntuples
   
   bool OK1=myReader1.Next();
   Int_t ThisEvent1;
   NEvents=0;
   while (OK1){
     ThisEvent1=*myEvent1;
     NEvents++;
     if((Nevs>0) && (NEvents>Nevs))break; 
     Nhit=0;
/*
     IDvec.clear();
     Evec.clear();
     Tvec.clear();
*/

     while(ThisEvent1==*myEvent1){
		 
		 if(MaxNhitObs>=MaxNhits){cout<<"Error: MaxNhit Observed >= MaxNhits hits"<<endl;break;}
		 DetIDhit[Nhit]=*myDetectorID1;
		 Ehit[Nhit]=*myEnergy1;
		 Thit[Nhit]=*myTime1;
		 PDGchit[Nhit]=*myPDGcode1;
		 Nhit++;
		 if(Nhit>MaxNhitObs)MaxNhitObs=Nhit;
/*		 
		 it = find (IDvec.begin(), IDvec.end(), *myDetectorID);
		 if (it == IDvec.end()){
			 IDvec.push_back(*myDetectorID);
			 Evec.push_back(*myEnergy);
			 Tvec.push_back(*myTime);
		 }
		 else {
			 int d=distance(IDvec.begin(),it);
			 Evec[d]+=*myEnergy;		 
			 if(*myTime<Tvec[d]) Tvec[d]=*myTime;
		 }
*/		 
		 
		 OK1=myReader1.Next();
		 if(!OK1) break;
    }		 
    
/*
 		 for (it = IDvec.begin() ; it != IDvec.end(); ++it){
			 DetectorID[id]=IDvec[id];
			 Energy[id]=Evec[id];
			 Time[id]=Tvec[id];
		 }
		 Fold=IDvec.size();	 
*/		 
		 tree->Fill();
 
   }
	tree->Write();

	outputFile.Close();
	myFile1->Close();
	cout<<"Number of events found: "<<NEvents<<endl;
	cout<<"Max Nhit= "<<MaxNhitObs<<endl;
	return;
} 


