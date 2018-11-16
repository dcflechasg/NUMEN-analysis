#include "TVector3.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "functions.cxx"



class TreeTrack{
	private:
		TTree *tree;
		TTree *tracked;	
		TFile* InFile;		
		TFile* OutFile;
		
		// in tree event specs
			Int_t Nin;
			Int_t Din[1000];
			Double_t Ein[1000];
			Double_t Tin[1000];
			Int_t Pin[1000];
			Double_t Exin;				
			Int_t NRin;
			Int_t kRin;
					
		// gamma array specs
		Int_t Ndets;
//		vector <TVector3> ArrayList;
		Double_t TimeResolution=1.;

		
		// event data summary 
		vector <Int_t> DetID;
		vector <Double_t> E;
		vector <Double_t> T;
		
		// tracking event summary
		Int_t Ngam;
		vector <Double_t> Egam;
		vector <Double_t> Tgam;
		vector <TVector3> Vgam;
		
		// tracked tree event specs
			Int_t Ng;			
			Double_t Eg[400];
			Double_t Tg[400];
			Double_t Thg[400];
			Double_t Phg[400];
			Double_t Ex;	
			Int_t kR;			

		void TrackGammas_9pix();
		void TrackGammas_1pix();
		//void TrackGammas_NeighbourCouple();
		void TrackGammas(Int_t ttype);
		Int_t AddPixelEnergies();
		void BlurrData(Int_t npix);
		void ClearVecs();// clear pixel data
		bool HitFilter(Int_t ihit);
		void OutputEvent();
		void CloseFiles();
		Int_t Ntracked=0;// number of tracked events
							
	public:
		Int_t DetectorPosTable(const char* nametab);
		void OpenTree(const char* namefile);
		void OpenTrackedFile(const char* namefile);
		void ProcessTree(Int_t TrackType);
		void ProcessTreeTest(Int_t nevlim, Int_t iev0);
		void ListEvent();
		Int_t DetType=1;//1=LYSO
		// 5keV (HitFilter parameters)
		Double_t Tfilter=50.;
		Double_t Ethresh=0.005;
		Double_t GetTheta(Int_t detID);
		Double_t GetPhi(Int_t detID);
		vector <TVector3> ArrayList;
};

void TreeTrack::TrackGammas(Int_t ttype){
	switch (ttype){
		case 1:
			TrackGammas_1pix();
			break;
		case 9:
			TrackGammas_9pix();
			break;
		default:
			TrackGammas_1pix();
		}
	return;
}

void TreeTrack::TrackGammas_1pix(){

	Int_t npix=AddPixelEnergies();// number of fired pixels // totl E per pix is added
	Ngam=npix; // no tracking, actually
	Egam.clear();
	Tgam.clear();
	Vgam.clear();		
	
	vector <Int_t> ClusterList;
	vector <Int_t>::iterator it;
	for (Int_t ipix=0;ipix<npix;ipix++){ // loop on fired pixels
		Int_t ClusterID=DetID[ipix];// each pix is a "cluster"
		Int_t icenter=ClusterID;		
		//it=find (ClusterList.begin(), ClusterList.end(), ClusterID);
			 ClusterList.push_back(ClusterID);
			 Egam.push_back(E[ipix]);
			 Tgam.push_back(T[ipix]);
			 Vgam.push_back(ArrayList[icenter]);
			 //Ngam++;		
	 }	
	 return;
}


Double_t TreeTrack::GetTheta(Int_t detID){
	return ArrayList[detID].Theta();
}

Double_t TreeTrack::GetPhi(Int_t detID){
	return ArrayList[detID].Phi();
}


void TreeTrack::ProcessTreeTest(Int_t nevlim, Int_t iev0){
	Int_t nev=tree->GetEntries();
	if(nev<1){cout<<"No events in tree or not open\n";return;}
	Int_t nlist=200;
	if(Ndets<nlist)nlist=Ndets;
	cout<<"First "<<nlist<<" detectors:"<<endl;
	cout<<"i th[i] ph[i] (degrees)\n";
	Double_t deg=180./acos(-1);
	for(int i=0;i<nlist;i++){
		cout<<i<<" "<<ArrayList[i].Theta()*deg<<" "<<ArrayList[i].Phi()*deg<<endl;
	}
	cout<<nev<<" events in input tree"<<endl;
	cout<<"Ng (Eg,Tg,Thg,Phg) Ex kR"<<endl;
	 
	for (Int_t iev=iev0;iev<(iev0+nevlim);iev++){
		tree->GetEvent(iev);
		TrackGammas(9);
		OutputEvent();//fills tree here
		ListEvent();
		}
	tracked->Write();
	cout<<Ntracked<<" tracked events"<<endl;
	CloseFiles();
	return;
}

void TreeTrack::ListEvent(){

	cout<<Ng<<" ";
	for (int i=0; i<Ng; i++){
		cout<<"("<<Eg[i]<<","<< Tg[i]<<","<< Thg[i]<<","<<Phg[i]<<") ";
	}
	cout<<Ex<<" "<<kR<<endl;
	cout<<"------------------------------------------"<<endl;
}

void TreeTrack::OpenTrackedFile(const char* namefile) {

	OutFile =  TFile::Open(namefile,"RECREATE");
	if(OutFile->IsOpen())cout<<"Output File  "<<namefile <<" is open"<<endl;
	
	tracked= new TTree("tracked","tracked");
	
	tracked->Branch("Ng",&Ng,"Ng/I");
	tracked->Branch("Eg",Eg,"Eg[Ng]/D");
	tracked->Branch("Tg",Tg,"Tg[Ng]/D");
	tracked->Branch("Thg",Thg,"Thg[Ng]/D");
	tracked->Branch("Phg",Phg,"Phg[Ng]/D");
	tracked->Branch("Ex",&Ex,"Ex/D");	
	tracked->Branch("kR",&kR,"kR/I");	
	return;
}

void TreeTrack::OutputEvent(){
	Ng=Ngam;
	for(Int_t i=0;i<Ng;i++){
		Eg[i]=Egam[i];
		Tg[i]=Tgam[i];
		Thg[i]=Vgam[i].Theta();
		Phg[i]=Vgam[i].Phi();
	}
		Ex=Exin;
		kR=kRin;
		tracked->Fill();
		Ntracked++;
		return;
}

void TreeTrack::ProcessTree(Int_t TrackType){
	Int_t nev=tree->GetEntries();
	if(nev<1){cout<<"No events in tree or not open\n";return;}
	cout<<nev<<" events in input tree"<<endl;
	for (Int_t iev=0;iev<nev;iev++){
		tree->GetEvent(iev);
		TrackGammas(TrackType);
		OutputEvent();
	}
	tracked->Write();
	cout<<Ntracked<<" tracked events"<<endl;
	CloseFiles();
	return;
}

void TreeTrack::TrackGammas_9pix(){

	Int_t npix=AddPixelEnergies();// number of fired pixels // totl E per pix is added
	Ngam=0;
	Egam.clear();
	Tgam.clear();
	Vgam.clear();		
	
	vector <Int_t> ClusterList;
	vector <Int_t>::iterator it;
	for (Int_t ipix=0;ipix<npix;ipix++){ // loop on fired pixels
		Int_t ClusterID=DetID[ipix]/9;// packs of 9 dets
		Int_t icenter=9*ClusterID;
		it=find (ClusterList.begin(), ClusterList.end(), ClusterID);
		if (it == ClusterList.end()){
			 ClusterList.push_back(ClusterID);
			 Egam.push_back(E[ipix]);
			 Tgam.push_back(T[ipix]);
			 Vgam.push_back(ArrayList[icenter]);
			 Ngam++;
		}
		else {
			 int d=distance(ClusterList.begin(),it);
			 Egam[d]+=E[ipix];		 
			 if(T[ipix]<Tgam[d]) Tgam[d]=T[ipix];
		}
	 }	
	 return;
}



Int_t TreeTrack::AddPixelEnergies(){	
	ClearVecs();
	vector <Int_t>::iterator it;
	Int_t npix=0;
	
	for (Int_t ihit=0;ihit<Nin;ihit++){ // loop on hits
		if(HitFilter(ihit)){;//applies hit filter (e.g. particle type, time...)
			it = find (DetID.begin(), DetID.end(), Din[ihit]);
			if (it == DetID.end()){
				DetID.push_back(Din[ihit]);
				E.push_back(Ein[ihit]);
				T.push_back(Tin[ihit]);
				npix++;
			}
			else {
				int d=distance(DetID.begin(),it);
				E[d]+=Ein[ihit];		 
				if(Tin[ihit]<T[d]) T[d]=Tin[ihit];
			}
		}
	 }
	 BlurrData(npix);
	 return npix;
}

void TreeTrack::BlurrData(Int_t npix){
	for(Int_t ipix=0;ipix<npix;ipix++){
		Double_t sigmaE=SigmaRes(E[ipix],DetType);
		Double_t sigmaT=TimeResolution;
		E[ipix]+=gRandom->Gaus(0.,sigmaE);
		T[ipix]+=gRandom->Gaus(0.,sigmaT);
	}
	return;
}

void TreeTrack::ClearVecs(){//clear pixel data
	E.clear();
	T.clear();
	DetID.clear();
	return;
}

bool  TreeTrack::HitFilter(Int_t ihit){
	
	if(!((Pin[ihit]==22)||(Pin[ihit]=2112))) return false;// only g or n
	if(Tin[ihit]>Tfilter) return false;  
	if(Ein[ihit]<Ethresh) return false;
	// skips to next ihit if none above, else stays the same
	return true;
}

Int_t TreeTrack::DetectorPosTable(const char* nametab){
	string line;
	Double_t x,y,z;
	Double_t th,ph,psi;
	TVector3 vpos;
	Int_t n,ndets=0;
	ifstream myfile (nametab);
	if (myfile.is_open()){
		cout<<"Detector position file "<<nametab<<" is open"<<endl;
		myfile>>n;
		while ( 1 ){
			myfile>>x>>y>>z>>ph>>th>>psi;
			if(!myfile.good())break;
			vpos.SetXYZ(x,y,z);
			ArrayList.push_back(vpos);
			ndets++;
		}
		cout<<"last detector read:"<<endl;
			Int_t last=ndets-1;
			cout<<"x,y,z,ph,th\n";
			cout<<ArrayList[last].X()<<","<< ArrayList[last].Y()<<","<< ArrayList[last].Z()<<","<<
				ArrayList[last].Phi()<<","<<ArrayList[last].Theta()<<endl;
			cout<<"\n";
		myfile.close();
	}
	if(n!=ndets)cout<<"error reading detector positions file"<<endl;
	cout<<ndets<<" detector positions read from "<< nametab<<endl;
	return ndets;
}


void TreeTrack::OpenTree(const char* namefile) {
		InFile=TFile::Open(namefile);
		if(InFile->IsOpen())cout<<"File  "<<namefile <<" is open"<<endl;
		
		tree = (TTree*) InFile->Get("tree");	
		tree->SetBranchAddress("Nhit",&Nin);
		tree->SetBranchAddress("DetIDhit",Din);
		tree->SetBranchAddress("Ehit",Ein);
		tree->SetBranchAddress("Thit",Tin);
		tree->SetBranchAddress("PDGchit",Pin);	
		tree->SetBranchAddress("ExDCE",&Exin);	
		tree->SetBranchAddress("NRhit",&NRin);	
		tree->SetBranchAddress("kR",&kRin);	
		return;

}

void TreeTrack::CloseFiles(){
	InFile->Close();
	OutFile->Close();
}

