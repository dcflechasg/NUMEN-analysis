#include "TVector3.h"

//General functions

Double_t Hip(Double_t Cat1, Double_t Cat2){
	return sqrt( pow(Cat1,2)+pow(Cat2,2) );
	
}

Double_t Ang(Double_t th1, Double_t ph1, Double_t th2, Double_t ph2){
	TVector3 v1(sin(th1)*cos(ph1),sin(th1)*sin(ph1),cos(th1));
	TVector3 v2(sin(th2)*cos(ph2),sin(th2)*sin(ph2),cos(th2));
	return acos(v1.Dot(v2));
}

	
Double_t RelFWHM(Double_t EG, Int_t IDtype){ // gamma det FWHM %, E in MeV
	Double_t res=1.;
	Double_t a,b,c,d,y;
	switch (IDtype){
		case 1:// LYSO intrinsic resolution parameters * 1.2
			a=5.04*1.2;
			b=-0.48*1.2;
			y=a+b*log(EG*1000.);
			res=exp(y);
			break;
		case 2://LaBr3 = LYSO intrinsic resolution parameters * 1.2/4
			a=5.04*1.2;
			b=-0.48*1.2;
			y=a+b*log(EG*1000.);
			res=exp(y)/4.;
			break;
	}
	return res;
}

Double_t SigmaRes(Double_t EG, Int_t IDtype){
	return RelFWHM(EG,IDtype)*EG/235.;
}

Double_t DopplerCorrection(Double_t betarec, Double_t thetarec, Double_t phirec, Double_t thetagam, Double_t phigam){
	gSystem->Load("libPhysics.so");
	TVector3 Vr;
	Vr.SetMagThetaPhi(1.,thetarec,phirec);
	TVector3 Vg;
	Vg.SetMagThetaPhi(1.,thetagam,phigam);
	Double_t fDopp=(1.+betarec*Vr.Dot(Vg));
	return fDopp;
}

Double_t Beta(Double_t Arec, Double_t Erec){
	Double_t UMA=931.5;//MeV
	//Double_t beta=sqrt(2.*Erec/(Arec*UMA);// non relativistic
	Double_t Kred=Erec/(Arec*UMA);
	Double_t beta=sqrt(Kred*(2.+Kred))/(1+Kred);// relativistic
	return beta;
}
