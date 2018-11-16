#include "functions.cxx"

Double_t beta(){
	gSystem->Load("libPhysics.so");
	return Beta(20.,292.);
}

Double_t fDop(Double_t th){
	return (1.+beta()*cos(th));
}
	
	
