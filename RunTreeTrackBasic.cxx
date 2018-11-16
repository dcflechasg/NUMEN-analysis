#include "TreeTrack.cxx"

void RunTreeTrack(Int_t CrystalTypeID=1, Int_t TrackType=9, const char* infile="TreeRad2MeV.root", 
                  const char* outfile="TrackRad2MeV.root", const char * dettable="lplace-radial-1set.dat" ){
	gSystem->Load("libPhysics.so");
	TreeTrack TrTr;
	TrTr.DetType=CrystalTypeID;// 1 LYSO, 2 LaBr3
	TrTr.OpenTree(infile);
	TrTr.OpenTrackedFile(outfile);
	TrTr.DetectorPosTable(dettable);
	TrTr.ProcessTree(TrackType);
	return;
}
