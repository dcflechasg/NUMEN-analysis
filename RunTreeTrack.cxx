#include "TreeTrack.cxx"

void RunTreeTrack(Int_t CrystalTypeID=1, Int_t TrackType=1, const char* infile="DCEk1.0.root", 
                  const char* outfile="Tracked.root", const char * dettable="lplace-radialC.dat" ){
	gSystem->Load("libPhysics.so");
	TreeTrack TrTr;
	TrTr.DetType=CrystalTypeID;// 1 LYSO, 2 LaBr3
	TrTr.OpenTree(infile);
	TrTr.OpenTrackedFile(outfile);
	TrTr.DetectorPosTable(dettable);
	TrTr.ProcessTree(TrackType);
	return;
}
