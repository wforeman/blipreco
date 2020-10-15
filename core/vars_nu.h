
// AnaTree variables
const int kMax = 10000;
double    _NuEnergy;
int       _event;
int       _no_hits_stored;
int       _geant_list_size;
int       _inTPCActive[kMax];
int			  _pdg[kMax];
//double		_Mass[kMax];
double     _Eng[kMax];
double     _EndE[kMax];
double     _Px[kMax];
double     _Py[kMax];
double     _Pz[kMax];
double     _StartT[kMax];
double     _EndT[kMax];
double     _StartPointx[kMax];
double     _StartPointy[kMax];
double     _StartPointz[kMax];
double     _EndPointx[kMax];
double     _EndPointy[kMax];
double     _EndPointz[kMax];
double     _pathlen[kMax];
int       _NumberDaughters[kMax];
int       _Mother[kMax];
int       _TrackId[kMax];
int       _process_primary[kMax];
std::vector<std::string> *_processname = 0;

void setBranches(TTree *tree){
  TBranch *br = 0;
  tree->SetBranchAddress("enu_truth",&_NuEnergy);
  tree->SetBranchAddress("event",&_event);
  tree->SetBranchAddress("geant_list_size",&_geant_list_size);
  tree->SetBranchAddress("pdg",&_pdg);
  tree->SetBranchAddress("Eng",&_Eng);
  tree->SetBranchAddress("Px",&_Px);
  tree->SetBranchAddress("Py",&_Py);
  tree->SetBranchAddress("Pz",&_Pz);
  //tree->SetBranchAddress("EndE",&_EndE);
  //tree->SetBranchAddress("StartT",&_StartT);
  //tree->SetBranchAddress("EndT",&_EndT);
  tree->SetBranchAddress("StartPointx",&_StartPointx);
  tree->SetBranchAddress("StartPointy",&_StartPointy);
  tree->SetBranchAddress("StartPointz",&_StartPointz);
  tree->SetBranchAddress("EndPointx",&_EndPointx);
  tree->SetBranchAddress("EndPointy",&_EndPointy);
  tree->SetBranchAddress("EndPointz",&_EndPointz);
  //tree->SetBranchAddress("pathlen",&_pathlen);
  tree->SetBranchAddress("NumberDaughters",&_NumberDaughters);
  tree->SetBranchAddress("Mother",&_Mother);
  tree->SetBranchAddress("TrackId",&_TrackId);
  tree->SetBranchAddress("G4Process",&_processname,&br);
  //tree->SetBranchAddress("processname",&_processname,&br);
}
