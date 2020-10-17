
// AnaTree variables
const int kMax = 10000;
double    _NuEnergy;
int       _event;
int       _no_hits_stored;
int       _geant_list_size;
int       _inTPCActive[kMax];
int			  _pdg[kMax];
//float     _Mass[kMax];
float     _Eng[kMax];
float     _EndE[kMax];
float     _Px[kMax];
float     _Py[kMax];
float     _Pz[kMax];
float     _StartT[kMax];
float     _EndT[kMax];
float     _StartPointx[kMax];
float     _StartPointy[kMax];
float     _StartPointz[kMax];
float     _EndPointx[kMax];
float     _EndPointy[kMax];
float     _EndPointz[kMax];
float     _pathlen[kMax];
int       _NumberDaughters[kMax];
int       _Mother[kMax];
int       _TrackId[kMax];
int       _process_primary[kMax];
std::vector<std::string> *_processname = 0;

void setBranches(TTree *tree){
  TBranch *br = 0;
  tree->SetBranchAddress("event",&_event);
  tree->SetBranchAddress("geant_list_size",&_geant_list_size);
  tree->SetBranchAddress("inTPCActive",&_inTPCActive);
  tree->SetBranchAddress("pdg",&_pdg);
//  tree->SetBranchAddress("Mass",&_Mass);
  tree->SetBranchAddress("Eng",&_Eng);
  tree->SetBranchAddress("EndE",&_EndE);
  tree->SetBranchAddress("Px",&_Px);
  tree->SetBranchAddress("Py",&_Py);
  tree->SetBranchAddress("Pz",&_Pz);
  tree->SetBranchAddress("StartT",&_StartT);
  tree->SetBranchAddress("EndT",&_EndT);
  tree->SetBranchAddress("StartPointx",&_StartPointx);
  tree->SetBranchAddress("StartPointy",&_StartPointy);
  tree->SetBranchAddress("StartPointz",&_StartPointz);
  tree->SetBranchAddress("EndPointx",&_EndPointx);
  tree->SetBranchAddress("EndPointy",&_EndPointy);
  tree->SetBranchAddress("EndPointz",&_EndPointz);
  tree->SetBranchAddress("pathlen",&_pathlen);
  tree->SetBranchAddress("NumberDaughters",&_NumberDaughters);
  tree->SetBranchAddress("Mother",&_Mother);
  tree->SetBranchAddress("TrackId",&_TrackId);
  tree->SetBranchAddress("processname",&_processname,&br);
}
