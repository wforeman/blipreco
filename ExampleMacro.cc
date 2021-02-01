//////////////////////////////////////////////////////////////////
// 
//  Must have ROOT data file containing analysistree/anatree(_blip)
//  in same directory.
//
//  To run:
//  > root -l ExampleMacro.cc
//
////////////////////////////////////////////////////////////////// 

// Input and output file names
std::string   fFileName       = "AnaTree_BlipReco_positron_1MeV.root";
std::string   fOutFileName    = "output.root";

// Declare function to be defined later
float     Mass(int);

// Declare TFile/TTree objects, defined as pointers (*) for now
TFile*    fOutFile;
TTree*    fTree;

// Declare ROOT histograms, defined as pointers (*) for now
TH1D*     h_BlipEnergy;
TH1D*     h_SumBlipEnergy;

// AnaTree variables (I append these with "_" to signify
// that they are special tree variables, to distinguish them 
// from any local variables I might create in my code)
const int kMax = 10000;
double    _NuEnergy;
int       _event;
int       _no_hits_stored;
int       _geant_list_size;
int       _inTPCActive[kMax];
int			  _pdg[kMax];
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
int       _blip_list_size;
int       _blipTrackId[kMax];
float     _blipEnergy[kMax];
float     _blipLocX[kMax];
float     _blipLocY[kMax];
float     _blipLocZ[kMax]; 


//=====================================================
void ExampleMacro(){
  
  std::cout<<"Running macro...\n";

  // Open the file and get the TTree
  std::cout<<"Reading input file\n";
  TFile* file = new TFile(fFileName.c_str(),"READ");
  fTree = (TTree*)file->Get("anatree_blip");
  
  // Make output file to store plots
  fOutFile = new TFile(fOutFileName.c_str(),"recreate");
  fOutFile->cd();
  
  // Set up the branch structure (each branch should
  // correspond to a variable already defined above).
  std::cout<<"Defining branches\n";
  TBranch *br = 0;
  fTree->SetBranchAddress("event",&_event);
  fTree->SetBranchAddress("geant_list_size",&_geant_list_size);
  fTree->SetBranchAddress("inTPCActive",&_inTPCActive);
  fTree->SetBranchAddress("pdg",&_pdg);
  fTree->SetBranchAddress("Eng",&_Eng);
  fTree->SetBranchAddress("EndE",&_EndE);
  fTree->SetBranchAddress("Px",&_Px);
  fTree->SetBranchAddress("Py",&_Py);
  fTree->SetBranchAddress("Pz",&_Pz);
  fTree->SetBranchAddress("StartT",&_StartT);
  fTree->SetBranchAddress("EndT",&_EndT);
  fTree->SetBranchAddress("StartPointx",&_StartPointx);
  fTree->SetBranchAddress("StartPointy",&_StartPointy);
  fTree->SetBranchAddress("StartPointz",&_StartPointz);
  fTree->SetBranchAddress("EndPointx",&_EndPointx);
  fTree->SetBranchAddress("EndPointy",&_EndPointy);
  fTree->SetBranchAddress("EndPointz",&_EndPointz);
  fTree->SetBranchAddress("pathlen",&_pathlen);
  fTree->SetBranchAddress("NumberDaughters",&_NumberDaughters);
  fTree->SetBranchAddress("Mother",&_Mother);
  fTree->SetBranchAddress("TrackId",&_TrackId);
  fTree->SetBranchAddress("processname",&_processname,&br);
  fTree->SetBranchAddress("blip_list_size",&_blip_list_size);
  fTree->SetBranchAddress("BlipTrackId",&_blipTrackId);
  fTree->SetBranchAddress("BlipEnergy",&_blipEnergy);
  fTree->SetBranchAddress("BlipLocX",&_blipLocX);
  fTree->SetBranchAddress("BlipLocY",&_blipLocY);
  fTree->SetBranchAddress("BlipLocZ",&_blipLocZ);

  // Create the histograms that we'll fill later.
  h_BlipEnergy            = new TH1D("BlipEnergy",    "Individual Blip Energy;Energy (MeV)",100,0,10.);
  h_SumBlipEnergy         = new TH1D("SumBlipEnergy", "Summed Blip Energy;Energy (MeV)",100,0,10.);

  // --------------------------------------------
  // Loop over the events in the tree
  for(int iEvent=0; iEvent<fTree->GetEntries(); iEvent++){
    
    // The GetEntry call retrieves all the information
    // for this particular event and saves it into the
    // variables we defined above.
    fTree->GetEntry(iEvent);
    
    std::cout
    <<"--------------------------------------------------------------------------\n"
    <<"Event "<<iEvent<<" has "<<_geant_list_size<<" particles and "<<_blip_list_size<<" visible energy deposits\n";


    // ---------------------------------------
    // Loop over particles
    for(int i=0; i<_geant_list_size; i++){
      
      // do some calculations on the data from this particle, or
      // save data into simpler local variables.
      TVector3 loc(_StartPointx[i],_StartPointy[i],_StartPointz[i]);
      TVector3 locEnd(_EndPointx[i],_EndPointy[i],_EndPointz[i]);
      std::string proc = _processname->at(i);
      int   pdg     = _pdg[i];
      int   trackId = _TrackId[i]; 
      float E       = 1e3*_Eng[i]; // Converts GeV --> MeV
      float endE    = 1e3*_EndE[i]; // Converts GeV --> MeV
      float mass    = Mass(i);
      float KE      = E - mass;
      float KEf     = endE - mass;
      int   mother  = _Mother[i];
      int   nD      = _NumberDaughters[i];
      float startT  = _StartT[i];
      float endT    = _EndT[i];
      float dL      = _pathlen[i];
      
      // calculates the visible energy deposited
      //float edep    = CalcEnergyDep(i);
      
      // print out each particle
      printf("  %3i PDG: %10i,  dL=%8.3f, KE0=%8.3f,  KEf=%8.3f, T0=%8.3f, Tf=%8.3f, moth=%3i, %12s, Ndaught=%i\n",
        trackId,
        pdg,
        dL,
        KE,
        KEf,
        //edep,
        startT,
        endT,
        mother,
        proc.c_str(),
        nD
        );
    }//>> end PARTICLE loop


    // Energy depositions in the event (calculated from "CalcEnergyDep()" above)
    // eventually appear as isolated hits in the reconstruction, which we call
    // blips. This blip info has already been saved into the TTree.
    float total_blipEnergy  = 0.;
    
    std::cout<<"Visible energy depositions created in Event "<<iEvent<<":\n";
    for(int j=0; j<_blip_list_size; j++){
      
      std::cout<<"  blip "<<j<<"-->  E= "<<_blipEnergy[j]<<" MeV, initiated by particle with track ID "<<_blipTrackId[j]<<"\n";
      
      // fill histo
      h_BlipEnergy->Fill(_blipEnergy[j]);

      // add to total blip energy
      total_blipEnergy += _blipEnergy[j];

    }//>> end BLIP loop


    // Fill total energy histogram
    h_SumBlipEnergy->Fill(total_blipEnergy);
  
  }//>> end EVENT loop

  // Save histograms to output and close file
  h_BlipEnergy      ->Write();
  h_SumBlipEnergy   ->Write();
  fOutFile          ->Close();
  return;
}


//=============================================================
// Calculates mass of particle 'i'
float Mass( int i ) {
  if( _pdg[i] == 22 ) return 0.;
  if( abs(_pdg[i])==12 || abs(_pdg[i])==14 || abs(_pdg[i])==16 ) return 0.;
  float p = 1e3*sqrt( pow(_Px[i],2) + pow(_Py[i],2) + pow(_Pz[i],2) );
  float E = 1e3*_Eng[i];
  return sqrt( pow(E,2) - pow(p,2) );
}
