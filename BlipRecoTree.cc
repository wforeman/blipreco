//////////////////////////////////////////////////////////////////
// 
//  Basic blip reco alg. Does the standard blip reconstruction
//  for all depositions in the event and saves it into a new 
//  TTree alongside the G4 truth information from the original
//  anatree file.
//
//  Must have ROOT data file containing analysistree/anatree in 
//  same directory.
//
//  To run:
//  > root -l [macro-file]
//
////////////////////////////////////////////////////////////////// 
#include "core/vars.h"
#include "core/blip.h"
#include "core/tools.h"

// ===================   Parameters ===========================
std::string             fFileName       = "../mcfiles/AnaTree_elMinus_1MeV_10kEvents.root";
std::string             fTreeName       = "analysistree/anatree";
std::string             fOutFileName    = "plots.root";
float                   fThreshold      = 0.075; // 75 keV
float                   fSmear          = 0.0; // 50 keV
float                   fMinSep         = 0.2;   // 0.2 cm

//===============================================================

void      BlipRecoTree();
void      configure();
void      reco();

TFile*    fOutFile;
TTree*    fTree;

TH1D*     h_PrimaryEnergy;
TH1D*     h_NumBlips;
TH1D*     h_BlipEnergy;
TH1D*     h_SumBlipEnergy;
TH1D*     h_TotalDepEnergy;

TTree*    fTree2; // new tree w/blip info saved
int       _numBlips;
int       _blipTrackId[kMax];
float     _blipEnergy[kMax];
float     _blipLocX[kMax];
float     _blipLocY[kMax];
float     _blipLocZ[kMax]; 

// ============================================================
void configure(){

  // open the file and set up the TTree
  TFile* file = new TFile(fFileName.c_str(),"READ");
  fTree = (TTree*)file->Get(fTreeName.c_str());
  setBranches(fTree);
  
  // make output file to store plots
  fOutFile = new TFile(fOutFileName.c_str(),"recreate");
  fOutFile->cd();

  // copy tree to tree2 and add some new branches
  //fTree2  = new TTree();
  fTree2  = fTree->CloneTree();
  fTree2  ->SetName("anatree_blip");
  fTree2  ->Reset();
  fTree2  ->Branch("blip_list_size",&_numBlips,"blip_list_size/I");
  fTree2  ->Branch("BlipTrackId",_blipTrackId,"BlipTrackId[blip_list_size]/F");
  fTree2  ->Branch("BlipEnergy",_blipEnergy,"BlipEnergy[blip_list_size]/F");
  fTree2  ->Branch("BlipLocX",_blipLocX,"BlipLocX[blip_list_size]/F");
  fTree2  ->Branch("BlipLocY",_blipLocY,"BlipLocY[blip_list_size]/F");
  fTree2  ->Branch("BlipLocZ",_blipLocZ,"BlipLocZ[blip_list_size]/F");
  
  // make the histograms
  float maxE              = 10.;
  h_PrimaryEnergy         = new TH1D("PrimaryEnergy", "Primary Energy;Energy (MeV)",100,0,maxE);
  h_NumBlips              = new TH1D("NumBlips",      "Number of Blips",  20,0,20);
  h_BlipEnergy            = new TH1D("BlipEnergy",    "Individual Blip Energy;Energy (MeV)",100,0,maxE);
  h_SumBlipEnergy         = new TH1D("SumBlipEnergy", "Summed Blip Energy;Energy (MeV)",100,0,maxE);
  h_TotalDepEnergy        = new TH1D("TotalDepEnergy","Total Energy Deposited;Deposited Energy (MeV)",100,0,maxE);
  
}


// =============================================================
void BlipRecoTree(){

  configure();
  reco();
 
  /*  
  h_PrimaryEnergy ->Write();
  h_NumBlips      ->Write();
  h_BlipEnergy    ->Write();
  h_SumBlipEnergy ->Write();
  h_TotalDepEnergy->Write();
  */

  fTree2->Write();

  fOutFile->Close();
  return;
}



// ============================================================
void reco(){
  
  // --------------------------------------------
  // loop over the events
  for(int iEvent=0; iEvent<fTree->GetEntries(); iEvent++){
    fTree->GetEntry(iEvent);

    // reset blip variables 
    std::fill( std::begin(_blipEnergy),   std::end(_blipEnergy),  0);
    std::fill( std::begin(_blipTrackId),  std::end(_blipTrackId), 0);
    std::fill( std::begin(_blipLocX),     std::end(_blipLocX),    0);
    std::fill( std::begin(_blipLocY),     std::end(_blipLocY),    0);
    std::fill( std::begin(_blipLocZ),     std::end(_blipLocZ),    0);
    _numBlips               = 0;
    float total_edep        = 0.;
    float total_blipEnergy  = 0.;
    float total_primEnergy  = 0.;
    
    // Make a vector of EnergyDeposit objects to fill,
    std::vector<EnergyDeposit> v_blips;
    
    
    // Loop over particles and fill in the blips.
    int nParticles = _geant_list_size;
    for(int i=0; i<nParticles; i++){

      TVector3 loc(_StartPointx[i],_StartPointy[i],_StartPointz[i]);
      TVector3 locEnd(_EndPointx[i],_EndPointy[i],_EndPointz[i]);
      std::string proc = _processname->at(i);
      int   pdg     = _pdg[i];
      int   trackId = _TrackId[i]; 
      float p       = 1e3*sqrt( pow(_Px[i],2) + pow(_Py[i],2) + pow(_Pz[i],2) );
      float E       = 1e3*_Eng[i];
      float endE    = 1e3*_EndE[i];
      float mass    = sqrt( pow(E,2) - pow(p,2) );
      float KE      = E - mass;
      float KEf     = endE - mass;
      int   mother  = _Mother[i];
      int   nD      = _NumberDaughters[i];
      float startT  = _StartT[i];
      float endT    = _EndT[i];
      float dL      = _pathlen[i];

      // calculate visible energy from this particle
      float edep = CalcEnergyDep(i);
      total_edep += edep;
      
      if( proc == "primary" ) total_primEnergy += KE;

      // record blip
      if( edep ) AddNewBlip(v_blips,i);
//        EnergyDeposit b = MakeNewBlip(i);
//        v_blips.push_back(b);
//      }
      
      // debugging output
      printf("  %3i PDG: %10i,  dL=%8.3f, KE0=%8.3f,  KEf=%8.3f, Edep=%8.3f, T0=%8.3f, Tf=%8.3f, moth=%3i, %12s, Ndaught=%i\n",
        trackId,
        pdg,
        dL,
        KE,
        KEf,
        edep,
        startT,
        endT,
        mother,
        proc.c_str(),
        nD
        );
    
    }//>> end particle loop

    // .................................
    // merging, thresholding, and smearing
    MergeBlips    (v_blips,fMinSep    );
    ThresholdBlips(v_blips,fThreshold );
    SmearBlips    (v_blips,fSmear     );
    printf("Event %5d -- smear: %6.2f keV, blips found: %d\n",iEvent,fSmear*1e3,(int)v_blips.size());
    
    // calculate total blip energy
    _numBlips = v_blips.size();
    for(size_t i=0; i<v_blips.size(); i++){ 
      h_BlipEnergy->Fill(v_blips.at(i).Energy);
      total_blipEnergy += v_blips.at(i).Energy;
      _blipTrackId[i] = v_blips.at(i).TrackId;
      _blipEnergy[i]  = v_blips.at(i).Energy;
      _blipLocX[i]    = v_blips.at(i).Location.X();
      _blipLocY[i]    = v_blips.at(i).Location.Y();
      _blipLocZ[i]    = v_blips.at(i).Location.Z();
    }
 
    // fill histograms
    h_PrimaryEnergy     ->Fill(total_primEnergy);
    h_TotalDepEnergy    ->Fill(total_edep);
    h_SumBlipEnergy     ->Fill(total_blipEnergy);
    h_NumBlips          ->Fill(v_blips.size());
    
    fTree2->Fill();

  }//>> end loop over events

}

