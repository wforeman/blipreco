//////////////////////////////////////////////////////////////////
// 
//  BlipRecoExample
//
//  This is a simple example of how to reconstruct "blips" using
//  the Geant4-level truth information stored in the AnaTree.
//
//  Must have ROOT data file containing analysistree/anatree in 
//  same directory.
//
//  To run:
//  > root -l BlipRecoExample.cc
//
////////////////////////////////////////////////////////////////// 
#include "core/vars.h"
#include "core/blip.h"
#include "core/tools.h"

// ===================   Parameters ===========================
std::string             fFileName     = "../../mcfiles/anatree_neutron_10MeV_larworld.root";
std::string             fTreeName     = "analysistree/anatree";
std::string             fOutFileName  = "plots.root";

float                   fThreshold    = 0.0750; //> Blip threshold (75 keV default)
float                   fSmear        = 0.00;   //> Smearing (50 keV usually)
float                   fMinSep       = 0.20;   //> Min blip separation used 
                                                //  during blip merging stage

//===============================================================
void BlipRecoExample();
void configure();

TFile*    fOutFile;
TTree*    fTree;
TH1D*     h_BlipEnergy;
TH1D*     h_BlipMult;

//================================================================
void configure(){
  // make the histograms
  h_BlipEnergy      = new TH1D("BlipEnergy","Individual Blip Energies;Individual Blip Energy (MeV);Entries",1000,0.0,10.0);
  h_BlipMult        = new TH1D("BlipMult","Blip Multiplicity;Number of Blips",100,0,100);
  // open the file and set up the TTree
  TFile* file = new TFile(fFileName.c_str(),"READ");
  fTree = (TTree*)file->Get(fTreeName.c_str());
  setBranches(fTree);
}


// =============================================================
void BlipRecoExample(){

  // configure histograms and TFile
  configure();
  
  // --------------------------------------------
  // loop over the events
  for(int iEvent=0; iEvent<fTree->GetEntries(); iEvent++){
    fTree->GetEntry(iEvent);
    
    // Make a vector of Blip objects to fill
    std::vector<EnergyDeposit> v_blips;

    // Loop over particles and fill in the blips.
    int nParticles = _geant_list_size;
    for(int i=0; i<nParticles; i++){

      // ignore any particles outside TPC active
      if( _inTPCActive[i] == 0 ) continue;
    
      TVector3 loc(_StartPointx[i],_StartPointy[i],_StartPointz[i]);
      TVector3 locEnd(_EndPointx[i],_EndPointy[i],_EndPointz[i]);
      int pdg       = _pdg[i];
      int trackId   = _TrackId[i]; 
      float KE      = 1e3*(_Eng[i]-_Mass[i]);
      float KEf     = 1e3*(_EndE[i]-_Mass[i]);
      std::string proc = _processname->at(i);
      int mother    = _Mother[i];
      int nD        = _NumberDaughters[i];
      float startT  = _StartT[i];
      float endT    = _EndT[i];
      float dL      = _pathlen[i];
     
      // calc energy dep 
      float edep = CalcEnergyDep(i);

      // if electron, make blip
      if( fabs(pdg) == 11 ) {
        EnergyDeposit b;
        b.Energy    = edep;
        b.Location  = loc;
        b.isGrouped = false;
        v_blips.push_back(b);
      }
      
      // (for debugging output -- keep commented out during normal running)
      /*
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
          nD);
      */

    }//>> end particle loop
    
    // .................................
    // merge close-together blips
    MergeBlips(v_blips,fMinSep);
   
    // fill energy histogram 
    for(int i=0; i<v_blips.size(); i++){
      h_BlipEnergy->Fill(v_blips.at(i).Energy);
    }

    // .................................
    // do thresholding
    ThresholdBlips(v_blips,fThreshold);
    h_BlipMult->Fill(v_blips.size());

    // ................................
    // do smearing
    SmearBlips(v_blips,fSmear);

    printf("Event %5d -- smear: %6.2f keV, blips found: %d\n",iEvent,fSmear*1e3,(int)v_blips.size());

  }//>> end loop over events

  
  // ------------------------------------------------- 
  // make output file to store plots
  fOutFile = new TFile(fOutFileName.c_str(),"recreate");
  fOutFile->cd();
  
  // save all the lower-level histograms
  h_BlipEnergy->Write();
  h_BlipMult->Write();
  fOutFile->Close();
  return;
}

