//////////////////////////////////////////////////////////////////
// 
//  Must have ROOT data file containing analysistree/anatree in 
//  same directory.
//
//  To run:
//  > root -l BlipReco_AnaMacro.cc
//
////////////////////////////////////////////////////////////////// 
#include "blip.h"
#include "tools.h"

// ============================================================
// ===================   Parameters ===========================
#define nconfigs_SphereR 3
//std::string             fFileName     = "../mcfiles/anatree_gamma1_larworld.root";
//std::string             fFileName    = "../mcfiles/anatree_gamma2_larworld.root";
std::string             fFileName    = "../mcfiles/anatree_neutrons_1eV.root";
//std::string             fFileName    = "../mcfiles/anatree_neutron_0to20MeV_larworld.root";
std::string             fTreeName     = "analysistree/anatree";
std::string             fOutFileName  = "plots.root";
float                   fThreshold    = 0.0750; // 75 keV
float                   fSmear        = 0.00; // 50 keV
float                   fMinSep       = 0.2; // cm
float                   fSphereR      = 60.;

//===============================================================

void NeutronReco_AnaMacro();
void configure();
void reco();
void loopTheBlips(std::vector<EnergyDeposit> blips );
void makePlots();

float     fNeutronEnergy;
float     fRecoEnergy;
TVector3  fStartPoint;

TFile*    fOutFile;
TTree*    fTree;
TH2D*     h_Energy_TrueVsReco;

void configure(){
  
  // open the file and set up the TTree
  TFile* file = new TFile(fFileName.c_str(),"READ");
  fTree = (TTree*)file->Get(fTreeName.c_str());
  setBranches(fTree);
  
  // make output file to store plots
  fOutFile = new TFile(fOutFileName.c_str(),"recreate");
  fOutFile->cd();

  // make histograms
  h_Energy_TrueVsReco   = new TH2D("Energy_TrueVsReco",
                                   "R=60cm;True Neutron Energy (MeV);Reconstructed Energy (MeV)",
                                  100,0,20,100,0,20);
  h_Energy_TrueVsReco->SetOption("COLZ");

}


// =============================================================
void NeutronReco_AnaMacro(){
  
  // configure histograms and TFile
  configure();
  reco();
  makePlots();

  // write histos
  h_Energy_TrueVsReco->Write();
  fOutFile->Close();
  return;
}



// ============================================================
void reco(){
  
  // --------------------------------------------
  // loop over the events
  for(int iEvent=0; iEvent<fTree->GetEntries(); iEvent++){
    fTree->GetEntry(iEvent);
  
    // Make a vector of blip objects to fill,
    std::vector<EnergyDeposit> v_blips;
    float total_edep = 0.;
    
    // Loop over particles and fill in the blips.
    int nParticles = _geant_list_size;
    for(int i=0; i<nParticles; i++){

      TVector3 loc(_StartPointx[i],_StartPointy[i],_StartPointz[i]);
      TVector3 locEnd(_EndPointx[i],_EndPointy[i],_EndPointz[i]);
      std::string proc = _processname->at(i);
      int   pdg     = _pdg[i];
      int   trackId = _TrackId[i]; 
      float KE      = 1e3*(_Eng[i]-_Mass[i]);
      float KEf     = 1e3*(_EndE[i]-_Mass[i]);
      int   mother  = _Mother[i];
      int   nD      = _NumberDaughters[i];
      float startT  = _StartT[i];
      float endT    = _EndT[i];
      float dL      = _pathlen[i];

      float edep = CalcEnergyDep(i);
      total_edep += edep;
      
      if( proc == "primary" ) {
        fNeutronEnergy  = KE;
        fStartPoint     = loc;
      }

      // is it an electron/positron?
      if(     fabs(pdg) == 11 
          ||  (pdg==2212 && KE > 3.) ){
        EnergyDeposit b;
        b.Location = loc;
        b.Energy = edep;
        b.PathLength = dL;
        v_blips.push_back(b);
      }
      
      // debugging output
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
        nD
        );
      */
    
    }//>> end particle loop

    // .................................
    // merging, thresholding, and smearing
    MergeBlips    (v_blips,fMinSep    );
    ThresholdBlips(v_blips,fThreshold );
    SmearBlips    (v_blips,fSmear     );

    printf("Event %5d -- smear: %6.2f keV, blips found: %d\n",iEvent,fSmear*1e3,(int)v_blips.size());
  
    fRecoEnergy = 0.;
    for(size_t i=0; i<v_blips.size(); i++){
      EnergyDeposit thisBlip = v_blips.at(i);
      float d = (thisBlip.Location - fStartPoint).Mag();
      if( d < fSphereR || fSphereR < 0 ) {
        fRecoEnergy += thisBlip.Energy;
      }
    }

    //if( fRecoEnergy > (fNeutronEnergy) ) std::cout<<"WTF\n";

    // fill histograms
    if( fNeutronEnergy > 0.2 ) 
    h_Energy_TrueVsReco->Fill(fNeutronEnergy,fRecoEnergy);

  }

}

// ============================================================
void makePlots(){

  TCanvas* c0 = new TCanvas("c0","c0",50,50);

  TCanvas* c1 = new TCanvas("c1","c1",700,450);
  gPad->SetLogz(1);
  gStyle->SetOptStat(0);
  FormatAxes(h_Energy_TrueVsReco,0.045,0.045,1,1);
  h_Energy_TrueVsReco->DrawCopy();
  c0->cd();
  
  TCanvas* c2 = new TCanvas("c2","c2",700,450);
  TH1D* h;
  //find bin for 10 MeV
  int bin = h_Energy_TrueVsReco->GetXaxis()->FindBin(10.);
  h = h_Energy_TrueVsReco->ProjectionY("yproj_10MeV",bin-1,bin);
  h->DrawCopy();
  c0->cd();

  
  TF1 f_gaus("f_gaus","gaus(0)");
  
  TCanvas* c_slice;
  for(size_t i=3; i<18; i++){
    float E = float(i);
    int bin = h_Energy_TrueVsReco->GetXaxis()->FindBin(E);
    c_slice = new TCanvas(Form("slice_%3.1fMeV",E),Form("slice_%3.1fMeV",E),500,500);
    h = h_Energy_TrueVsReco->ProjectionY(Form("yproj_%3.1fMeV",E),bin-1,bin);
    h ->DrawCopy();
    f_gaus.SetRange(0.5,20);
    h->Fit("f_gaus","RQ");
    std::cout<<"Slice "<<E<<" MeV --> overall gaus fit mu: "<<f_gaus.GetParameter(1)<<"  sig: "<<f_gaus.GetParameter(1)<<"\n";
    c0->cd();
  }


}


