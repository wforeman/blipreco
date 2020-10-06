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
//std::string             fFileName    = "../mcfiles/anatree_gamma2_larworld.root";
std::string             fFileName    = "../mcfiles/anatree_neutrons_1eV.root";
//std::string             fFileName    = "../mcfiles/anatree_neutron_0to20MeV_larworld.root";
//std::string             fFileName    = "../mcfiles/anatree_neutron_10MeV_larworld.root";
std::string             fTreeName     = "analysistree/anatree";
std::string             fOutFileName  = "plots.root";
float                   fThreshold    = 0.100; // MeV (default 0.075)
float                   fSmear        = 0.00; // 50 keV
float                   fMinSep       = 0.2; // cm
std::vector<float>      fSphereR      = {30., 60.};
#define NSphereR 2

//===============================================================

void NeutronReco_AnaMacro();
void configure();
void reco();
void makePlots();

void BlipGrouping(std::vector<EnergyDeposit> blips );

TVector3  fStartPoint;
float     fNeutronEnergy;
float     fRecoEnergy[NSphereR];
int     fRecoMult[NSphereR];
float     fNCaptures;

TFile*    fOutFile;
TTree*    fTree;

TH1D*     h_Mult;
TH1D*     h_BlipE;
TH1D*     h_MaxBlipE;
TH1D*     h_AveBlipE;
TH1D*     h_TotBlipE;
TH1D*     h_BlipDist;
TH1D*     h_BlipTime;
TH1D*     h_NCaptPerEvt;
TH1D*     h_NCaptTime;
TH1D*     h_NCaptDist;
TH1D*     h_NCaptGammaE;
TH1D*     h_NCaptGammaESum;
TH1D*     h_NScatGammaE;
TH2D*     h_Energy_TrueVsReco[NSphereR];
TH1D*     h_RecoEnergy[NSphereR];
TH1D*     h_RecoMult[NSphereR];

void configure(){
  
  // open the file and set up the TTree
  TFile* file = new TFile(fFileName.c_str(),"READ");
  fTree = (TTree*)file->Get(fTreeName.c_str());
  setBranches(fTree);
  
  // make output file to store plots
  fOutFile = new TFile(fOutFileName.c_str(),"recreate");
  fOutFile->cd();

  // make histograms
  
  h_BlipE           = new TH1D("BlipE",     "Individual blip energies;Energy (MeV)", 200,0,10);
  h_BlipDist        = new TH1D("BlipDist",  "Individual blips;Distance (cm)",2000,0,50000);
  h_BlipTime        = new TH1D("BlipTime",  "Individual blips;Time (#mus)",1000,0,3000);
  h_MaxBlipE        = new TH1D("MaxBlipE",  "Maximum blip energy (``trunk'' candidate);Energy (MeV);Events",100,0.0,20.0);
  h_AveBlipE        = new TH1D("AveBlipE",  "Average blip energy;Energy (MeV);Events",100,0.0,5.);
  h_TotBlipE        = new TH1D("TotBlipE",  "Summed blip energy;Energy (MeV);Events", 300,0,30);
  h_Mult            = new TH1D("Mult",  "Blip multiplicity;Number of blips;Events", 100,0,100);
  h_NCaptPerEvt     = new TH1D("NCaptPerEvt","Neutron capture multiplicity;Number of captures;Events",10,0,10);
  h_NCaptTime       = new TH1D("NCaptTime",";Time of neutron capture (#mus);Events",1000,0,3000);
  h_NCaptDist       = new TH1D("NCaptDist",";Distance of neutron capture (cm);Events",100,0,50000);
  h_NCaptGammaE     = new TH1D("NCaptGammaE","Individual neutron capture #gamma's;Energy (MeV)",280,0,7);
  h_NCaptGammaESum  = new TH1D("NCaptGammaESum","Summed neutron capture #gamma's;Energy (MeV)",280,0,7);
  h_NScatGammaE     = new TH1D("NScatGammaE","Single gammas from neutron inelastic scatters;Energy (MeV)",280,0,7);

  for(int i=0; i<fSphereR.size(); i++){
    h_Energy_TrueVsReco[i]   = new TH2D(Form("Energy_TrueVsReco_%2.0fcm",fSphereR.at(i)),Form("R=%2.0f;True Neutron Energy (MeV);Reconstructed Energy (MeV)",fSphereR.at(i)),
                                  100,0,20,100,0,20);
    h_Energy_TrueVsReco[i]   ->SetOption("COLZ");
    
    h_RecoEnergy[i]   = new TH1D(Form("RecoEnergy_%2.0fcm",fSphereR.at(i)),Form("R=%2.0fcm;Reconstructed Energy (MeV)",fSphereR.at(i)),100,0,20);
    h_RecoMult[i]   = new TH1D(Form("RecoMult_%2.0fcm",fSphereR.at(i)),Form("R=%2.0fcm;Multiplicity",fSphereR.at(i)),50,0,50);
  }


}


// =============================================================
void NeutronReco_AnaMacro(){
  
  // configure histograms and TFile
  configure();
  reco();
  makePlots();

  // write histos
  h_BlipE         ->Write();
  h_BlipDist      ->Write();
  h_BlipTime      ->Write();
  h_MaxBlipE      ->Write();
  h_AveBlipE      ->Write();
  h_TotBlipE      ->Write();
  h_Mult          ->Write();
  h_NCaptPerEvt   ->Write();
  h_NCaptTime     ->Write();
  h_NCaptDist     ->Write();
  h_NCaptGammaE   ->Write();
  h_NCaptGammaESum->Write();
  h_NScatGammaE   ->Write();

  for(size_t i=0; i<fSphereR.size(); i++) h_Energy_TrueVsReco[i]->Write();
  for(size_t i=0; i<fSphereR.size(); i++) h_RecoEnergy[i]->Write();
  for(size_t i=0; i<fSphereR.size(); i++) h_RecoMult[i]->Write();

  fOutFile->Close();
  return;
}



// ============================================================
void reco(){
  
  // --------------------------------------------
  // loop over the events
  for(int iEvent=0; iEvent<fTree->GetEntries(); iEvent++){
    fTree->GetEntry(iEvent);
 
    // Reset class variables (event-wise)
    fNeutronEnergy  = -999.;

    // Counters
    int nCapturesPerEvent = 0;
    float nCaptGammaESum = 0.;

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
      

      // -------------------------------
      // Neutron checks
      if( pdg == 2112 ) {
        
        // Primary?
        if( proc == "primary" ) {
          fNeutronEnergy  = KE;
          fStartPoint     = loc;
        }
        
        // Does it capture?
        for(int ii=i; ii<nParticles; ii++){
          if( _Mother[ii] == trackId && _processname->at(ii) == "nCapture" ) {
            nCapturesPerEvent++;
            fNCaptures++;
            TVector3 capturePoint;
            capturePoint.SetXYZ(_StartPointx[ii],_StartPointy[ii],_StartPointz[ii]);
            h_NCaptDist->Fill( (capturePoint-fStartPoint).Mag() );
            h_NCaptTime->Fill( _StartT[ii]*1e-3 );
            break;
          }
        }

      }//<-- end neutron stuff


      // ------------------------------
      // Blip construction (electrons/positrons)
      if( fabs(pdg) == 11 ){
        EnergyDeposit b;
        b.Location = loc;
        b.Energy = edep;
        b.PathLength = dL;
        b.Time  = startT;
        v_blips.push_back(b);
        h_BlipTime->Fill(startT*1e-3);
      }

      // ----------------------------------------
      // Check for neutron capture gammas
      if( pdg == 22 && proc == "nCapture" ) {
        h_NCaptGammaE->Fill(KE);
        nCaptGammaESum += KE;
      }

      // ---------------------------------------
      // Check for neutron inelastic scatter gammas
      if( pdg == 22 && proc == "neutronInelastic" ) {
        h_NScatGammaE->Fill(KE);
      }
      
      // debugging output
      if(0){
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
      } 
    }//>> end particle loop

    // ---------------------------------
    // Per event neutron plots
    h_NCaptPerEvt->Fill(nCapturesPerEvent);
    if( nCapturesPerEvent == 1 ) h_NCaptGammaESum->Fill(nCaptGammaESum);

    // ---------------------------------
    // merging, thresholding, and smearing
    MergeBlips    (v_blips,fMinSep    );
    ThresholdBlips(v_blips,fThreshold );
    SmearBlips    (v_blips,fSmear     );
    
    printf("Event %5d -- smear: %6.2f keV, blips found: %d\n",iEvent,fSmear*1e3,(int)v_blips.size());


    // ----------------------------------
    // Sphere-drawing reconstruction 
    for(size_t i=0; i<fSphereR.size();i++) {
      fRecoEnergy[i] = 0.;
      fRecoMult[i] = 0;
    }
    // Find the largest blip
    float maxBlipE  = -999.;
    for(size_t i=0; i<v_blips.size();i++){
      if( v_blips.at(i).Energy > maxBlipE ) {
        maxBlipE = v_blips.at(i).Energy;
        fStartPoint = v_blips.at(i).Location; // !!!!!!! temp 2020-10-05
      }
    }
    float sumE      = 0.;
    for(size_t i=0; i<v_blips.size(); i++){
      EnergyDeposit thisBlip = v_blips.at(i);
      float E = thisBlip.Energy;
      sumE += E;
      float d = (thisBlip.Location - fStartPoint).Mag();
      
      for(size_t j=0; j<fSphereR.size(); j++ ){
        if( d < fSphereR.at(j) || fSphereR.at(j) < 0 ) {
          fRecoEnergy[j] += thisBlip.Energy;
          fRecoMult[j]++;
        }
      }
      h_BlipDist->Fill(d);
      h_BlipE->Fill(E);
    }

    // fill histograms
    for(size_t i=0; i<fSphereR.size(); i++){
      if( fNeutronEnergy > 0.2 ) h_Energy_TrueVsReco[i]->Fill(fNeutronEnergy,fRecoEnergy[i]);
      h_RecoEnergy[i]->Fill(fRecoEnergy[i]);
      h_RecoMult[i]->Fill(fRecoMult[i]);
    }
   
    h_TotBlipE->Fill(sumE);
    h_AveBlipE->Fill(sumE/v_blips.size()); 
    h_MaxBlipE->Fill(maxBlipE);
    h_Mult->Fill(v_blips.size());
    
    //-----------------------------------------------
    // Iterative blip 
    //BlipGrouping( v_blips );
  
  }

}



// =============================================================
void BlipGrouping(std::vector<EnergyDeposit> blips ){
  // Now that we've saved all the interesting blips, we will perform
  // the iterative sphere-drawing method where for gamma reconstruction
  //   (a) Identify highest-E blip
  //   (b) Draw sphere centered on that blip, and
  //       add up energies of all blips inside sphere
  //   (3) Locate highest-E blip of those that were not 
  //       grouped into the first sphere, and repeat,
  //       only grouping previously un-grouped blips.
  
  // Loop through the sphere sizes
  for(size_t iR = 0; iR < fSphereR.size(); iR++){
   
    // sphere radius
    float R = fSphereR.at(iR);

    // "ungroup" all the blips
    for(size_t i=0; i<blips.size(); i++) blips.at(i).isGrouped = false;
    size_t blips_grouped = 0;

    while( blips_grouped < blips.size() ) { 
      
      // sum E
      float sumE          = 0.;

      // find largest un-grouped blip
      float E_highest = -9.;
      TVector3 sphereCenter;
      for(size_t i=0; i<blips.size(); i++){
        EnergyDeposit thisBlip = blips.at(i);
        if( thisBlip.isGrouped ) continue;
        if( thisBlip.Energy > E_highest ) {
          E_highest = thisBlip.Energy;
          sphereCenter = thisBlip.Location;
        }
      }
    
      // keep count of the blips that go into this sphere
      int nblips_sphere = 0;

      // now group blips that fall in this sphere
      for(size_t i=0; i<blips.size(); i++){
        EnergyDeposit thisBlip = blips.at(i);
        if( thisBlip.isGrouped ) continue;
        float d = (thisBlip.Location - sphereCenter).Mag();
        if( d < R ) {
          blips.at(i).isGrouped = true;
          blips_grouped++;
          nblips_sphere++;
          sumE += thisBlip.Energy;
        }
      }
      
      // fill the histos
      // for the neutron cases, we throw out any sphere whose center
      // is > 2.5 m from the primary generation point (Whit's cut)
      //if( !(primPDG == 2112) || (sphereCenter-primLoc).Mag() < 100. ) {
        //h_Energy[iR]      ->Fill(sumE);
        //h_Energy_zoom[iR] ->Fill(sumE);
        //h_BlipMultSphere[iR]->Fill(nblips_sphere);
      //}

    }//<< end while(ungrouped blips)

  }//<< end loop over sphere radii

}


// ============================================================
void makePlots(){

/*
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
*/

}
