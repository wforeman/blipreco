//////////////////////////////////////////////////////////////////
//  
//  NEUTRON ANALYSIS
//
//  This macro is specialized to look at the MC neutron events and
//  look at various diagnostics as well as resolution slices.
//
//  To run:
//  > root -l macro.cc
//
////////////////////////////////////////////////////////////////// 
#include "core/blip.h"
#include "core/tools.h"

// ============================================================
// ===================   Parameters ===========================
//std::string             fFileName    = "../mcfiles/anatree_neutrons_1eV.root";
//std::string             fFileName    = "../mcfiles/anatree_neutron_0to20MeV_larworld.root";
std::string             fFileName    = "../mcfiles/AnaTree_neutron_10MeV.root";

std::string             fTreeName     = "analysistree/anatree";
std::string             fOutFileName  = "plots.root";
float                   fThreshold    = 0.100; // MeV (default 0.075)
float                   fSmear        = 0.00; // 50 keV
float                   fMinSep       = 0.2; // cm
std::vector<float>      fSphereR      = {30., 60.};
#define NSphereR 2

//===============================================================

void AnaMacro_Neutron();
void configure();
void reco();
void makePlots();

TVector3  fPrimaryStartPoint;
float     fNeutronEnergy;
float     fRecoEnergy[NSphereR];
int     fRecoMult[NSphereR];

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
TH1D*     h_NScatPerEvt;
TH2D*     h_Energy_vs_NScatPerEvt;
TH2D*     h_IncEnergy_vs_NScatGammaESum;
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
  h_NScatPerEvt     = new TH1D("NScatPerEvt","Number of inelastic neutron scatters per event;Number of scatters",10,0,10);
  h_NCaptTime       = new TH1D("NCaptTime",";Time of neutron capture (#mus);Events",1000,0,3000);
  h_NCaptDist       = new TH1D("NCaptDist",";Distance of neutron capture (cm);Events",1000,0,50000);
  h_NCaptGammaE     = new TH1D("NCaptGammaE","Individual neutron capture #gamma's;Energy (MeV)",280,0,7);
  h_NCaptGammaESum  = new TH1D("NCaptGammaESum","Summed neutron capture #gamma's;Energy (MeV)",280,0,7);
  h_NScatGammaE     = new TH1D("NScatGammaE","Single gammas from neutron inelastic scatters;Energy (MeV)",280,0,7);
  h_Energy_vs_NScatPerEvt = new TH2D("Energy_vs_NScatPerEvt",";Primary neutron energy (MeV);Number of inelastic scatters",100,0,20,10,0,10);
  h_IncEnergy_vs_NScatGammaESum  = new TH2D("IncEnergy_vs_GammaESum","Neutron inelastic scatters;Incident neutron energy (MeV);Summed E_{#gamma} from interaction (MeV)",100,0,20,200,0,20);
  
  h_Energy_vs_NScatPerEvt->SetOption("colz");
  h_IncEnergy_vs_NScatGammaESum->SetOption("colz");

  for(int i=0; i<fSphereR.size(); i++){
    h_Energy_TrueVsReco[i]   = new TH2D(Form("Energy_TrueVsReco_%2.0fcm",fSphereR.at(i)),Form("R=%2.0f;True Neutron Energy (MeV);Reconstructed Energy (MeV)",fSphereR.at(i)),
                                  100,0,20,100,0,20);
    h_Energy_TrueVsReco[i]   ->SetOption("COLZ");
    
    h_RecoEnergy[i]   = new TH1D(Form("RecoEnergy_%2.0fcm",fSphereR.at(i)),Form("R=%2.0fcm;Reconstructed Energy (MeV)",fSphereR.at(i)),100,0,20);
    h_RecoMult[i]   = new TH1D(Form("RecoMult_%2.0fcm",fSphereR.at(i)),Form("R=%2.0fcm;Multiplicity",fSphereR.at(i)),50,0,50);
  }


}


// =============================================================
void AnaMacro_Neutron(){
  
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
  h_NScatPerEvt   ->Write();
  h_NCaptTime     ->Write();
  h_NCaptDist     ->Write();
  h_NCaptGammaE   ->Write();
  h_NCaptGammaESum->Write();
  h_NScatGammaE   ->Write();
  h_Energy_vs_NScatPerEvt->Write();
  h_IncEnergy_vs_NScatGammaESum->Write();

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
 
    // Event-wise values to keep track of
    fNeutronEnergy        = -999.;
    int nCapturesPerEvent = 0;
    int nScattersPerEvent = 0;

    // Make a vector of blip objects to fill,
    std::vector<EnergyDeposit> v_blips;
    float total_edep = 0.;
    
    // Loop over particles and fill in the blips.
    int nParticles = _geant_list_size;
    for(int i=0; i<nParticles; i++){

      TVector3 loc(_StartPointx[i],_StartPointy[i],_StartPointz[i]);
      //TVector3 locEnd(_EndPointx[i],_EndPointy[i],_EndPointz[i]);
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
      


      // =====================================================
      // Neutron checks
      
      if( pdg == 2112 ) {

        // Primary?
        if( proc == "primary" ) {
          fNeutronEnergy      = KE;
          fPrimaryStartPoint  = loc;
        }

        bool  isScatter = false;
        bool  isCapture = false;
        float gammaEsum_capt = 0;
        float gammaEsum_scat = 0;

        // Note that at low enough energy, < ~2.9 MeV, the SAME neutron can
        // both scatter and capture, without a new particle being created in
        // the list.

        // Figure out what the neutron eventually does
        // by looking for daughters
        for(int ii=i; ii<nParticles; ii++){
          if( _pdg[ii] == 22 && _Mother[ii] == trackId ) {

            float Eg = 1e3*_Eng[ii];

            // CAPTURE?
            if( _processname->at(ii) == "nCapture" ) {
              gammaEsum_capt += Eg;
              h_NCaptGammaE->Fill(Eg);
              if( !isCapture ) {
                isCapture = true;
                nCapturesPerEvent++;
                TVector3 capturePoint;
                capturePoint.SetXYZ(_StartPointx[ii],_StartPointy[ii],_StartPointz[ii]);
                h_NCaptDist->Fill( (capturePoint-fPrimaryStartPoint).Mag() );
                h_NCaptTime->Fill( _StartT[ii]*1e-3 );
              }
            }

            //INELASTIC SCATTER?
            if( _processname->at(ii) == "neutronInelastic" ) {
              gammaEsum_scat += Eg;
              h_NScatGammaE->Fill(Eg);
              if( !isScatter ) {
                isScatter = true;
                nScattersPerEvent++;
              }
            }
            
          }
        }//< end determination of capture/scatter

        if( isCapture ) h_NCaptGammaESum->Fill(gammaEsum_capt);
        if( isScatter ) h_IncEnergy_vs_NScatGammaESum->Fill(KE,gammaEsum_scat);

      }//<-- end neutron stuff
      // =====================================================




      // =====================================================
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
      // =====================================================
      

      // debugging output
      if(0){
        printf("  %3i PDG: %10i,  dL=%8.3f, d0=%8.3f, KE0=%8.3f, Edep=%8.3f, T0=%8.3f, Tf=%8.3f, moth=%3i, %12s, Ndaught=%i\n",
        trackId,
        pdg,
        dL,
        (loc-fPrimaryStartPoint).Mag(),
        KE,
        //KEf,
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
    h_NScatPerEvt->Fill(nScattersPerEvent);
    h_Energy_vs_NScatPerEvt->Fill(fNeutronEnergy,nScattersPerEvent);

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

    TVector3 sphereCenter = fPrimaryStartPoint;

    // Find the largest blip
    float maxBlipE  = -999.;
    for(size_t i=0; i<v_blips.size();i++){
      if( v_blips.at(i).Energy > maxBlipE ) {
        maxBlipE = v_blips.at(i).Energy;
        sphereCenter = v_blips.at(i).Location;
      }
    }
    float sumE      = 0.;
    for(size_t i=0; i<v_blips.size(); i++){
      EnergyDeposit thisBlip = v_blips.at(i);
      float E = thisBlip.Energy;
      sumE += E;
      float d = (thisBlip.Location - sphereCenter).Mag();
      
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
    
  }

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
