//////////////////////////////////////////////////////////////////
// 
//  BlipReco_SNnu_vs_nCapture
//
//  Macro for looking at neutrino CC and ES files from Ivan (their
//  variables differ from the ones from the SBND ana tree).
//
//  To run:
//  > root -l BlipReco_SNnu_vs_nCapture.cc
//
////////////////////////////////////////////////////////////////// 
#include "blip_nu.h"
#include "tools.h"

// ===================   Parameters ===========================
//std::string             fFileName     = "../mcfiles/anatree_neutrons_1eV.root";
std::string             fFileName     = "../mcfiles/AnaTree_SNnueCC.root";
std::string             fTreeName     = "analysistree/anatree";
std::string             fOutFileName  = "plots.root";

float                   fThreshold    = 0.0750; //> Blip threshold (75 keV default)
float                   fSmear        = 0.00;   //> Smearing (50 keV usually)
float                   fMinSep       = 0.20;   //> Min blip separation used 
                                                //  during blip merging stage

//===============================================================
void SNnu();
void configure();

TFile*    fOutFile;
TTree*    fTree;
TH1D*     h_NuEnergy;
TH1D*     h_BlipEnergy;
TH1D*     h_BlipMult;
TH1D*     h_MaxBlipE;
TH2D*     h_NuEnergy_vs_BlipMult;
TH2D*     h_NuEnergy_vs_AveBlipE;
TH2D*     h_NuEnergy_vs_MaxBlipE;

TH1D*     h_NeutronEnergy;
  
TH2D*     h_EnergyTrk_TrueVsReco;
TH2D*     h_EnergyTrk_TrueVsRes;
TH2D*     h_EnergyTotal_TrueVsReco;
TH2D*     h_EnergyTotal_TrueVsRes;

TH1D*     h_Energy_vs_Res_Trk;
TH1D*     h_Energy_vs_Res_Total;

TH1D*     h_Energy_vs_RMS_Trk;
TH1D*     h_Energy_vs_RMS_Total;

//================================================================
void configure(){
  // make the histograms
  h_NuEnergy        = new TH1D("NuEnergy","Neutrino Energy;Energy (MeV);Entries",60,0.0,60.0);
  h_BlipEnergy      = new TH1D("BlipEnergy","Individual Blip Energies;Individual Blip Energy (MeV);Entries",3000,0.0,30.0);
  h_MaxBlipE        = new TH1D("MaxBlipE","Maximum Blip Energy (``trunk'' candidate);Energy (MeV);Entries",3000,0.0,30.0);
  h_BlipMult        = new TH1D("BlipMult","Blip Multiplicity;Number of Blips Per Event",100,0,100);
  h_NuEnergy_vs_BlipMult = new TH2D("NuEnergy_vs_BlipMult",";Neutrino Energy (MeV);Blip Multiplicity per Event",60,0,60,50,0,50);
  h_NuEnergy_vs_AveBlipE = new TH2D("NuEnergy_vs_AveBlipE",";Neutrino Energy (MeV);Average Blip Energy (MeV)",120,0,60,100,0,10);
  h_NuEnergy_vs_MaxBlipE = new TH2D("NuEnergy_vs_MaxBlipE",";Neutrino Energy (MeV);Max Blip Energy (MeV)",120,0,60,120,0,60);
  h_NuEnergy_vs_BlipMult->SetOption("colz");
  h_NuEnergy_vs_AveBlipE->SetOption("colz");
  h_NuEnergy_vs_MaxBlipE->SetOption("colz");
  
  float resbins_energy    = 24;
  float resbins           = 400;
  h_EnergyTrk_TrueVsReco  = new TH2D("EnergyTrk_TrueVsReco","Electron Track;True Energy (MeV);Reconstructed Track Energy (MeV)",
                                      //120,0,60,120,0,60);
                                      resbins_energy,0,60,120,0,60);
  h_EnergyTrk_TrueVsRes   = new TH2D("EnergyTrk_TrueVsRes","Electron Track;True Energy (MeV);Resolution",
                                      resbins_energy,0,60,resbins,-1.,1.);
  h_EnergyTrk_TrueVsReco  ->SetOption("colz");
  h_EnergyTrk_TrueVsRes   ->SetOption("colz");
  h_EnergyTotal_TrueVsReco  = new TH2D("EnergyTotal_TrueVsReco","Electron Track + Blips (< 30cm);True Energy (MeV);Reconstructed Energy (MeV)",
                                      //120,0,60,120,0,60);
                                      resbins_energy,0,60,120,0,60);
  h_EnergyTotal_TrueVsRes   = new TH2D("EnergyTotal_TrueVsRes","Electron Track + Blips (< 30cm);True Energy (MeV);Resolution",
                                      resbins_energy,0,60,resbins,-1.,1.);
  h_EnergyTotal_TrueVsReco  ->SetOption("colz");
  h_EnergyTotal_TrueVsRes   ->SetOption("colz");

  h_Energy_vs_Res_Trk     = new TH1D("Energy_vs_Res_Trk","Electron Track;True Neutrino Energy (MeV);RMS Resolution",
                                      resbins_energy,0,60);
  h_Energy_vs_Res_Total     = new TH1D("Energy_vs_Res_Total","Electron Track;True Neutrino Energy (MeV);RMS Resolution",
                                      resbins_energy,0,60);
  h_Energy_vs_RMS_Trk     = new TH1D("Energy_vs_RMS_Trk","Electron Track;True Neutrino Energy (MeV);Reconstructed Energy RMS",
                                      resbins_energy,0,60);
  h_Energy_vs_RMS_Total     = new TH1D("Energy_vs_RMS_Total","Electron Track;True Neutrino Energy (MeV);Reconstructed Energy RMS",
                                      resbins_energy,0,60);

    
  // open the file and set up the TTree
  TFile* file = new TFile(fFileName.c_str(),"READ");
  fTree = (TTree*)file->Get(fTreeName.c_str());
  setBranches(fTree);
}


// =============================================================
void SNnu(){

  // configure histograms and TFile
  configure();
  
  // --------------------------------------------
  // loop over the events
  for(int iEvent=0; iEvent<fTree->GetEntries(); iEvent++){
    fTree->GetEntry(iEvent);
    
    float nuEnergy = _NuEnergy*1e3;
    TVector3 nuVert;

    // Make a vector of Blip objects to fill
    std::vector<EnergyDeposit> v_blips;
      
    float elEnergy    = -999.;
    float elEnergyDep = -999.;
    float totalBlipE  = 0.;
    float maxBlipE    = -999.;

    // Loop over particles and fill in the blips.
    int nParticles = _geant_list_size;
    for(int i=0; i<nParticles; i++){

      TVector3 loc(_StartPointx[i],_StartPointy[i],_StartPointz[i]);
      TVector3 locEnd(_EndPointx[i],_EndPointy[i],_EndPointz[i]);
      int pdg       = _pdg[i];
      int trackId   = _TrackId[i]; 
      std::string proc = _processname->at(i);
      int mother    = _Mother[i];
      int nD        = _NumberDaughters[i];
      //float dL      = _pathlen[i];
     
      // if electron, make blip
      if( fabs(pdg) == 11 ) {
        float edep = CalcEnergyDep(i);
      
        if( proc == "primary" ) {
          elEnergy = _Eng[i]*1e3;
          elEnergyDep = edep;
          nuVert = loc;
        } 
        else { 
        EnergyDeposit b;
          b.Energy    = edep;
          b.Location  = loc;
          b.isGrouped = false;
          v_blips.push_back(b);
        }
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
   
    // .................................
    // do thresholding
    ThresholdBlips(v_blips,fThreshold);

    float totalE_30cm = elEnergyDep;
    
    h_NuEnergy->Fill(nuEnergy);    
    h_BlipMult->Fill(v_blips.size());
    int N = v_blips.size();
    for(int i=0; i<N; i++){
      float E = v_blips.at(i).Energy;
      h_BlipEnergy->Fill(E);
      totalBlipE += E;
      if( E > maxBlipE ) maxBlipE = E;
      if( (nuVert-v_blips.at(i).Location).Mag() < 30 ) totalE_30cm += E;
    }
   
    h_MaxBlipE->Fill(maxBlipE); 
    h_NuEnergy_vs_AveBlipE->Fill(nuEnergy,totalBlipE/v_blips.size());
    h_NuEnergy_vs_BlipMult->Fill(nuEnergy,v_blips.size());
    h_NuEnergy_vs_MaxBlipE->Fill(nuEnergy,maxBlipE);

    h_EnergyTrk_TrueVsReco->Fill( nuEnergy,elEnergyDep);
    h_EnergyTrk_TrueVsRes ->Fill( nuEnergy, (elEnergyDep-nuEnergy)/nuEnergy);
    h_EnergyTotal_TrueVsReco->Fill( nuEnergy,totalE_30cm);
    h_EnergyTotal_TrueVsRes ->Fill( nuEnergy, (totalE_30cm-nuEnergy)/nuEnergy);

    // ................................
    // do smearing
    SmearBlips(v_blips,fSmear);

    printf("Event %5d -- smear: %6.2f keV, blips found: %d\n",iEvent,fSmear*1e3,(int)v_blips.size());

    // ===================================
    // Make histograms



  }//>> end loop over events

  

  // ------------------------------------------------
  // RMS resolution histograms
  TH1D* h_pfx_trk = h_EnergyTrk_TrueVsRes->ProfileX("pfx_trk",1,-1,"s");
  TH1D* h_pfx_tot = h_EnergyTotal_TrueVsRes->ProfileX("pfs_tot",1,-1,"s");

  for(size_t i=0; i<h_EnergyTrk_TrueVsRes->GetXaxis()->GetNbins();i++){
    // only plot 5 to 50 MeV
    float En = h_EnergyTrk_TrueVsRes->GetXaxis()->GetBinCenter(i);
    if( En < 5 || En > 50 ) continue;

    std::cout<<"Looking at bin "<<i<<" corresponding to energy "<<En<<" MeV\n";
    
    TH1D* h_pjy_trk = h_EnergyTrk_TrueVsRes->ProjectionY("blah1",i,i);
    TH1D* h_pjy_tot = h_EnergyTotal_TrueVsRes->ProjectionY("blah2",i,i);
    

    float bias_trk = h_pfx_trk->GetBinContent(i);
    float bias_tot = h_pfx_tot->GetBinContent(i);
    float res_trk = h_pjy_trk->GetRMS(); //h_pfx_trk->GetBinError(i);
    float res_tot = h_pjy_tot->GetRMS(); //h_pfx_tot->GetBinError(i);
    h_Energy_vs_Res_Trk->SetBinContent(i,res_trk);
    h_Energy_vs_Res_Total->SetBinContent(i,res_tot);
  
    TH1D* h_p_trk = h_EnergyTrk_TrueVsReco->ProjectionY("blah3",i,i);
    TH1D* h_p_tot = h_EnergyTotal_TrueVsReco->ProjectionY("blah4",i,i);
    TCanvas* c1 = new TCanvas(Form("bin %lu",i),Form("bin %lu",i));
    h_p_tot->DrawCopy();
    h_Energy_vs_RMS_Trk   ->SetBinContent(i,h_p_trk->GetRMS()/h_p_trk->GetMean());
    h_Energy_vs_RMS_Total ->SetBinContent(i,h_p_tot->GetRMS()/h_p_tot->GetMean());

  }




  
  // ------------------------------------------------- 
  // make output file to store plots
  fOutFile = new TFile(fOutFileName.c_str(),"recreate");
  fOutFile->cd();
  
  // save all the lower-level histograms
  h_NuEnergy->Write();
  h_BlipEnergy->Write();
  h_BlipMult->Write();
  h_MaxBlipE->Write();
  h_NuEnergy_vs_AveBlipE->Write(); 
  h_NuEnergy_vs_BlipMult->Write();
  h_NuEnergy_vs_MaxBlipE->Write();
  h_EnergyTrk_TrueVsReco->Write();
  h_EnergyTrk_TrueVsRes->Write();
  h_EnergyTotal_TrueVsReco->Write();
  h_EnergyTotal_TrueVsRes->Write();
  h_Energy_vs_Res_Trk->Write();
  h_Energy_vs_Res_Total->Write();
  h_Energy_vs_RMS_Trk->Write();
  h_Energy_vs_RMS_Total->Write();
  fOutFile->Close();
  return;
}

