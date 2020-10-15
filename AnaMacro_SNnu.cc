//////////////////////////////////////////////////////////////////
// 
//  Macro for looking at neutrino CC and ES files from Ivan (their
//  variables differ from the ones from the SBND ana tree).
//
//  To run:
//  > root -l BlipReco_SNnu_vs_nCapture.cc
//
////////////////////////////////////////////////////////////////// 
#include "core/blip_nu.h"
//#include "core/tools.h"

// ===================   Parameters ===========================
std::string             fFileName     = "AnaTree_SNnue_CC.root";
//std::string             fFileName     = "AnaTree_SNnu_ES.root";
std::string             fTreeName     = "analysistree/anatree";
std::string             fOutFileName  = "plots.root";
std::string             fMCFile       = "../mcfiles/"+fFileName;         

float                   fThreshold    = 0.300; //0.075; //> Blip threshold (75 keV default)
float                   fSmear        = 0.00;   //> Smearing (50 keV usually)
float                   fMinSep       = 0.20;   //> Min blip separation used (0 = no merging)
float                   fSphereR      = 60.;    //> cm

//===============================================================
void AnaMacro_SNnu();
void configure();

int       n_sel_1MeV        = 0;
int       n_sel_4blips      = 0;
int       n_sel_4blips_1MeV = 0;
int       n_sel_2MeV        = 0;
int       n_sel_5blips      = 0;
int       n_sel_5blips_2MeV = 0;
int       n_sel_2blips      = 0;
int       n_sel_2blips_1MeV = 0;
int       n_sel_1blip       = 0;

TFile*    fOutFile;
TTree*    fTree;
TH1D*     h_NuEnergy;
TH1D*     h_BlipEnergy;
TH1D*     h_BlipL;
TH1D*     h_BlipMult;
TH1D*     h_MaxBlipE;
TH2D*     h_NuEnergy_vs_BlipMult;
TH2D*     h_NuEnergy_vs_MaxBlipE;
TH2D*     h_NuEnergy_vs_AveBlipE;

TH1D*     h_NeutronEnergy;
TH2D*     h_EnergyTrk_TrueVsReco;
TH2D*     h_EnergyTrk_TrueVsRes;
TH2D*     h_EnergyTotal_TrueVsReco;
TH2D*     h_EnergyTotal_TrueVsRes;
TH1D*     h_BlipMult_Sphere;
TH2D*     h_Mult_vs_Energy;
TH2D*     h_Mult_vs_MaxBlipE;
TH2D*     h_Mult_vs_AveBlipE;
TH2D*     h_Mult_vs_Ratio;
TH2D*     h_Mult_vs_ElEnergy;
TH1D*     h_Energy_vs_Res_Trk;
TH1D*     h_Energy_vs_Res_Total;
TH1D*     h_Energy_vs_RMS_Trk;
TH1D*     h_Energy_vs_RMS_Total;

//================================================================
void configure(){
  // make the histograms
  h_NuEnergy        = new TH1D("NuEnergy","Neutrino Energy;Energy (MeV);Entries",60,0.0,60.0);
  h_BlipEnergy      = new TH1D("BlipEnergy","Individual Blip Energies;Individual Blip Energy (MeV);Entries",3000,0.0,30.0);
  h_BlipL        = new TH1D("BlipL","Electron path length;Path Length (cm);Entries",500,0.0,5.0);
  h_MaxBlipE        = new TH1D("MaxBlipE","Maximum Blip Energy (``trunk'' candidate);Energy (MeV);Entries",3000,0.0,30.0);
  h_BlipMult        = new TH1D("BlipMult","Blip Multiplicity;Number of Blips Per Event",100,0,100);
  h_BlipMult_Sphere = new TH1D("BlipMult_Sphere","Blip Multiplicity (within 30cm sphere);Number of Blips Per Event",100,0,100);
  h_NuEnergy_vs_BlipMult = new TH2D("NuEnergy_vs_BlipMult",";Neutrino Energy (MeV);Blip Multiplicity per Event",60,0,60,50,0,50);
  h_NuEnergy_vs_AveBlipE = new TH2D("NuEnergy_vs_AveBlipE",";Neutrino Energy (MeV);Average Blip Energy (MeV)",120,0,60,100,0,5);
  h_NuEnergy_vs_MaxBlipE = new TH2D("NuEnergy_vs_MaxBlipE",";Neutrino Energy (MeV);Max Blip Energy (MeV)",120,0,60,120,0,60);

  h_Mult_vs_Energy    = new TH2D("Mult_vs_Energy",";Blip Multiplicity;Summed Blip Energy (MeV)",25,0,25,40,0,20.);
  h_Mult_vs_MaxBlipE  = new TH2D("Mult_vs_MaxBlipE",";Blip Multiplicity;Max Blip Energy (MeV)",25,0,25,25,0,5.);
  h_Mult_vs_AveBlipE  = new TH2D("Mult_vs_AveBlipE",";Blip Multiplicity;Average Blip Energy (MeV)",25,0,25,25,0,5.);
  h_Mult_vs_Ratio = new TH2D("Mult_vs_Ratio",";Blip Multiplicity;Ratio of E_{blips} to E_e",25,0,25,25,0,1.);
  h_Mult_vs_ElEnergy = new TH2D("Mult_vs_ElEnergy",";Blip Multiplicity;Primary Electron Energy Deposited (MeV)",25,0,25,100,0,50.);
  h_NuEnergy_vs_BlipMult->SetOption("colz");
  h_NuEnergy_vs_AveBlipE->SetOption("colz");
  h_NuEnergy_vs_MaxBlipE->SetOption("colz");
  h_Mult_vs_Energy ->SetOption("colz");
  h_Mult_vs_MaxBlipE->SetOption("colz");
  h_Mult_vs_AveBlipE->SetOption("colz");
  h_Mult_vs_Ratio->SetOption("colz");
  h_Mult_vs_ElEnergy->SetOption("colz");
  
  float resbins_energy    = 24;
  float resbins           = 400;
  h_EnergyTrk_TrueVsReco  = new TH2D("EnergyTrk_TrueVsReco","Electron Track;True Energy (MeV);Reconstructed Track Energy (MeV)",
                                      //120,0,60,120,0,60);
                                      resbins_energy,0,60,120,0,60);
  h_EnergyTrk_TrueVsRes   = new TH2D("EnergyTrk_TrueVsRes","Electron Track;True Energy (MeV);Resolution",
                                      resbins_energy,0,60,resbins,-1.,1.);
  h_EnergyTotal_TrueVsReco  = new TH2D("EnergyTotal_TrueVsReco","Electron Track + Blips (< 30cm);True Energy (MeV);Reconstructed Energy (MeV)",
                                      //120,0,60,120,0,60);
                                      resbins_energy,0,60,120,0,60);
  h_EnergyTotal_TrueVsRes   = new TH2D("EnergyTotal_TrueVsRes","Electron Track + Blips (< 30cm);True Energy (MeV);Resolution",
                                      resbins_energy,0,60,resbins,-1.,1.);
  h_Energy_vs_Res_Trk     = new TH1D("Energy_vs_Res_Trk","Electron Track;True Neutrino Energy (MeV);RMS Resolution",
                                      resbins_energy,0,60);
  h_Energy_vs_Res_Total     = new TH1D("Energy_vs_Res_Total","Electron Track;True Neutrino Energy (MeV);RMS Resolution",
                                      resbins_energy,0,60);
  h_Energy_vs_RMS_Trk     = new TH1D("Energy_vs_RMS_Trk","Electron Track;True Neutrino Energy (MeV);Reconstructed Energy RMS",
                                      resbins_energy,0,60);
  h_Energy_vs_RMS_Total     = new TH1D("Energy_vs_RMS_Total","Electron Track;True Neutrino Energy (MeV);Reconstructed Energy RMS",
                                      resbins_energy,0,60);
  h_EnergyTrk_TrueVsReco  ->SetOption("colz");
  h_EnergyTrk_TrueVsRes   ->SetOption("colz");
  h_EnergyTotal_TrueVsReco  ->SetOption("colz");
  h_EnergyTotal_TrueVsRes   ->SetOption("colz");

    
  // open the file and set up the TTree
  TFile* file = new TFile(fMCFile.c_str(),"READ");
  fTree = (TTree*)file->Get(fTreeName.c_str());
  setBranches(fTree);
}



// =============================================================
void AnaMacro_SNnu(){

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
      float E       = 1e3*_Eng[i];
      float mass    = 1e3*Mass(i);
      float dL      = (loc-locEnd).Mag();
      
      // calculate the energy deposited by this particle
      // (includes depositions from contiguous daughters
      // like "eIoni" electrons)   
      float edep = CalcEnergyDep(i);

      // record primary electron info
      if( fabs(pdg) == 11 && proc == "primary" ) {
          elEnergy    = E;
          elEnergyDep = edep;
          nuVert      = loc;
      }

      // make blip if electron or proton with edep < 3 MeV
      if(     proc != "primary" && edep > 0 
          &&  ( fabs(pdg) == 11 || (fabs(pdg) == 2212 && edep < 3. )) ){
          EnergyDeposit b;
          b.Energy    = edep;
          b.Location  = loc;
          b.PathLength= dL;
          b.isGrouped = false;
          v_blips.push_back(b);
      }
      
      // (for debugging output -- keep commented out during normal running)
      if(1){
        printf("  %3i PDG: %10i, dL=%5.2f, E=%8.3f,Edep=%8.3f, moth=%3i, %12s, Ndaught=%3i, blips=%3lu\n",
          trackId,
          pdg,
          dL,
          E-mass,
          edep,
          mother,
          proc.c_str(),
          nD,
          v_blips.size());
      }
      

    }//>> end particle loop
    
    // .................................
    // merge close-together blips
    MergeBlips(v_blips,fMinSep);
   
    // .................................
    // do thresholding
    ThresholdBlips(v_blips,fThreshold);
    
    h_NuEnergy->Fill(nuEnergy);    
    h_BlipMult->Fill(v_blips.size());


    // ===========================================
    // Sphere drawing stuff
    //
    float totalBlipE        = 0.;
    float totalBlipE_sphere = 0.; //elEnergyDep;
    float maxBlipE          = 0.;
    int   nBlips            = v_blips.size();
    int   nBlips_sphere     = 0;
    for(int i=0; i<nBlips; i++){
      float E = v_blips.at(i).Energy;
      h_BlipEnergy  ->Fill(E);
      h_BlipL       ->Fill(v_blips.at(i).PathLength);
      totalBlipE    += E;
      TVector3 loc  = v_blips.at(i).Location;
      TVector3 d    = (loc-nuVert);
      if( d.Mag() < fSphereR && d.Mag() > 0.5 ){
        if( E > maxBlipE ) maxBlipE = E;
        nBlips_sphere++;
        totalBlipE_sphere += E;
      }
    }
    float aveBlipE = 0.;
    if( nBlips_sphere ) aveBlipE = totalBlipE_sphere / nBlips_sphere;

    /*
    // -------------------------
    // Ar39 smearing
    int dN = 0;
    float dE = 0.;
    if( fSphereR == 30.  ) {
      dN = std::max(0,(int)(fRand->Gaus(0.3,0.6)+0.5));
      dE = std::max(0.,fRand->Gaus(0.1,0.15));
    } else
    if( fSphereR == 60. ) {
      dN = std::max(0,(int)(fRand->Gaus(2.5,1.6)+0.5));
      dE = std::max(0.,fRand->Gaus(0.65,0.45));
    }
      nBlips_sphere += dN;
      totalBlipE_sphere += dE;
    */

    // channel ID cuts
    float nn = nBlips_sphere;
    float ee = totalBlipE_sphere;
    /*
    if( nn < 4 )            n_sel_4blips++;
    if( nn < 4 && ee < 1. ) n_sel_4blips_1MeV++;
    if( ee < 1. )           n_sel_1MeV++;
    if( ee < 2. )           n_sel_2MeV++;
    if( nn < 5 )            n_sel_5blips++;
    if( nn < 5 && ee < 2.)  n_sel_5blips_2MeV++;
    if( nn < 2 )            n_sel_2blips++;
    if( nn < 2 && ee < 1. ) n_sel_2blips_1MeV++;
    if( nn < 1 )            n_sel_1blip++;
    */
    if( nn < 3 )            n_sel_4blips++;
    if( nn < 3 && ee < 1. ) n_sel_4blips_1MeV++;
    if( ee < 1. )           n_sel_1MeV++;
    if( ee < 1.5 )           n_sel_2MeV++;
    if( nn < 4 )            n_sel_5blips++;
    if( nn < 4 && ee < 1.5)  n_sel_5blips_2MeV++;
    if( nn < 2 )            n_sel_2blips++;
    if( nn < 2 && ee < 1. ) n_sel_2blips_1MeV++;
    if( nn < 1 )            n_sel_1blip++;
    
    
    h_BlipMult_Sphere->Fill(nBlips_sphere);
    h_MaxBlipE->Fill(maxBlipE); 
    h_NuEnergy_vs_AveBlipE->Fill(nuEnergy,totalBlipE/v_blips.size());
    h_NuEnergy_vs_BlipMult->Fill(nuEnergy,v_blips.size());
    h_NuEnergy_vs_MaxBlipE->Fill(nuEnergy,maxBlipE);
    h_EnergyTrk_TrueVsReco->Fill( nuEnergy,elEnergyDep);
    h_EnergyTrk_TrueVsRes ->Fill( nuEnergy, (elEnergyDep-nuEnergy)/nuEnergy);
    h_EnergyTotal_TrueVsReco->Fill( nuEnergy,totalBlipE_sphere+elEnergyDep);
    h_EnergyTotal_TrueVsRes ->Fill( nuEnergy, (totalBlipE_sphere+elEnergyDep-nuEnergy)/nuEnergy);
    h_Mult_vs_Energy->Fill(nBlips_sphere,totalBlipE_sphere);
    h_Mult_vs_MaxBlipE->Fill(nBlips_sphere,maxBlipE); 
    h_Mult_vs_AveBlipE->Fill(nBlips_sphere,aveBlipE); 
    h_Mult_vs_Ratio->Fill(nBlips_sphere,totalBlipE_sphere/elEnergyDep);
    h_Mult_vs_ElEnergy->Fill(nBlips_sphere,elEnergyDep);

    printf("MC file: %s, event %5d, blips found: %d\n",fFileName.c_str(), iEvent,(int)v_blips.size());

  }//>> end loop over events


  // ------------------------------------------------
  // Channel ID counts
  /*
  std::cout
  <<"=====================================\n"
  <<" N< --     E< 1 MeV    Nsel = "<<n_sel_1MeV<<"\n"
  <<" N< 4      E< -----    Nsel = "<<n_sel_4blips<<"\n"
  <<" N< 4      E< 1 MeV    Nsel = "<<n_sel_4blips_1MeV<<"\n"
  <<" ------------------------------------\n"
  <<" N< --     E< 2 MeV    Nsel = "<<n_sel_2MeV<<"\n"
  <<" N< 5      E< -----    Nsel = "<<n_sel_5blips<<"\n"
  <<" N< 5      E< 2 MeV    Nsel = "<<n_sel_5blips_2MeV<<"\n"
  <<" ------------------------------------\n"
  <<" N< --     E< 1 MeV    Nsel = "<<n_sel_1MeV<<"\n"
  <<" N< 2      E< -----    Nsel = "<<n_sel_2blips<<"\n"
  <<" N< 2      E< 1 MeV    Nsel = "<<n_sel_2blips_1MeV<<"\n"
  <<" N< 1      E< -----    Nsel = "<<n_sel_1blip<<"\n"
  <<"=====================================\n";
  */
  std::cout
  <<"=====================================\n"
  <<" N< --     E< 1 MeV    Nsel = "<<n_sel_1MeV<<"\n"
  <<" N< 3      E< -----    Nsel = "<<n_sel_4blips<<"\n"
  <<" N< 3      E< 1 MeV    Nsel = "<<n_sel_4blips_1MeV<<"\n"
  <<" ------------------------------------\n"
  <<" N< --     E< 1.5 MeV  Nsel = "<<n_sel_2MeV<<"\n"
  <<" N< 4      E< -----    Nsel = "<<n_sel_5blips<<"\n"
  <<" N< 4      E< 1.5 MeV  Nsel = "<<n_sel_5blips_2MeV<<"\n"
  <<" ------------------------------------\n"
  <<" N< --     E< 1 MeV    Nsel = "<<n_sel_1MeV<<"\n"
  <<" N< 2      E< -----    Nsel = "<<n_sel_2blips<<"\n"
  <<" N< 2      E< 1 MeV    Nsel = "<<n_sel_2blips_1MeV<<"\n"
  <<" N< 1      E< -----    Nsel = "<<n_sel_1blip<<"\n"
  <<"=====================================\n";


  TCanvas* c1 = new TCanvas("c1","c1",600,450);
  gStyle->SetOptStat(0);
  gPad->SetMargin(0.10, 0.12, 0.12, 0.08);
  h_Mult_vs_Energy->GetXaxis()->CenterTitle();
  h_Mult_vs_Energy->GetYaxis()->CenterTitle();
  h_Mult_vs_Energy->GetXaxis()->SetLabelSize(0.045);
  h_Mult_vs_Energy->GetXaxis()->SetTitleSize(0.045);
  h_Mult_vs_Energy->GetYaxis()->SetLabelSize(0.045);
  h_Mult_vs_Energy->GetYaxis()->SetTitleSize(0.045);
  h_Mult_vs_Energy->GetZaxis()->SetLabelSize(0.045);
  h_Mult_vs_Energy->GetZaxis()->SetTitleSize(0.045);
  h_Mult_vs_Energy->Draw();
  

  /*
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

  */



  
  // ------------------------------------------------- 
  // make output file to store plots
  fOutFile = new TFile(fOutFileName.c_str(),"recreate");
  fOutFile->cd();
  
  // save all the lower-level histograms
  h_NuEnergy->Write();
  h_BlipEnergy->Write();
  h_BlipL->Write();
  h_BlipMult->Write();
  h_BlipMult_Sphere->Write();
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
  h_Mult_vs_Energy->Write();
  h_Mult_vs_MaxBlipE->Write();
  h_Mult_vs_AveBlipE->Write();
  h_Mult_vs_Ratio->Write();
  h_Mult_vs_ElEnergy->Write();
  fOutFile->Close();
  return;
}

