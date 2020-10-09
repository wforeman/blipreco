//////////////////////////////////////////////////////////////////
// 
//  ElectronReco
//
//  Must have ROOT data file containing analysistree/anatree in 
//  same directory.
//
//  To run:
//  > root -l ElectronReco_AnaMacro.cc
//
////////////////////////////////////////////////////////////////// 
#include "blip.h"
#include "tools.h"

// ============================================================
// ===================   Parameters ===========================
#define nconfigs_SphereR 2
std::string             fFileName_0to60 = "../mcfiles/AnaTree_electron_0to60MeV_50kEvents.root";
std::string             fFileName_100   = "../mcfiles/AnaTree_electron_100MeV_10kEvents.root";
std::string             fFileName_200   = "../mcfiles/AnaTree_electron_200MeV_10kEvents.root";
std::string             fFileName_500   = "../mcfiles/AnaTree_electron_500MeV_10kEvents.root";
//std::string             fFileName_snnue = "../mcfiles/AnaTree_SNnueCC.root";
std::string             fFileName       = fFileName_0to60;

std::string             fTreeName       = "analysistree/anatree";
std::string             fOutFileName    = "plots.root";
float                   fThreshold      = 0.075; // 75 keV
float                   fSmear          = 0.0; // 50 keV
float                   fMinSep         = 0.2;   // 0.2 cm
std::vector<float>      v_SphereR       = { 60., 100. }; // cm

//===============================================================

void ElectronReco_AnaMacro();
void configure();
void reco();
void showerReco(std::vector<EnergyDeposit>& blips);
void makePlots();
void MakeSlices( const TH2D*,std::vector<TH1D*>&,std::vector<float>&,float,float);
void FillResolutionGraphs( TGraph*, TGraph*, std::vector<TH1D*>&,std::vector<float>&,std::string,bool,bool);

TFile*    fOutFile;
TTree*    fTree;

TVector3  fSphereCen;
float     fTrueEnergy;
float     fTrkEnergyDep;
//float     fShwEnergyDep[nconfigs_SphereR];
float     fTotalEnergyDep[nconfigs_SphereR];

TH1D*     h_EnergyBudget;
TH1D*     h_BlipEnergy;
TH1D*     h_MaxBlipEnergy;
TH1D*     h_BlipPathLength;

TH1D*     h_EnergyTrk;
TH2D*     h_EnergyTrk_TrueVsReco;
TH2D*     h_EnergyTrk_TrueVsRes;

//TH1D*     h_EnergyShw[nconfigs_SphereR];
//TH2D*     h_EnergyShw_TrueVsReco[nconfigs_SphereR];
//TH2D*     h_EnergyShw_TrueVsRes[nconfigs_SphereR];

TH1D*     h_EnergyTotal[nconfigs_SphereR];
TH2D*     h_EnergyTotal_TrueVsReco[nconfigs_SphereR];
TH2D*     h_EnergyTotal_TrueVsRes[nconfigs_SphereR];

TH2D*     h_EnergyTrk_Frac;
TH2D*     h_EnergyTotal_Frac[nconfigs_SphereR];

// ============================================================
void configure(){

  // check that sphere size vectors are correct 
  if( v_SphereR.size() != nconfigs_SphereR ) {
    std::cout<<"nconfigs_SphereR doesn't match number of sphere radii in vector!\n";
    exit(0);
  }
  
  // open the file and set up the TTree
  TFile* file = new TFile(fFileName.c_str(),"READ");
  fTree = (TTree*)file->Get(fTreeName.c_str());
  setBranches(fTree);
  
  // make output file to store plots
  fOutFile = new TFile(fOutFileName.c_str(),"recreate");
  fOutFile->cd();
  
  // make the histograms
  float resbins_energy    = 12;
  float resbins           = 1600;
  h_BlipEnergy            = new TH1D("BlipEnergy","Blip Energy (pre-smearing);Energy (MeV)",1000,0,10.);
  h_MaxBlipEnergy            = new TH1D("MaxBlipEnergy","Max Blip Energy per Event;",300,0,30.);
  h_BlipPathLength        = new TH1D("BlipPathLength","Blip Path Length;Path Length (cm)",100,0,5.);
  h_EnergyBudget          = new TH1D("EnergyBudget",        ";E_{reco} - E_{0} (MeV)",1000,-25,25);
  h_EnergyTrk             = new TH1D("EnergyTrk",           "Electron Track;Reconstructed Track Energy (MeV)",3000,0,600);
  h_EnergyTrk_TrueVsReco  = new TH2D("EnergyTrk_TrueVsReco","Electron Track;Electron Energy (MeV);Reconstructed Track Energy (MeV)",
                                      120,0,60,120,0,60);
  h_EnergyTrk_TrueVsRes   = new TH2D("EnergyTrk_TrueVsRes","Electron Track;Electron Energy (MeV);Resolution",
                                      resbins_energy,0,60,resbins,-1.,1.);
  h_EnergyTrk_Frac        = new TH2D("EnergyTrk_Frac","Electron Track;Electron Energy (MeV);Mean E_{reco} / E_{true}",resbins_energy,0,60,resbins,0,1);
  h_EnergyTrk_TrueVsReco  ->SetOption("colz");
  h_EnergyTrk_TrueVsRes   ->SetOption("colz");
  h_EnergyTrk_Frac        ->SetOption("colz");
    
  for(int j=0; j<nconfigs_SphereR; j++){
    float R = v_SphereR.at(j);
    //float E = fThreshold*1e3;
    //float s = fSmear*1e3;
    
    //h_EnergyShw[j]    = new TH1D( Form("EnergyShw_%.0fcm",R),
    //                              Form("R = %0.0fcm;Reconstructed Shower Energy (MeV)",R),
    //                              3000,0,600);
    //h_EnergyShw_TrueVsReco[j]    = new TH2D( Form("EnergyShw_%.0fcm_TrueVsReco",R),
    //                              Form("R = %0.0fcm;True Shower-Only Energy (MeV);Reconstructed Shower Energy (MeV)",R),
    //                              120,0,60,120,0,60);
    //h_EnergyShw_TrueVsRes[j]    = new TH2D( Form("EnergyShw_%.0fcm_TrueVsRes",R),
    //                              Form("Shower, R = %0.0fcm;True Shower-Only Energy (MeV);Resolution",R),
    //                              resbins_energy,0,60,resbins,-1.,1.);

    h_EnergyTotal[j]    = new TH1D( Form("EnergyTotal_%.0fcm",R),
                                  Form("R = %0.0fcm;Reconstructed Shower Energy (MeV)",R),
                                  6000,0,600);
    h_EnergyTotal_TrueVsReco[j]    = new TH2D( Form("EnergyTotal_%.0fcm_TrueVsReco",R),
                                  Form("R = %0.0fcm;Electron Energy (MeV);Reconstructed Energy (MeV)",R),
                                  120,0,60,120,0,60);
    h_EnergyTotal_TrueVsRes[j]    = new TH2D( Form("EnergyTotal_%.0fcm_TrueVsRes",R),
                                  Form("Total, R = %0.0fcm;Electron Energy (MeV);Resolution",R),
                                  resbins_energy,0,60,resbins,-1.,1.);
    h_EnergyTotal_Frac[j]    = new TH2D( Form("EnergyTotal_%.0fcm_Frac",R),
                                  Form("Total, R = %0.0fcm;Electron Energy (MeV);Mean E_{reco}/E_{true}",R),
                                  resbins_energy,0,60,resbins,-1.,1.);
    
    h_EnergyTotal_TrueVsRes[j]    ->SetOption("colz");
    h_EnergyTotal_TrueVsReco[j]   ->SetOption("colz");
    h_EnergyTotal_Frac[j]     ->SetOption("colz");
  }
  
  
}


// =============================================================
void ElectronReco_AnaMacro(){

  configure();
  reco();
  makePlots();

  // save all the lower-level histograms
  h_EnergyBudget        ->Write();
  h_BlipEnergy          ->Write();
  h_MaxBlipEnergy       ->Write();
  h_BlipPathLength      ->Write();
  h_EnergyTrk           ->Write();
  h_EnergyTrk_TrueVsReco->Write();
  h_EnergyTrk_TrueVsRes ->Write();
  h_EnergyTrk_Frac->Write();
  for(int id=0; id<nconfigs_SphereR; id++) {
    h_EnergyTotal[id]           ->Write();
    h_EnergyTotal_TrueVsReco[id]->Write();
    h_EnergyTotal_TrueVsRes[id] ->Write();
    //h_EnergyShw[id]             ->Write();
    //h_EnergyShw_TrueVsReco[id]  ->Write();
    //h_EnergyShw_TrueVsRes[id]   ->Write();
    h_EnergyTotal_Frac[id]      ->Write();
  }

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
    
    fTrkEnergyDep = 0.;
    float nuEnergy = _NuEnergy*1e3;

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
      
      // is it an electron/positron?
      if( fabs(pdg) == 11 ){
        
        // if primary electron, record energy dep
        if( proc == "primary" ) {
          fTrueEnergy = KE;
          fTrkEnergyDep   = edep; 
          fSphereCen      = loc;
          if( nuEnergy > 0 ) fTrueEnergy = nuEnergy;
        } 

        // otherwise, record a new blip
        else {
          EnergyDeposit b;
          b.Location = loc;
          b.Energy = edep;
          b.PathLength = dL;
          v_blips.push_back(b);
        }
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

   // if( fabs(total_edep -fTrueEnergy) > 3. ) std::cout<<"!!! ERROR !!!\n";
   // std::cout
   // <<"Total blips      : "<<v_blips.size()<<"\n"
   // <<"Initial energy   : "<<fTrueEnergy<<" MeV\n"
   // <<"Total energy dep : "<<total_edep<<" MeV\n"
   // <<"Difference       : "<<total_edep - fTrueEnergy<<" MeV\n";
    h_EnergyBudget->Fill((total_edep-fTrueEnergy));

    // .................................
    // merging, thresholding, and smearing
    MergeBlips    (v_blips,fMinSep    );

    // plot blip diagnostics
    float maxBlipE = -9.;
    for(size_t i=0; i<v_blips.size(); i++) {
      h_BlipEnergy->Fill(v_blips.at(i).Energy);
      h_BlipPathLength->Fill(v_blips.at(i).PathLength);
      if( v_blips.at(i).Energy > maxBlipE ) maxBlipE = v_blips.at(i).Energy;
    }
    h_MaxBlipEnergy->Fill(maxBlipE);
    
    ThresholdBlips(v_blips,fThreshold );

    // ..................................
    // smearing
    SmearBlips    (v_blips,fSmear     );

    printf("Event %5d -- smear: %6.2f keV, blips found: %d\n",iEvent,fSmear*1e3,(int)v_blips.size());

    //....................................
    // reconstruct the electron shower products
    showerReco(v_blips); 
    
    // fill histograms
    h_EnergyTrk             ->Fill(fTrkEnergyDep);
    h_EnergyTrk_TrueVsReco  ->Fill(fTrueEnergy, fTrkEnergyDep);
    h_EnergyTrk_TrueVsRes   ->Fill(fTrueEnergy,(fTrkEnergyDep-fTrueEnergy)/fTrueEnergy);
    h_EnergyTrk_Frac        ->Fill(fTrueEnergy,fTrkEnergyDep/fTrueEnergy);
    for(size_t iR = 0; iR < v_SphereR.size(); iR++){
      //float true_shower_edep = fTrueEnergy - fTrkEnergyDep;
      //h_EnergyShw[iR]           ->Fill(fShwEnergyDep[iR]);
      //h_EnergyShw_TrueVsReco[iR]->Fill(true_shower_edep,fShwEnergyDep[iR]);
      //h_EnergyShw_TrueVsRes[iR] ->Fill(true_shower_edep,(fShwEnergyDep[iR]-true_shower_edep)/true_shower_edep);
      h_EnergyTotal[iR]           ->Fill(fTotalEnergyDep[iR]);
      h_EnergyTotal_TrueVsReco[iR]->Fill(fTrueEnergy,fTotalEnergyDep[iR]);
      h_EnergyTotal_TrueVsRes[iR] ->Fill(fTrueEnergy,(fTotalEnergyDep[iR]-fTrueEnergy)/fTrueEnergy);
      h_EnergyTotal_Frac[iR]      ->Fill(fTrueEnergy,fTotalEnergyDep[iR]/fTrueEnergy);
    }

  }//>> end loop over events

}



// =============================================================
void showerReco(std::vector<EnergyDeposit>& blips ){
  // Loop through the sphere sizes
  for(size_t iR = 0; iR < v_SphereR.size(); iR++){
    // sphere radius
    float R = v_SphereR.at(iR);
      //fShwEnergyDep[iR]     = 0.;
      fTotalEnergyDep[iR]   = fTrkEnergyDep;
      for(size_t i=0; i<blips.size(); i++){
        EnergyDeposit thisBlip = blips.at(i);
        float d = (thisBlip.Location - fSphereCen).Mag();
        if( d < R || R < 0 ) {
          fTotalEnergyDep[iR] += thisBlip.Energy;
          //fShwEnergyDep[iR]   += thisBlip.Energy;
        }
      }
  }
}



// ============================================================
void makePlots(){
  
  // Make a default canvas
  TCanvas* cdefault = new TCanvas("default","default",200,200);
  cdefault->cd();

  TGraph* gr_trk_frac;
  TGraph* gr_total_frac[nconfigs_SphereR];

  TGraph* gr_trk_sig;
  TGraph* gr_trk_RMS;
  TGraph* gr_shw_sig[nconfigs_SphereR];
  TGraph* gr_shw_RMS[nconfigs_SphereR];
  TGraph* gr_total_sig[nconfigs_SphereR];
  TGraph* gr_total_RMS[nconfigs_SphereR];
  
  // ---------------------------------------------
  std::vector<TH1D*> vec_h1d;
  std::vector<float> vec_E;
  TH2D* hist2D;

  // get RMS histogram for e track
  gr_trk_RMS = new TGraph(); gr_trk_RMS->SetTitle("Trunk-only, RMS");
  gr_trk_sig = new TGraph(); gr_trk_sig->SetTitle("Trunk-only, Peak Fit");
  gr_trk_frac = new TGraph(); gr_trk_frac->SetTitle("Trunk-only");
  hist2D = h_EnergyTrk_TrueVsRes;
  MakeSlices(hist2D,vec_h1d,vec_E,5.,55.);
  FormatTGraph(gr_trk_RMS,kRed,kRed,0,9,0,3);
  FillResolutionGraphs( gr_trk_sig,gr_trk_RMS,
                        vec_h1d, vec_E,
                        "trk",
                        false,false);

  // fractional energy
  hist2D = h_EnergyTrk_Frac;
  MakeSlices(hist2D,vec_h1d,vec_E,5.,55.);
  FormatTGraph(gr_trk_frac,kRed,kRed,0,9,0,3);
  for(size_t i=0; i<vec_h1d.size(); i++){
    gr_trk_frac->SetPoint(gr_trk_frac->GetN(), vec_E.at(i), vec_h1d.at(i)->GetMean());
  }


  std::vector<Color_t> v_coloption= {kGreen+2,kBlue};
  for(size_t iR=0; iR<v_SphereR.size(); iR++){
    float R = v_SphereR.at(iR);
    gr_shw_sig[iR] = new TGraph();    gr_shw_sig[iR]  ->SetTitle(Form("Trunk + Shower, Peak Fit (R = %3.0f cm)",R));
    gr_shw_RMS[iR] = new TGraph();    gr_shw_RMS[iR]  ->SetTitle(Form("Trunk + Shower, RMS (R = %3.0f cm)",R));
    gr_total_sig[iR] = new TGraph();  gr_total_sig[iR]->SetTitle(Form("Trunk + Shower, Peak Fit (R = %3.0f cm)",R));
    gr_total_RMS[iR] = new TGraph();  gr_total_RMS[iR]->SetTitle(Form("Trunk + Shower, RMS (R = %3.0f cm)",R));
    gr_total_frac[iR]= new TGraph();  gr_total_frac[iR]->SetTitle(Form("Trunk + Shower, R = %3.0f cm",R));
    FormatTGraph(gr_total_sig[iR], v_coloption.at(iR), v_coloption.at(iR), 0, 1, 0, 3);
    FormatTGraph(gr_total_RMS[iR], v_coloption.at(iR), v_coloption.at(iR), 0, 9, 0, 3);
    FormatTGraph(gr_shw_sig[iR], v_coloption.at(iR), v_coloption.at(iR), 0, 1, 0, 3);
    FormatTGraph(gr_shw_RMS[iR], v_coloption.at(iR), v_coloption.at(iR), 0, 9, 0, 3);
    FormatTGraph(gr_total_frac[iR], v_coloption.at(iR),v_coloption.at(iR),0,9,0,3);

    hist2D  = h_EnergyTotal_TrueVsRes[iR];
    MakeSlices(hist2D,vec_h1d,vec_E,5.,55.);
    std::cout
    <<"-----------------------------------------\n"
    <<"Beginning fits for "<<hist2D->GetTitle()<<"\n";
    FillResolutionGraphs( gr_total_sig[iR],gr_total_RMS[iR],
                          vec_h1d, vec_E,
                          Form("total_%3.0fcm",R),
                          true, true);
    
    //hist2D  = h_EnergyShw_TrueVsRes[iR];
    //MakeSlices(hist2D,vec_h1d,vec_E,5.,30.);
    //std::cout
    //<<"-----------------------------------------\n"
    //<<"Beginning fits for "<<hist2D->GetTitle()<<"\n";
    //FillResolutionGraphs( gr_shw_sig[iR],gr_shw_RMS[iR],
    //                      vec_h1d, vec_E,
    //                      Form("shower_%3.0fcm",R),
    //                      true, false);
  
    // fractional energy
    hist2D = h_EnergyTotal_Frac[iR];
    MakeSlices(hist2D,vec_h1d,vec_E,5.,55.);
    for(size_t i=0; i<vec_h1d.size(); i++){
      gr_total_frac[iR]->SetPoint(gr_total_frac[iR]->GetN(), vec_E.at(i), vec_h1d.at(i)->GetMean());
    }



    //cdefault->cd();
  }
  
  TCanvas* cfrac = new TCanvas("cfrac","cfrac",700,450);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.05);
  auto mgf = new TMultiGraph();
  mgf->Add(gr_trk_frac,     "APL");
  mgf->Add(gr_total_frac[0], "APL");
  mgf->Add(gr_total_frac[1], "APL");
  mgf->Draw("a");
  //mg->SetTitle("Total Energy");
  FormatAxes(mgf,0.045,0.045,1,1);
  mgf->GetXaxis()->SetTitle("True Electron Energy (MeV)");
  mgf->GetYaxis()->SetTitle("Mean Reconstructed Energy Fraction");
  mgf->GetXaxis()->CenterTitle();
  mgf->GetYaxis()->CenterTitle();
  cfrac->Update();
 
  TCanvas* c = new TCanvas("c","c",700,450);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.05);
  auto mg = new TMultiGraph();
  mg->Add(gr_trk_RMS,      "APL");
  mg->Add(gr_total_RMS[0], "APL");
  mg->Add(gr_total_RMS[1], "APL");
  mg->Add(gr_total_sig[0], "APL");
  mg->Add(gr_total_sig[1], "APL");
  mg->Draw("a");
  FormatAxes(mg,0.045,0.045,1,1);
  //mg->SetTitle("Total Energy");
  mg->GetXaxis()->SetTitle("True Electron Energy (MeV)");
  mg->GetYaxis()->SetTitle("Fractional Resolution");
  mg->GetXaxis()->CenterTitle();
  mg->GetYaxis()->CenterTitle();
  c->Update();

  /*
  TCanvas* c2 = new TCanvas("c2","c2",700,500);
  auto mg2 = new TMultiGraph();
  mg2->Add(gr_shw_sig[0], "APL");
  mg2->Add(gr_shw_sig[1], "APL");
  mg2->Add(gr_shw_RMS[0], "APL");
  mg2->Add(gr_shw_RMS[1], "APL");
  mg2->Draw("a");
  mg2->SetTitle("Shower-Only Energy");
  mg2->GetXaxis()->SetTitle("Shower-Only Energy (MeV)");
  mg2->GetYaxis()->SetTitle("Fractional Resolution");
  c2->Update();
  */

  delete cdefault;

}

// ####################################################################################
// Helper function for the resolution slice methods below. This function takes in a 
// 2D histogram and slices it up into some number of 1D projections.
void MakeSlices( const TH2D* h2, 
                std::vector<TH1D*>& vec_h1, 
                std::vector<float>& vec_x,
                float x_min, float x_max)
{
  vec_h1.clear();
  vec_x.clear(); 
  float         bwidth    = h2->GetXaxis()->GetBinWidth(1);
  const size_t  nbins     = h2->GetXaxis()->GetNbins();
  for(size_t i=1; i<nbins; i++){
    float x       = h2->GetXaxis()->GetBinCenter(i);
    float lowEdge = h2->GetXaxis()->GetBinLowEdge(i);
    float highEdge= lowEdge + bwidth;
    if( lowEdge >= x_min && highEdge <= x_max ){
      TH1D* h = h2->ProjectionY(Form("yproj_%6.2f_MeV",x),i,i);
      h->SetTitle(h2->GetTitle());
      vec_h1.push_back(h);
      vec_x.push_back(x);
    }
  }
}

//#######################################################################################
void FillResolutionGraphs( TGraph* gr_sig, TGraph* gr_rms, 
                          std::vector<TH1D*>& vec_h1d,
                          std::vector<float>& vec_x,
                          std::string name,
                          bool verbose = false, 
                          bool plotSliceFits = false )
{
    TF1 f_gaus("f_gaus","gaus(0)");
    TCanvas* c = new TCanvas("tmp","default",200,200);
    c->cd();
    TCanvas* c_slice;
    for(size_t i=0; i<vec_h1d.size(); i++){
      float E = vec_x.at(i);
      TH1D* h = vec_h1d.at(i);
      float rms = h->GetRMS();
      float peakE = h->GetXaxis()->GetBinCenter(h->GetMaximumBin());
      float max   = h->GetBinContent(h->GetMaximumBin());
      float mu, sig;
      if(verbose) std::cout<<"Slice E= "<<E<<" MeV ---\n";
      // 1st fit attempt
      f_gaus.SetParameter(0,max);
      f_gaus.SetParLimits(0,max*0.9,max*1.1);
      f_gaus.SetParameter(1,peakE);
      //f_gaus.SetParLimits(1,peakE-fabs(peakE*0.1),peakE+fabs(peakE*0.1));
      f_gaus.SetParLimits(1,-0.5,0.1);
      f_gaus.SetParameter(2,0.008);
      f_gaus.SetParLimits(2,0.0005,0.2);
      h->Fit("f_gaus","BQN WW");
      mu  = f_gaus.GetParameter(1);
      sig = f_gaus.GetParameter(2);
      std::cout<<"  1st pass fit --> mean= "<<mu<<"  sig= "<<sig<<"\n";
      // 2nd attempt
      f_gaus.SetRange(mu-5.*sig,mu+5*sig);
      h->Fit("f_gaus","BQRN WW");
      mu  = f_gaus.GetParameter(1);
      sig = f_gaus.GetParameter(2);
      std::cout<<"  2nd pass fit --> mean= "<<mu<<"  sig= "<<sig<<"\n";
      // 3rd attempt
      f_gaus.SetRange(mu-2*sig,mu+2*sig);
      h->Fit("f_gaus","BQR WW");
      mu  = f_gaus.GetParameter(1);
      sig = f_gaus.GetParameter(2);
      //// 4th attempt
      //f_gaus.SetRange(mu-1.5*sig,mu+1.5*sig);
      //h->Fit("f_gaus","BQR");
      // if chi^2 still bad, try one more time
      if( plotSliceFits ){
        c_slice = new TCanvas(Form("slice_%s_%3.1fMeV",name.c_str(),E),Form("slice_%s_%3.1fMeV",name.c_str(),E),500,500);
        c_slice->cd();
        gStyle->SetOptFit(1);
        h->DrawCopy();
        c->cd();
      }
      mu  = f_gaus.GetParameter(1);
      sig = f_gaus.GetParameter(2);
      float chi2_v = f_gaus.GetChisquare()/f_gaus.GetNDF();
      std::cout<<"  final fit --> mean= "<<mu<<"  sig= "<<sig<<"   Chi2/NDF = "<<f_gaus.GetChisquare()/f_gaus.GetNDF()<<"\n";
      //if( chi2_v < 500. )
      AddPointToGraph(gr_sig,E,sig);
      AddPointToGraph(gr_rms,E,rms);
    }
    delete c;
}

