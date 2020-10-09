//////////////////////////////////////////////////////////////////
// 
//  BlipReco
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
std::string             fFileName     = "../mcfiles/anatree_gamma1_larworld.root";
//std::string             fFileName    = "../mcfiles/anatree_gamma2_larworld.root";
//std::string             fFileName    = "../mcfiles/anatree_neutrons_1eV.root";
//std::string             fFileName    = "../mcfiles/anatree_neutron_10MeV_larworld.root";
std::string             fTreeName     = "analysistree/anatree";
std::string             fOutFileName  = "plots.root";
float                   fThreshold    = 0.0750; // 75 keV
float                   fSmear        = 0.00; // 50 keV
float                   fMinSep       = 0.20; // cm
std::vector<float>      v_SphereR   = { 20., 30., 60. }; // cm

bool printOut = false;

// initial guesses for Gaus fit
float fMean = 5.6;
float fSigma = 0.3;

//===============================================================

void AnaMacro_Gamma();
void configure();
void loopTheBlips(std::vector<EnergyDeposit> blips );
void makePlots();
void doFits(TH1D*);

TFile*    fOutFile;
TTree*    fTree;

TH1D*     h_ElEnergy_fromEl;
TH1D*     h_ElEnergy_fromPh;
TH1D*     h_ElEnergy_fromBrem;
TH1D*     h_BremEnergy;

TH1D*     h_BlipEnergy;
TH1D*     h_BlipMult;
TH1D*     h_BlipSeparation;
TH1D*     h_BlipMultSphere[nconfigs_SphereR];
TH1D*     h_Energy_zoom[nconfigs_SphereR];
TH1D*     h_Energy[nconfigs_SphereR];

void configure(){
  
  // make the histograms
  
  int   bins = 1000;
  float emax = 2.; // MeV
  h_ElEnergy_fromEl   = new TH1D("ElEnergy_fromEl",   "Electrons produced by an electron;Kinetic Energy (keV);Entries",     bins,0.0,emax*1e3);
  h_ElEnergy_fromPh   = new TH1D("ElEnergy_fromPh",   "Electrons produced by a photon;Kinetic Energy (keV);Entries",        bins,0.0,emax*1e3);
  h_ElEnergy_fromBrem = new TH1D("ElEnergy_fromBrem", "Electrons produced by a bremm. photon;Kinetic Energy (keV);Entries", bins,0.0,emax*1e3);
  h_BremEnergy        = new TH1D("BremEnergy",        "Bremm. photons;Kinetic Energy (keV);Entries",                        bins,0.0,emax*1e3);
  h_BlipEnergy      = new TH1D("BlipEnergy","Individual Blip Energies;Individual Blip Energy (MeV);Entries",1000,0.0,10.0);
  h_BlipMult        = new TH1D("BlipMult","Blip Multiplicity;Number of Blips",100,0,100);
  h_BlipSeparation  = new TH1D("BlipSeparation","Inter-Blip Separation;Separation Distance [cm]",2000,0,100);

  for(int j=0; j<nconfigs_SphereR; j++){
    float R = v_SphereR.at(j);
    float E = fThreshold*1e3;
    float s = fSmear*1e3;

    h_Energy_zoom[j]  = new TH1D(Form("Energy_%.0f_zoom",R),
                                 Form("R = %0.f cm;Energy (MeV);Entries",R),
                                 320,0.1,3.3);
    h_Energy[j]       = new TH1D(Form("Energy_%.0f",R),
                                 Form("R = %0.f cm;Energy (MeV);Entries",R),
                                 200,0.,20);
    h_BlipMultSphere[j] = new TH1D(Form("BlipMult_%.0f",R),
                                 Form("R = %0.f cm;Number of Blips;Entries",R),
                                 100,0,100);
  }
  
  // open the file and set up the TTree
  TFile* file = new TFile(fFileName.c_str(),"READ");
  fTree = (TTree*)file->Get(fTreeName.c_str());
  setBranches(fTree);
  
}


// =============================================================
void AnaMacro_Gamma(){
  
  // check that sphere size vectors are correct 
  if( v_SphereR.size() != nconfigs_SphereR ) {
    std::cout<<"nconfigs_SphereR doesn't match number of sphere radii in vector!\n";
    exit(0);
  }

  // configure histograms and TFile
  configure();
  
  // --------------------------------------------
  // loop over the events
  for(int iEvent=0; iEvent<fTree->GetEntries(); iEvent++){
    fTree->GetEntry(iEvent);
   
    // Make a vector of Blip objects to fill
    std::vector<EnergyDeposit> v_blips;
    float total_edep = 0.;
    float primary_energy = 0.;

    // Loop over particles and fill in the blips.
    int nParticles = _geant_list_size;
    for(int i=0; i<nParticles; i++){

      // ignore any particles outside TPC active
      //if( _inTPCActive[i] == 0 ) continue;
    
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
      total_edep += edep;

      if( proc == "primary" ) {
        primary_energy += KE;
      }

      // if bremm photon, record energy
      if( pdg == 22 && _processname->at(i) == "eBrem" ) h_BremEnergy->Fill(KE*1e3);

      // if electron, record the energy
      if( fabs(pdg) == 11 ) {
        int motherPDG = PdgOfMother(trackId);
        float KE_keV = KE*1e3;
        if( fabs(motherPDG) == 11 ) h_ElEnergy_fromEl->Fill(KE_keV);
        if( fabs(motherPDG) == 22 ) {
          h_ElEnergy_fromPh->Fill(KE_keV);
          if( ProcessOfMother(trackId) == "eBrem" )
            h_ElEnergy_fromBrem->Fill(KE_keV);
        }
      }

      // if electron, make blip
      if( fabs(pdg) == 11 ) {
        EnergyDeposit b;
        b.Energy = edep;
        b.Location = loc;
        b.isGrouped = false;
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
          nD);
      */

    }//>> end particle loop
    
    //if( (total_edep - primary_energy) > 0.01 ) std::cout<<"!!! ERROR !!!\n";
    //std::cout
    //<<"Total blips      : "<<v_blips.size()<<"\n"
    //<<"Initial energy   : "<<primary_energy<<" MeV\n"
    //<<"Total energy dep : "<<total_edep<<" MeV\n";
    
    // .................................
    // merge close-together blips
    MergeBlips(v_blips,fMinSep);
    
    //float flagE0 = 3.08;
    //float flagE1 = 3.09;
    for(int i=0; i<v_blips.size(); i++){
      float E = v_blips.at(i).Energy;
      h_BlipEnergy->Fill(E);
      //if( E >= flagE0 && E <= flagE1 ) {
      //  std::cout<<"BLIP ENERGY FLAG ... E = "<<v_blips.at(i).Energy<<" MeV\n";
      //}
    }

    // .................................
    // do thresholding
    ThresholdBlips(v_blips,fThreshold);
    h_BlipMult->Fill(v_blips.size());

    // ................................
    // do smearing
    SmearBlips(v_blips,fSmear);

    printf("Event %5d -- smear: %6.2f keV, blips found: %d\n",iEvent,fSmear*1e3,(int)v_blips.size());

    // first look at all the inter-blip separations
    for(int i=0; i<v_blips.size(); i++){
      for(int j=0; j<v_blips.size(); j++){
        if( j <= i ) continue;
        float d = (v_blips.at(i).Location - v_blips.at(j).Location).Mag();
        h_BlipSeparation->Fill(d);
      }
    }

    // now that we've saved all the interesting blips, we will perform
    // the iterative sphere-drawing method where for gamma reconstruction
    //   (a) Identify highest-E blip
    //   (b) Draw sphere centered on that blip, and
    //       add up energies of all blips inside sphere
    //   (3) Locate highest-E blip of those that were not 
    //       grouped into the first sphere, and repeat,
    //       only grouping previously un-grouped blips.
    loopTheBlips( v_blips );
  
  }//>> end loop over events

  // ------------------------------------------------- 
  // make output file to store plots
  fOutFile = new TFile(fOutFileName.c_str(),"recreate");
  fOutFile->cd();
  
  // make Ebias and resolution TGraphs
  makePlots();

  //doFits(h_EnergySmeared[1]);

  // save all the lower-level histograms
  h_BremEnergy->Write();
  h_ElEnergy_fromEl->Write();
  h_ElEnergy_fromPh->Write();
  h_ElEnergy_fromBrem->Write();
  h_BlipEnergy->Write();
  h_BlipMult->Write();
  h_BlipSeparation->Write();
  for(int id=0; id<nconfigs_SphereR; id++) {
    h_Energy_zoom[id]->Write();
    h_Energy[id]->Write();
    h_BlipMultSphere[id]->Write();
  }

  fOutFile->Close();
  return;
}


// =============================================================
void loopTheBlips(std::vector<EnergyDeposit> blips ){
  
  // Get primary location (to serve as sphere center)
  int primPDG = _pdg[0];
  TVector3 primLoc(_StartPointx[0],_StartPointy[0],_StartPointz[0]);

  /*
  // -----------------------------------------
  // TEMPORARY simple smearing
  //std::cout<<"using simple 1-sphere clsutering, ceneterd at primary start-point\n";
  for(size_t iR = 0; iR < v_SphereR.size(); iR++){
    // sphere radius
    float R = v_SphereR.at(iR);
    float sumE          = 0.;
    for(size_t i=0; i<blips.size(); i++){
      EnergyDeposit thisBlip = blips.at(i);
      float d = (thisBlip.Location - primLoc).Mag();
      if( d < R ) {
        blips.at(i).isGrouped = true;
        sumE += thisBlip.Energy;
//        if( R > 999. ) std::cout<<"  adding blip at distance "<<d<<" of energy "<<thisBlip.Energy<<"\n";
      }
    }
  //  if( R > 999 ) std::cout<<"  sumE total = "<<sumE<<"\n";
    h_Energy[iR]      ->Fill(sumE);
    h_Energy_zoom[iR] ->Fill(sumE);
    //h_BlipMultSphere[iR]->Fill(nblips_sphere);
  }
  // ------------------------------------------
  */
  
  // Now do sphere grouping...
  // Loop through the sphere sizes
  for(size_t iR = 0; iR < v_SphereR.size(); iR++){
   
    // sphere radius
    float R = v_SphereR.at(iR);

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
        h_Energy[iR]      ->Fill(sumE);
        h_Energy_zoom[iR] ->Fill(sumE);
        h_BlipMultSphere[iR]->Fill(nblips_sphere);
      //}

    }//<< end while(ungrouped blips)

  }//<< end loop over sphere radii

}

void doFits(TH1D* h){
  
  /*
  // par0 = N
  // par1 = mu
  // par2 = sigma

  std::cout<<h->GetTitle()<<"\n";

  TF1 fit("fit","[0]*exp(-0.5*pow((x-[1])/[2],2))");
  fit.SetParameter(0,h->GetMaximum());
  fit.SetParameter(1,fMean);
  fit.SetParameter(2,fSigma);
  h->Fit("fit","");
  */

}

// ============================================================
void makePlots(){
  
  // ---------------------------------------------------------
  // canvas dimensions (x, y)
  int wx = 600;
  int wy = 400;

  std::vector<float>      v_SphereR   = { 20., 30., 60. }; // cm
 
  TCanvas* ce = new TCanvas( "BlipEnergy","BlipEnergy",wx,wy);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.12);
  FormatAxes(h_BlipEnergy,0.045,0.045,1.1,1.4);
  h_BlipEnergy->SetLineColor(kBlack);
  h_BlipEnergy->SetLineWidth(2);
  h_BlipEnergy  ->SetTitle("");
  h_BlipEnergy  ->GetXaxis()->SetTitle("Individual Blip Energy (MeV)");
  h_BlipEnergy  ->GetYaxis()->SetTitle("Entries");
  h_BlipEnergy  ->GetXaxis()->CenterTitle();
  h_BlipEnergy  ->GetYaxis()->CenterTitle();
  h_BlipEnergy  ->Draw();
  gStyle->SetOptStat(0);
  
  /*
  TCanvas* cc = new TCanvas( "Energy60","Energy60",wx,wy);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.12);
  FormatAxes(h_EnergySmeared_SingleGamma,0.045,0.045,1.1,1.4);
  h_EnergySmeared_SingleGamma->SetLineColor(kBlack);
  h_EnergySmeared_SingleGamma->SetLineWidth(2);
  h_EnergySmeared_SingleGamma->SetFillColor(kCyan-8);
  h_EnergySmeared_SingleGamma  ->SetTitle("");
  h_EnergySmeared_SingleGamma  ->GetXaxis()->SetTitle("Energy (MeV)");
  h_EnergySmeared_SingleGamma  ->GetYaxis()->SetTitle("Entries");
  h_EnergySmeared_SingleGamma  ->GetXaxis()->CenterTitle();
  h_EnergySmeared_SingleGamma  ->GetYaxis()->CenterTitle();
  h_EnergySmeared_SingleGamma  ->Draw();
  gStyle->SetOptStat(0);
  
  TCanvas* cc2 = new TCanvas( "Energy30","Energy30",wx,wy);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.12);
  FormatAxes(h_EnergySmeared_SingleGamma_30,0.045,0.045,1.1,1.4);
  h_EnergySmeared_SingleGamma_30->SetLineColor(kBlack);
  h_EnergySmeared_SingleGamma_30->SetLineWidth(2);
  h_EnergySmeared_SingleGamma_30->SetFillColor(kCyan-8);
  h_EnergySmeared_SingleGamma_30  ->SetTitle("");
  h_EnergySmeared_SingleGamma_30  ->GetXaxis()->SetTitle("Energy (MeV)");
  h_EnergySmeared_SingleGamma_30  ->GetYaxis()->SetTitle("Entries");
  h_EnergySmeared_SingleGamma_30  ->GetXaxis()->CenterTitle();
  h_EnergySmeared_SingleGamma_30  ->GetYaxis()->CenterTitle();
  h_EnergySmeared_SingleGamma_30  ->Draw();
  gStyle->SetOptStat(0);
*/

  for(size_t i=0; i<v_SphereR.size(); i++){
  
    float R = v_SphereR.at(i);

    TCanvas* cz = new TCanvas( Form("cz_%.0f",R),Form("Energy_%.0f_zoom",R),wx,wy);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.12);
    FormatAxes(h_Energy_zoom[i],0.045,0.045,1.1,1.4);
    h_Energy_zoom[i]->SetFillColor(kCyan-8);
    h_Energy_zoom[i]->SetLineColor(kBlack);
    h_Energy_zoom[i]->SetLineWidth(2);
    h_Energy_zoom[i]->SetTitle(Form("Energy (R = %.0f cm)",R));
    h_Energy_zoom[i]->GetXaxis()->SetTitle("Energy (MeV)");
    h_Energy_zoom[i]->GetYaxis()->SetTitle("Entries");
    h_Energy_zoom[i]->GetXaxis()->CenterTitle();
    h_Energy_zoom[i]->GetYaxis()->CenterTitle();
    h_Energy_zoom[i]->Draw();
    gStyle->SetOptStat(0);
    
    TCanvas* c = new TCanvas( Form("c_%.0f",R),Form("Energy_%.0f",R),wx,wy);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.12);
    FormatAxes(h_Energy[i],0.045,0.045,1.1,1.4);
    h_Energy[i]->SetFillColor(kCyan-8);
    h_Energy[i]->SetLineColor(kBlack);
    h_Energy[i]->SetLineWidth(2);
    h_Energy[i]->SetTitle(Form("Energy (R = %.0f cm)",R));
    h_Energy[i]->GetXaxis()->SetTitle("Energy (MeV)");
    h_Energy[i]->GetYaxis()->SetTitle("Entries");
    h_Energy[i]->GetXaxis()->CenterTitle();
    h_Energy[i]->GetYaxis()->CenterTitle();
    h_Energy[i]->Draw();
    gStyle->SetOptStat(0);

    // ---------------------------------------------------
    // Fitting
    TH1D* h = h_Energy_zoom[i];
    
    // find FWHM
    float max = -999.;
    float imax = -9;
    int bin1 = h->FindBin(1.3);
    int bin2 = h->FindBin(1.5);
    for(int i=bin1; i<bin2; i++){
      float v = h->GetBinContent(i);
      if( v > max ) {
        max = v;
        imax = i;
      }
    }
    float Epeak = h->GetXaxis()->GetBinCenter(imax);
    float x1, x2;
    for(int i=imax; i<h->GetXaxis()->GetNbins(); i++){
      if( h->GetBinContent(i) <= max/2. ){
        x2 = h->GetXaxis()->GetBinCenter(i);
        break;
      }
    }
    for(int i=imax; i>0; i--){
      float v = h->GetBinContent(i);
      if( h->GetBinContent(i) <= max/2. ){
        x1 = h->GetXaxis()->GetBinCenter(i);
        break;
      }
    }
    float FWHM = x2-x1;
    float sig = FWHM/2.355;
    float mid = (x1+x2)/2;
    std::cout
    <<"R = "<<R<<" --> FWHM= "<<FWHM<<" MeV,   Ep= "<<Epeak<<" MeV,   midpt= "<<mid<<" MeV,   sig= "<<sig<<" MeV,  sig/Ep = "<<100.*sig/Epeak<<" %, sig/mid = "<<100.*sig/mid<<" % \n";


  }


  /*
  TCanvas* c_Energy = new TCanvas("Energy","Energy",wx,wy);
  h_Energy[ 
    mg->Add(gr_Ebias[i], "APL");
  mg->Draw("a");  
  FormatAxes(mg, 0.045, 0.045, 1.1, 1);
  mg->GetXaxis()->SetTitle("Sphere Radius (cm)");
  mg->GetYaxis()->SetTitle("Mean Energy Bias (MeV)");
  mg->GetXaxis()->CenterTitle();
  mg->GetYaxis()->CenterTitle();
  mg->GetYaxis()->SetRangeUser(0.001,30);

  TCanvas* c_Erms = new TCanvas("Erms","Energy RMS",wx,wy);
  gPad->SetGrid();
  gPad->SetLogy();
  auto mg2 = new TMultiGraph();
  for(int i=0; i<nconfigs_EThresh; i++ )
    mg2->Add(gr_Erms[i], "APL");
  mg2->Draw("a");  
  FormatAxes(mg2, 0.045, 0.045, 1.1, 1);
  mg2->GetXaxis()->SetTitle("Sphere Radius (cm)");
  mg2->GetYaxis()->SetTitle("Energy RMS (MeV)");
  mg2->GetXaxis()->CenterTitle();
  mg2->GetYaxis()->CenterTitle();
  mg2->GetYaxis()->SetRangeUser(0.03,6);
  
  TCanvas* c_NBlips = new TCanvas("NBlips","Number of Blips",wx,wy);
  gPad->SetGrid();
  gPad->SetLogy();
  auto mg3 = new TMultiGraph();
  for(int i=0; i<nconfigs_EThresh; i++ )
    mg3->Add(gr_Nblips[i], "APL");
  mg3->Draw("a");  
  FormatAxes(mg3, 0.045, 0.045, 1.1, 1);
  mg3->GetXaxis()->SetTitle("Sphere Radius (cm)");
  mg3->GetYaxis()->SetTitle("Mean Blip Multiplicity");
  mg3->GetXaxis()->CenterTitle();
  mg3->GetYaxis()->CenterTitle();
  mg3->GetYaxis()->SetRangeUser(0.001,100);
  
  TCanvas* c_NBlipsRMS = new TCanvas("NBlipsRMS","RMS of Number of Blips",600,600);
  gPad->SetGrid();
  auto mg4 = new TMultiGraph();
  for(int i=0; i<nconfigs_EThresh; i++ )
    mg4->Add(gr_NblipsRMS[i], "APL");
  mg4->Draw("a");  
  mg4->GetXaxis()->SetTitle("Sphere Radius [cm]");
  mg4->GetYaxis()->SetTitle("RMS Number of Blips");
  mg4->GetXaxis()->SetTitleOffset(1.1);
  mg4->GetYaxis()->SetTitleOffset(1.3);

  //mg->Draw();
  fOutFile->cd();
  c_Ebias->Write();
  c_Erms->Write();
  c_NBlips->Write();
  c_NBlipsRMS->Write();
  */

}


