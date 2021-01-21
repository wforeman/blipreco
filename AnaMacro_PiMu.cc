//////////////////////////////////////////////////////////////////
// 
//  To run:
//  > root -l <macro-name>.cc
//
//////////////////////////////////////////////////////////////////
#include "core/vars.h" 
#include "core/blip.h"

// ===================   Parameters ===========================
//std::string             fFileName     = "AnaTree_piMinus_1MeV.root";
std::string             fFileName     = "AnaTree_muMinus_1MeV.root";
float                   fThreshold    = 0.075; //0.075; //> Blip threshold (75 keV default)
float                   fSmear        = 0.00;   //> Smearing (50 keV usually)
float                   fMinSep       = 0.20;   //> Min blip separation used (0 = no merging)
float                   fSphereR      = 60.;    //> cm
//===============================================================

void AnaMacro_PiMu();
void configure();
void process();
void makeplots();

std::string fTreeName     = "analysistree/anatree";
std::string fOutFileName  = "output_plots.root";
std::string fMCFile       = "../mcfiles/"+fFileName;         

int       ntotal[2]         = {};
int       nsel_4MeV[2]       = {};
int       nsel_7[2]          = {};
int       nsel_7_4MeV[2]     = {};
int       nsel_14[2]         = {};
int       nsel_8MeV[2]       = {};
int       nsel_14_8MeV[2]    = {};
int       nsel_vtx5MeV[2]    = {};
int       nsel_7_4MeV_vtx5MeV[2] = {};
int       nsel_14_8MeV_vtx5MeV[2] = {};

TFile*    fOutFile;
TTree*    fTree;
TH1D*     h_BlipE;
TH1D*     h_BlipMult;
TH1D*     h_BlipMult_Sphere;
TH1D*     h_EventType;
TH1D*     h_VertexE;
TH1D*     h_NumProtons;
TH1D*     h_NumNeutrons;
TH2D*     h_Mult_vs_Energy_capture;
TH2D*     h_Mult_vs_Energy_decay;

// =============================================================
void AnaMacro_PiMu() {
  
  // configure histograms and TFile
  configure();
  
  // open the file and set up the TTree
  TFile* file = new TFile(fMCFile.c_str(),"READ");
  fTree = (TTree*)file->Get(fTreeName.c_str());
  setBranches(fTree);

  // run the main analyzer
  process();
  
  // make output file to store plots
  fOutFile = new TFile(fOutFileName.c_str(),"recreate");
  fOutFile->cd();

  // make plots
  makeplots();

  // save histograms
  h_BlipE           ->Write();
  h_EventType       ->Write();
  h_BlipMult        ->Write();
  h_BlipMult_Sphere ->Write();
  h_VertexE         ->Write();
  h_NumProtons      ->Write();
  h_NumNeutrons      ->Write();
  h_Mult_vs_Energy_capture->Write();
  h_Mult_vs_Energy_decay->Write();
  fOutFile->Close();
  
  return;

}

//================================================================
void configure(){
  float Emax  = 50;
  int   Ebins = 50;
  float Nmax  = 60; 
  h_EventType               = new TH1D("EventType","Capture = 0, Decay = 1",2,0,2);
  h_BlipE                   = new TH1D("BlipEnergy","Individual Blip Energies;Individual Blip Energy (MeV);Entries",1000,0.0,10.0);
  h_BlipMult                = new TH1D("BlipMult","Blip Multiplicity;Number of Blips Per Event",300,0,300);
  h_BlipMult_Sphere         = new TH1D("BlipMult_Sphere",Form("Blip Multiplicity (R < %3.0f cm);Number of Blips Per Event",fSphereR),100,0,100);
  h_VertexE                 = new TH1D("VertexEnergy","Vertex energy (protons > 5 MeV);Energy (MeV);Number of Events",300,0.,150.);
  h_NumProtons              = new TH1D("NumProtons","Number of protons emitted from capture;Number of Protons;Number of Events",8,0,8);
  h_NumNeutrons              = new TH1D("NumNeutrons","Number of neutrons emitted from capture;Number of Neutrons;Number of Events",8,0,8);
  h_Mult_vs_Energy_capture  = new TH2D("Mult_vs_Energy_capture","Capture;Blip Multiplicity;Summed Blip Energy (MeV)",Nmax,0,Nmax,Ebins,0,Emax);
  h_Mult_vs_Energy_decay    = new TH2D("Mult_vs_Energy_decay","Decay;Blip Multiplicity;Summed Blip Energy (MeV)", Nmax,0,Nmax,Ebins,0,Emax);
  h_Mult_vs_Energy_capture  ->SetOption("COLZ");
  h_Mult_vs_Energy_decay    ->SetOption("COLZ");
}

// =============================================================
void process() {
  
  // loop over the events
  for(int iEvent=0; iEvent<fTree->GetEntries(); iEvent++){
    fTree->GetEntry(iEvent);
    //if( iEvent >= 5000 ) break;

    // event info 
    bool      isDecay     = false;
    bool      isCapture   = false;
    int       primaryId   = 0;
    int       primaryPdg  = -9999;
    int       primarySign = 0;
    float     primaryTf   = -9999;
    float     primaryKEf  = -9999;
    TVector3  primaryEnd;
    float     vertexEnergy= 0;
    int       nProtons    = 0;
    int       nNeutrons   = 0;

    // vector of Blip objects to fill
    std::vector<EnergyDeposit> v_blips;
      
    // Loop over particles and fill in the blips.
    for(int i=0; i<_geant_list_size; i++){

      TVector3 loc(_StartPointx[i],_StartPointy[i],_StartPointz[i]);
      TVector3 locEnd(_EndPointx[i],_EndPointy[i],_EndPointz[i]);
      std::string proc = _processname->at(i);
      int pdg       = _pdg[i];
      int trackId   = _TrackId[i]; 
      int mother    = _Mother[i];
      int nD        = _NumberDaughters[i];
      float T0      = _StartT[i];
      float Tf      = _EndT[i];
      float E       = 1e3*_Eng[i];
      float mass    = 1e3*Mass(i);
      float KEf     = 1e3*_EndE[i] - mass;
      float dL      = (loc-locEnd).Mag();
      bool isCAR = ( proc.find("CaptureAtRest") != std::string::npos );

      // calculate the energy deposited by this particle
      // (includes depositions from contiguous daughters
      // like "eIoni" electrons)   
      float edep = CalcEnergyDep(i);

      // make blip if electron or proton with edep < 3 MeV
      if(     proc != "primary" && edep > 0 
          &&  ( abs(pdg) == 11 || (abs(pdg) == 2212 && edep < 3. )) ){
          EnergyDeposit b;
          b.Energy    = edep;
          b.Location  = loc;
          b.PathLength= dL;
          b.isGrouped = false;
          v_blips.push_back(b);
      }
     
      // .....................................................................
      // If primary particle is pion, end processes will be:
      //  - Decay
      //  - pi+/pi-Inelastic
      //  - hBertiniCaptureAtRest
      //
      //  If muon, then end processes are:
      //  - Decay
      //  - muMinusCaptureAtRest
      //
      //  For pions, the CaptureAtRest process implies the pion was absorbed
      //  by the nucleus, but this is not the case for muons. A mu- can have
      //  a CaptureAtRest end process but still undergo decay before it is 
      //  finally absorbed. So we must check for decay products.
      // .....................................................................
      if( proc == "primary" ){
        primaryId     = trackId;
        primaryPdg    = pdg;
        primaryTf     = Tf;
        primarySign   = -1*pdg/abs(pdg);
        primaryKEf    = KEf;
        primaryEnd    = locEnd;
      }


      // --------------------------------------------------------------
      // special case for mu-: look for either decay electron or capture products
      if( abs(primaryPdg)  == 13 && mother == primaryId ) {
        if( !isCapture && !isDecay && isCAR )                   { isCapture = true; isDecay = false; }
        if( !isDecay && fabs(pdg) == 11 && (T0-primaryTf > 1) ) { isCapture = false; isDecay = true; }
      }//< end mu- death type
      
      // --------------------------------------------------------------
      // special case for pi-: look for either decay muon or capture products
      if( abs(primaryPdg)  == 211 && mother == primaryId ) {
        if( !isCapture && !isDecay && isCAR ) { isCapture = true; isDecay = false; }
        if( !isDecay && fabs(pdg) == 13 )     { isCapture = false; isDecay = true; }
      }//< end pi- death type

      // ---------------------------------------------------------------
      // keep track of vertex energy (protons > 5 MeV)
      if( abs(pdg) == 2212 && mother == primaryId && edep > 5.0 ) 
        vertexEnergy += edep;

      // ---------------------------------------------------------------
      // if capture, count outgoing protons/neutrons
      if( isCapture && mother == primaryId && proc.find("CaptureAtRest") != std::string::npos ){
        if( pdg == 2212 ) nProtons++;
        if( pdg == 2112 ) nNeutrons++;
      }

      // (for debugging output -- keep commented out during normal running)
      if(0){ //&& i < 100){
        float d = (loc-primaryEnd).Mag();
        printf("  %3i PDG: %10i, dL=%8.2f,  E=%8.3f,  Edep=%8.3f,  T0=%7.2f, Tf=%7.2f, moth=%3i,  %12s,  Ndaught=%3i,  blips=%3lu\n",
          trackId,
          pdg,
          //d,
          dL,
          E-mass,
          edep,
          T0,
          Tf,
          mother,
          proc.c_str(),
          nD,
          v_blips.size());
      }

      

    }//>> end particle loop
    
    int eventType = -1;
    if( isCapture ) {eventType = 0; }
    if( isDecay   ) {eventType = 1; }
    h_EventType->Fill(eventType);

    if ( isCapture ) {
      h_NumProtons->Fill(nProtons);
      h_NumNeutrons->Fill(nNeutrons);
    }

//    if( !isCapture && !isDecay ) std::cout<<"FLAG\n";
//    std::cout<<"Decay "<<isDecay<<"    Capture "<<isCapture<<"\n";

    // .................................
    // merge close-together blips
    MergeBlips(v_blips,fMinSep);
   
    // .................................
    // do thresholding
    ThresholdBlips(v_blips,fThreshold);
    
    h_BlipMult->Fill(v_blips.size());
    h_VertexE->Fill(vertexEnergy);
  
    // ===========================================
    // Sphere drawing stuff
    //
    float totalE        = 0.;
    float totalE_sphere = 0.; //elEnergyDep;
    int   nBlips            = v_blips.size();
    int   nBlips_sphere     = 0;
    for(int i=0; i<nBlips; i++){
      float E = v_blips.at(i).Energy;
      h_BlipE  ->Fill(E);
      totalE    += E;
      TVector3 d    = (v_blips.at(i).Location-primaryEnd);
      if( d.Mag() > 0.5 ) {
        if( d.Mag() < fSphereR ) {
          nBlips_sphere++;
          totalE_sphere += E;
        }
      }
    }
    //float aveBlipE = 0.;
    //if( nBlips_sphere ) aveBlipE = totalBlipE_sphere / nBlips_sphere;

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
    
    h_BlipMult_Sphere->Fill(nBlips_sphere);
    if( isDecay   ) h_Mult_vs_Energy_decay  ->Fill(nBlips_sphere,totalE_sphere);
    if( isCapture ) h_Mult_vs_Energy_capture->Fill(nBlips_sphere,totalE_sphere);

    // selection cuts
    float nn = nBlips_sphere;
    float ee = totalE_sphere;
    float vt = vertexEnergy;
    //std::cout<<"eventType "<<eventType<<"   nn "<<nBlips_sphere<<"   ee "<<totalE_sphere<<"   vt "<<vertexEnergy<<"\n";
    if( eventType >= 0 ) {
    ntotal[eventType]++;
    if( nn > 7 )              nsel_7[eventType]++;
    if( ee >= 4. )            nsel_4MeV[eventType]++;
    if( nn > 7 && ee >= 4.)   nsel_7_4MeV[eventType]++;
    if( nn > 14 )             nsel_14[eventType]++;
    if( ee >= 8. )            nsel_8MeV[eventType]++;
    if( nn > 14 && ee >= 8.)  nsel_14_8MeV[eventType]++;
    if( vt > 5. ){            nsel_vtx5MeV[eventType]++;
      if( nn > 7 && ee >=4)   nsel_7_4MeV_vtx5MeV[eventType]++;
      if( nn > 14 && ee >=8)  nsel_14_8MeV_vtx5MeV[eventType]++; }
    }
   
    //if( isCapture && totalE_sphere > 15. ) std::cout<<"FLAG\n";
    printf("MC file: %s, event %5d, blips found: %d, total E: %5.1f MeV (%5.1f MeV in sphere)\n",fFileName.c_str(), iEvent,(int)v_blips.size(),totalE,totalE_sphere);

  }//>> end loop over events


  std::cout<<"totals "<<ntotal[0]<<"    "<<ntotal[1]<<"\n";

  // ------------------------------------------------
  // Channel ID counts
  for(size_t i=0; i<2; i++){
  if( ntotal[i] == 0 ) continue;
  if( i==0 ) std::cout<<"\n\nCAPTURE\n";
  if( i==1 ) std::cout<<"\n\nDECAY\n";
  std::cout
  <<"======================================================\n"
  <<" N> 7      E> -----                  Nsel = "<<float(nsel_7[i])/ntotal[i]<<"\n"
  <<" N> --     E> 4 MeV                  Nsel = "<<float(nsel_4MeV[i])/ntotal[i]<<"\n"
  <<" N> 7      E> 4 MeV                  Nsel = "<<float(nsel_7_4MeV[i])/ntotal[i]<<"\n"
  <<" -----------------------------------------------------\n"
  <<" N> 14     E> -----                  Nsel = "<<float(nsel_14[i])/ntotal[i]<<"\n"
  <<" N> --     E> 8 MeV                  Nsel = "<<float(nsel_8MeV[i])/ntotal[i]<<"\n"
  <<" N> 14     E> 8 MeV                  Nsel = "<<float(nsel_14_8MeV[i])/ntotal[i]<<"\n"
  <<" -----------------------------------------------------\n"
  <<" N> --     E> -----  E_vtx> 5 MeV    Nsel = "<<float(nsel_vtx5MeV[i])/ntotal[i]<<"\n"
  <<" N> 7      E> 4 MeV  E_vtx> 5 MeV    Nsel = "<<float(nsel_7_4MeV_vtx5MeV[i])/ntotal[i]<<"\n"
  <<" N> 14     E> 8 MeV  E_vtx> 5 MeV    Nsel = "<<float(nsel_14_8MeV_vtx5MeV[i])/ntotal[i]<<"\n"
  <<"======================================================\n";
  }

  /*
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
  */


}

void makeplots(){
  
  TCanvas* c1 = new TCanvas("c1","c1",600,450);
  gStyle->SetOptStat(0);
  gPad->SetMargin(0.10, 0.12, 0.12, 0.08);
  h_Mult_vs_Energy_capture->GetXaxis()->CenterTitle();
  h_Mult_vs_Energy_capture->GetYaxis()->CenterTitle();
  h_Mult_vs_Energy_capture->GetXaxis()->SetLabelSize(0.045);
  h_Mult_vs_Energy_capture->GetXaxis()->SetTitleSize(0.045);
  h_Mult_vs_Energy_capture->GetYaxis()->SetLabelSize(0.045);
  h_Mult_vs_Energy_capture->GetYaxis()->SetTitleSize(0.045);
  h_Mult_vs_Energy_capture->GetZaxis()->SetLabelSize(0.045);
  h_Mult_vs_Energy_capture->GetZaxis()->SetTitleSize(0.045);
  h_Mult_vs_Energy_capture->Draw();
  c1->Write();

  TCanvas* c2 = new TCanvas("c2","c2",600,450);
  gStyle->SetOptStat(0);
  gPad->SetMargin(0.10, 0.12, 0.12, 0.08);
  h_Mult_vs_Energy_decay->GetXaxis()->CenterTitle();
  h_Mult_vs_Energy_decay->GetYaxis()->CenterTitle();
  h_Mult_vs_Energy_decay->GetXaxis()->SetLabelSize(0.045);
  h_Mult_vs_Energy_decay->GetXaxis()->SetTitleSize(0.045);
  h_Mult_vs_Energy_decay->GetYaxis()->SetLabelSize(0.045);
  h_Mult_vs_Energy_decay->GetYaxis()->SetTitleSize(0.045);
  h_Mult_vs_Energy_decay->GetZaxis()->SetLabelSize(0.045);
  h_Mult_vs_Energy_decay->GetZaxis()->SetTitleSize(0.045);
  h_Mult_vs_Energy_decay->Draw();
  c2->Write();

}
