//////////////////////////////////////////////////////////////////
// 
//  Ar39 Analysis Macro
//
//  Must have ROOT data file containing analysistree/anatree in 
//  same directory.
//
//  To run:
//  > root -l Ar39_AnaMarco.cc
//
////////////////////////////////////////////////////////////////// 
#include "core/blip.h"
#include "core/tools.h"

// ============================================================
// ===================   Parameters ===========================
#define nconfigs_EThresh 2
#define nconfigs_SphereR 8
std::string             fileName    = "../mcfiles/AnaTree_ar39_DUNE.root";
std::string             treeName    = "analysistree/anatree";
std::string             outFile     = "plots.root";
std::vector<float>      v_SphereR   = { 10., 20., 40., 60., 80., 100., 125., 150. }; // cm
std::vector<float>      v_EThresh   = { 75., 300. }; // keV
std::vector<Color_t>    v_EThreshCol= { kRed,kBlue};
float                   fMinSep     = 0.2; 
bool                    fiducialCut = true;
int                     nMax        = 100; // # blips selected to draw spheres around 
                                      // (if < 0, uses all available)

// volume dimensions
float xlim[2] = {-390.,390.};
float ylim[2] = {-660.,660.};
float zlim[2] = {-60., 1460.};

//===============================================================

void AnaMacro_Ar39();
void configure();
bool isInActiveVol(TVector3 point, float margin);
TVector3 randomPointInVol(float margin);
void loopTheBlips(std::vector<EnergyDeposit> blips, int n );
void makePlots();

TFile*                fOutFile;
TTree*                fTree;
  
TH1D*     h_Edep_PerElectron; 
TH1D*     h_BlipE;
TH1D*     h_BlipE_thresh[nconfigs_EThresh];
TH1D*     h_Dist;           
TH1D*     h_DistClosest[nconfigs_EThresh];
TH1D*     h_NBlips[nconfigs_EThresh][nconfigs_SphereR];
TH1D*     h_Edep[nconfigs_EThresh][nconfigs_SphereR];

void configure(){
  
  // check that things are set up correctly
  if( v_SphereR.size() != nconfigs_SphereR || v_EThresh.size() != nconfigs_EThresh ) {
    printf("ERROR! Check that definition of nconfigs_XXX in header matches associated vectors\n");
    return;
  }


  // make the histograms
  h_Edep_PerElectron  = new TH1D("Edep_PerElectron","Energy deposited in LAr per electron;Deposited energy [MeV]",100,0.,1.); 
  h_Dist              = new TH1D("Dist","Distance between blips;Distance [cm]",100,0.,2000.);
  h_BlipE             = new TH1D("BlipE","Blip energies;^{39}Ar Blip Energy (MeV);Number of Blips",100,0,5.);
  for(int i=0; i<nconfigs_EThresh; i++) {
    float E = v_EThresh.at(i);
    h_DistClosest[i] = new TH1D(Form("DistClosest_%.0fkeV",E),Form("Closest blip, E > %.0f keV;Distance [cm]",E), 100,0.,100.);
    h_BlipE_thresh[i] = new TH1D(Form("BlipE_%.0fkeV",E),Form("Blip energies (E > %.0f keV);^{39}Ar Blip Energy (MeV);Number of Blips",E),100,0,5.);
    for(int j=0; j<nconfigs_SphereR; j++){
      float R = v_SphereR.at(j);
      h_NBlips[i][j]  = new TH1D(Form("NBlips_%.0fkeV_%.0fcm",E,R),Form("Number blips within R = %.0f cm, E > %.0f keV",R,E),100,0,100);
      h_Edep[i][j]    = new TH1D(Form("Edep_%.0fkeV_%.0fcm",E,R),Form("Total deposited energy within R=%.0f cm, E > %.0f keV;Deposited energy [MeV]",R,E),600,0.,30);
    }
  }
  
  // open the file and set up the TTree
  TFile* file = new TFile(fileName.c_str(),"READ");
  fTree = (TTree*)file->Get(treeName.c_str());
  setBranches(fTree);

}


// =============================================================
void AnaMacro_Ar39(){

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
      
      // make TVector of position
      TVector3 loc(_StartPointx[i],_StartPointy[i],_StartPointz[i]);
   
      // if electron, make a blip
      if( fabs(_pdg[i]) == 11 ){
        float edep = CalcEnergyDep(i);
        h_Edep_PerElectron->Fill(edep);
        EnergyDeposit b;
        b.Location = loc;
        b.Energy = edep;
        b.Time = _StartT[i];
        v_blips.push_back(b);
        //if( edep > 0.010 ) v_blips.push_back(b);
        // "golden" blips are those that are far enough
        // away from the edges to avoid edge effects when sphere-drawing
        //if( isInActiveVol(loc,v_SphereR.back()) ) {
        //  if( edep > 0.075 ) { nGolden++; v_blips_golden.push_back(b); }
        //  if( edep > 0.300 ) { nGolden_300++; }
        //}
      }

    

    }//>> end particle loop
   
    // ------------------------------------------------
    // Merge blips (SKIPPING THIS FOR NOW, SINCE IT TAKES
    // WAY TOO LONG ON ~9000 AR39 BLIPS!)
    //MergeBlips    (v_blips,fMinSep    );
    
    
    // fill blip energy histogram
    for(size_t j=0;j<v_blips.size();j++) h_BlipE->Fill(v_blips.at(j).Energy);

    // -----------------------------------------------
    // Make different thresholded subsets
    for(int ie=0; ie<nconfigs_EThresh; ie++){
      std::vector<EnergyDeposit> v_blips_2 = v_blips;
      ThresholdBlips(v_blips_2,1e-3*v_EThresh[ie]);
      for(size_t j=0;j<v_blips_2.size();j++){
        h_BlipE_thresh[ie]->Fill(v_blips_2.at(j).Energy);
      }
    }
    
    // -----------------------------------------------
    // Run the sphere-clustering routine, which cycles
    // through different sphere sizes and blip thresholds
    loopTheBlips( v_blips, nMax );

    printf("Event %5i -- number blips: %lu\n",iEvent,v_blips.size());

  }//>> end loop over events

  // ------------------------------------------------- 
  // make output file to store plots
  fOutFile = new TFile(outFile.c_str(),"recreate");
  fOutFile->cd();
  
  // make Ebias and resolution TGraphs
  makePlots();

  // save all the lower-level histograms
  h_Edep_PerElectron->Write();
  h_Dist              ->Write();
  h_BlipE     ->Write();
  for(int ie=0; ie<nconfigs_EThresh; ie++) h_DistClosest[ie]->Write();
  for(int ie=0; ie<nconfigs_EThresh; ie++) h_BlipE_thresh[ie]->Write();
  for(int ie=0; ie<nconfigs_EThresh; ie++){ 
    for(int id=0; id<nconfigs_SphereR; id++) h_NBlips[ie][id]->Write();
    for(int id=0; id<nconfigs_SphereR; id++) h_Edep[ie][id]->Write();
  }



  fOutFile->Close();
  return;
}



// ============================================================
void makePlots(){

  TGraph* gr_Ebias[nconfigs_EThresh];
  TGraph* gr_Erms[nconfigs_EThresh];
  TGraph* gr_Nblips[nconfigs_EThresh];
  TGraph* gr_NblipsRMS[nconfigs_EThresh];
  for(int i=0; i<nconfigs_EThresh; i++){
    gr_Ebias[i] = new TGraph();
    gr_Erms[i] = new TGraph();
    gr_Nblips[i] = new TGraph();
    gr_NblipsRMS[i] = new TGraph();
    FormatTGraph(gr_Ebias[i],v_EThreshCol[i],v_EThreshCol[i],20,1,0.7,2);
    FormatTGraph(gr_Erms[i],v_EThreshCol[i],v_EThreshCol[i],20,1,0.7,2);
    FormatTGraph(gr_Nblips[i],v_EThreshCol[i],v_EThreshCol[i],20,1,0.7,2);
    FormatTGraph(gr_NblipsRMS[i],v_EThreshCol[i],v_EThreshCol[i],20,1,0.7,2);
    for(int j=0; j<nconfigs_SphereR; j++){
      gr_Ebias[i]->SetPoint(gr_Ebias[i]->GetN(),v_SphereR.at(j),h_Edep[i][j]->GetMean());
      gr_Erms[i]->SetPoint(gr_Erms[i]->GetN(),v_SphereR.at(j),h_Edep[i][j]->GetRMS());
      gr_Nblips[i]->SetPoint(gr_Nblips[i]->GetN(),v_SphereR.at(j),h_NBlips[i][j]->GetMean());
      gr_NblipsRMS[i]->SetPoint(gr_NblipsRMS[i]->GetN(),v_SphereR.at(j),h_NBlips[i][j]->GetRMS());
    }
  }

  // canvas dimensions (x, y)
  int wx = 600;
  int wy = 400;

  TCanvas* c_Ebias = new TCanvas("Ebias","Energy Bias",wx,wy);
  gPad->SetGrid();
  gPad->SetLogy();
  auto mg = new TMultiGraph();
  for(int i=0; i<nconfigs_EThresh; i++ )
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
  
  TCanvas* c_NBlipsRMS = new TCanvas("NBlipsRMS","RMS of Number of Blips",wx,wy);
  gPad->SetGrid();
  gPad->SetLogy();
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

}

// =============================================================
void loopTheBlips(std::vector<EnergyDeposit> blips, int nmax ){

  // Choose a random point in the simulated volume, ensuring it is
  // at least the largest possible sphere-radius away from edge.
  for(int i=0; i<nmax; i++){
    
    TVector3 sphereCenter = randomPointInVol( v_SphereR.back() );

    float distClosest = 999999.;
    int   N[nconfigs_SphereR];
    float distClosestThresh[nconfigs_EThresh]={0};
    int   N_Thresh[nconfigs_EThresh][nconfigs_SphereR]={{0}};
    float Esum[nconfigs_EThresh][nconfigs_SphereR]={{0.}};
    
    // loop over ALL blips
    for(int j=0; j<blips.size(); j++){
      if( j == i ) continue;
      bool flag_D[nconfigs_SphereR]={false};
      bool flag_E[nconfigs_EThresh]={false};
      float d = (blips.at(j).Location - sphereCenter).Mag();
      float E = blips.at(j).Energy; // MeV
      h_Dist->Fill(d);
      if( d < distClosest ) distClosest = d;
 
      // flag this blip for: (a) in sphere, and (b) over threshold  
      for(int id=0; id<nconfigs_SphereR; id++){
        if( d < v_SphereR.at(id) ) { flag_D[id] = true; N[id]++; }
          for(int ie=0; ie<nconfigs_EThresh; ie++){
            if( id==0 ) {
              if( E > v_EThresh.at(ie)*1e-3 ) { flag_E[ie] = true; }
              if( flag_E[ie] && (d < distClosestThresh[ie] || distClosestThresh[ie] == 0) ) 
              distClosestThresh[ie]=d;
            }
          if( flag_E[ie] && flag_D[id] ) N_Thresh[ie][id]++;
        }
      }

    
      // now add up the energy
      for(int ie=0; ie<nconfigs_EThresh; ie++){
        for(int id=0; id<nconfigs_SphereR; id++){
          if( flag_D[id] && flag_E[ie] ) Esum[ie][id] += E;       
        }
      }
  

    }//<< endloop over *other* blips
  
    // fill histograms
    for(int ie=0; ie<nconfigs_EThresh; ie++){
      h_DistClosest[ie]->Fill(distClosestThresh[ie]);
      for(int id=0; id<nconfigs_SphereR; id++){
        h_NBlips[ie][id]->Fill(N_Thresh[ie][id]);
        h_Edep[ie][id]->Fill(Esum[ie][id]);
      }
    }

  }//<< endloop over *target* blip

}



//==============================================================
bool isInActiveVol(TVector3 p, float margin){
  if(     p.X() > (xlim[0]+margin) && p.X() < (xlim[1]-margin)
      &&  p.Y() > (ylim[0]+margin) && p.Y() < (ylim[1]-margin)
      &&  p.Z() > (zlim[0]+margin) && p.Z() < (zlim[1]-margin) ){
    return true; } else { return false; }
}

TVector3 randomPointInVol(float margin){
  TVector3 pt(0,0,0);
  float xmin = xlim[0]+margin;
  float xmax = xlim[1]-margin;
  float ymin = ylim[0]+margin;
  float ymax = ylim[1]-margin;
  float zmin = zlim[0]+margin;
  float zmax = zlim[1]-margin;
  pt.SetX(xmin + fRand->Rndm()*(xmax-xmin));
  pt.SetY(ymin + fRand->Rndm()*(ymax-ymin));
  pt.SetZ(zmin + fRand->Rndm()*(zmax-zmin));
  return pt;
}
