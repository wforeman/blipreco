#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom2.h"
#include "TMath.h"
#include <TMinuit.h>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TFitter.h"
#include "TVector3.h"
#include "TStyle.h"
#include "TLine.h"
#include <time.h>
#include <vector>
#include <math.h>

bool IsParticleDescendedFrom(int,int);
bool IsParticleDescendedFrom(int,int,bool);

// Random engine
TRandom2 *fRand = new TRandom2();

// Data structure for "blip" object
struct EnergyDeposit {
  TVector3  Location;
  float     Energy;
  float     PathLength;
  bool      isGrouped;
};

// AnaTree variables
const int kMax = 10000;
double     _NuEnergy;
int       _event;
int       _no_hits_stored;
int       _geant_list_size;
int			  _pdg[kMax];
double     _Eng[kMax];
//double     _EndE[kMax];
double     _StartT[kMax];
double     _EndT[kMax];
double     _StartPointx[kMax];
double     _StartPointy[kMax];
double     _StartPointz[kMax];
double     _EndPointx[kMax];
double     _EndPointy[kMax];
double     _EndPointz[kMax];
double     _pathlen[kMax];
int       _NumberDaughters[kMax];
int       _Mother[kMax];
int       _TrackId[kMax];
int       _process_primary[kMax];
std::vector<std::string> *_processname = 0;

void setBranches(TTree *tree){
  TBranch *br = 0;
  tree->SetBranchAddress("enu_truth",&_NuEnergy);
  tree->SetBranchAddress("event",&_event);
  tree->SetBranchAddress("geant_list_size",&_geant_list_size);
  tree->SetBranchAddress("pdg",&_pdg);
  tree->SetBranchAddress("Eng",&_Eng);
  //tree->SetBranchAddress("EndE",&_EndE);
  //tree->SetBranchAddress("StartT",&_StartT);
  //tree->SetBranchAddress("EndT",&_EndT);
  tree->SetBranchAddress("StartPointx",&_StartPointx);
  tree->SetBranchAddress("StartPointy",&_StartPointy);
  tree->SetBranchAddress("StartPointz",&_StartPointz);
  tree->SetBranchAddress("EndPointx",&_EndPointx);
  tree->SetBranchAddress("EndPointy",&_EndPointy);
  tree->SetBranchAddress("EndPointz",&_EndPointz);
  //tree->SetBranchAddress("pathlen",&_pathlen);
  tree->SetBranchAddress("NumberDaughters",&_NumberDaughters);
  tree->SetBranchAddress("Mother",&_Mother);
  tree->SetBranchAddress("TrackId",&_TrackId);
  tree->SetBranchAddress("G4Process",&_processname,&br);
  //tree->SetBranchAddress("processname",&_processname,&br);
}

//==============================================================
float CalcEnergyDepParticle(int iP){
  
  // if photon, no deposited energy
  if( _pdg[iP] == 22 ) return 0.;
  
  // at first approximation, Edep is the kinetic energy
  float Edep = 1e3*_Eng[iP] - 0.511;// - _EndE[iP]);

  //std::cout<<"  Edep0 = "<<Edep<<"\n";

  // if there are daughters, we may need to subtract their energy
  if( _NumberDaughters[iP] > 0 ) {
    for(size_t i=iP+1; i<_geant_list_size; i++){
      if( _Mother[i] == _TrackId[iP] ) {
       
        float E = 1e3*_Eng[i];

        if( _processname->at(i) == "eBrem" ) 
          Edep -= E;

        if( _processname->at(i) == "eIoni" ) 
          Edep -= (E - 0.511);

        // If we're dealing with positron annihilation, things are more 
        // complicated. If the annihilation happens at rest, the energy going
        // into the 2gamma pair is purely from mass energy so doesn't detract
        // from the blip. So we only subtract off the kinetic energy of 
        // annihilation photons that exceeds the electron mass energy.
        if( _processname->at(i) == "annihil" ) 
          Edep -= (E - 0.511);

      }
    }
  }
  return Edep;
}

//===========================================================
float CalcEnergyDep(int iP){
  
  // If this is an electron that came from another electron, 
  // it would have already been grouped as part of the
  // contiguous "blip" previously, so don't count it.
  if( _pdg[iP] == 11 && _processname->at(iP) == "eIoni" ) return 0.;

  // First calculate energy depoosited *directly* by this particle
  //std::cout<<"   Calculating Edep for particle "<<iP<<"... pdg = "<<_pdg[iP]<<"\n";
  float Edep = CalcEnergyDepParticle(iP);

  // We want to loop through any contiguous electrons (produced
  // with process "eIoni") and add the energy they deposit to this
  // energy deposition.
  if( _NumberDaughters[iP] > 0 ){
    for(size_t i=iP+1; i<_geant_list_size; i++){
      if( _processname->at(i) == "eIoni" ) {
        bool breakAtPhots=true;
        if( IsParticleDescendedFrom(_TrackId[i],_TrackId[iP],breakAtPhots) ){
          Edep += CalcEnergyDepParticle(i);
        }
      }
    }
  }

  return Edep;

}

//============================================================
// Boolean function which loops through the lineaage of a particle
bool IsParticleDescendedFrom( int particleID, int motherID, bool breakAtPhotons ){
  if( particleID < 0 || motherID < 0 ) return false;
  if( motherID == 0 || motherID == particleID ) return true;
  bool isDescended = false;
  std::vector<int> lineage;
  lineage.reserve(200);
//  lineage.push_back(motherID);
  for(size_t p=0; p<_geant_list_size; p++){
    // add the mother to the lineage once we come to it
    // (check if it's a photon too)
    if( _TrackId[p] == motherID ) {
      if( breakAtPhotons && _pdg[p] == 22 ) return false;
      lineage.push_back(_TrackId[p]); 
    }
    // update lineage
    for(size_t i=0; i<lineage.size(); i++){
      if( _Mother[p] == lineage[i] && (!breakAtPhotons || _pdg[p] != 22 ) ) {
        lineage.push_back(_TrackId[p]);
        if( _TrackId[p] == particleID ) isDescended = true;
      }
    }
    if( _TrackId[p] == particleID ) return isDescended;
  }
  return false;
}

bool IsParticleDescendedFrom( int particleID, int motherID ){
  return IsParticleDescendedFrom(particleID,motherID,false);
}

//=============================================================
void SmearBlips(std::vector<EnergyDeposit>&v, float sig){
  if( sig <= 0 ) return;
  for(int i=0; i<v.size(); i++) {
    float e = v.at(i).Energy;
    bool flag = true;
    while (flag){
      float es = e + fRand->Gaus(0.,sig);
      if( es > 0. ) {
        v.at(i).Energy = es;
        flag = false;
      }
    }
  }
}

//==============================================================
void ThresholdBlips(std::vector<EnergyDeposit> &v, float thresh){
  if(thresh <= 0 ) return;
  std::vector<EnergyDeposit> v_out;
  for(int i=0; i<v.size(); i++) {
    if( v.at(i).Energy > thresh ) v_out.push_back(v.at(i));
  }
  v.clear();
  v = v_out;
}

//===============================================================
void MergeBlips(std::vector<EnergyDeposit> &v, float sep){
  if( sep <= 0 ) return;
  std::vector<EnergyDeposit> v_out;
  for(int i=0; i<v.size(); i++) v.at(i).isGrouped = false;
  for(int i=0; i<v.size(); i++) {
    if( v.at(i).isGrouped ) continue;
    else v.at(i).isGrouped = true;
    EnergyDeposit newBlip;
    newBlip.Energy = v.at(i).Energy;
    newBlip.Location = v.at(i).Location;
    newBlip.PathLength = v.at(i).PathLength;
    for(int j=i+1; j<v.size(); j++){
      if( v.at(j).isGrouped ) continue;
      float d = (v.at(i).Location - v.at(j).Location).Mag();
      if( d < sep ) {
        newBlip.Energy += v.at(j).Energy;
        newBlip.PathLength += v.at(j).PathLength;
        v.at(j).isGrouped = true;
      }
    }
    v_out.push_back(newBlip);
  }
  v.clear();
  v = v_out;
  //return v_out;
}
