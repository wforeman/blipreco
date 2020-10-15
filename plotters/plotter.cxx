#include "core/blip_nu.h"

TH1D* h_cc_0to30_neutron = new TH1D("cc_0to30_neutron","",25,0,25);
TH1D* h_cc_30to60_neutron = new TH1D("cc_30to60_neutron","",25,0,25);
TH1D* h_cc_0to30_noneutron = new TH1D("cc_0to30_noneutron","",25,0,25);
TH1D* h_cc_30to60_noneutron = new TH1D("cc_30to60_noneutron","",25,0,25);
TH1D* h_es_0to30 = new TH1D("es_0to30","",25,0,25);
TH1D* h_es_30to60 = new TH1D("es_30to60","",25,0,25);

TH1D* h_ratio_cc_neutron = new TH1D("h_ratio_cc_neutron","",25,0,25);
TH1D* h_ratio_cc_noneutron = new TH1D("h_ratio_cc_noneutron","",25,0,25);

void plotter(){

  TFile* file_cc = new TFile("plots_cc.root","READ");
  TFile* file_es = new TFile("plots_es.root","READ");

  h_cc_0to30_neutron  = (TH1D*)file_cc->Get("Mult_0to30cm_Neutron"); 
  h_cc_30to60_neutron = (TH1D*)file_cc->Get("Mult_30to60cm_Neutron"); 
  h_cc_0to30_noneutron  = (TH1D*)file_cc->Get("Mult_0to30cm_noNeutron"); 
  h_cc_30to60_noneutron = (TH1D*)file_cc->Get("Mult_30to60cm_noNeutron"); 
  h_es_0to30  = (TH1D*)file_es->Get("Mult_0to30cm_noNeutron"); 
  h_es_30to60 = (TH1D*)file_es->Get("Mult_30to60cm_noNeutron"); 

  h_ratio_cc_neutron  = (TH1D*)file_cc->Get("Ratio_Neutron");
  h_ratio_cc_noneutron  = (TH1D*)file_cc->Get("Ratio_noNeutron");
  
  h_cc_0to30_neutron    ->Scale(1./h_cc_0to30_neutron->GetEntries());
  h_cc_30to60_neutron   ->Scale(1./h_cc_30to60_neutron->GetEntries());
  h_cc_0to30_noneutron  ->Scale(1./h_cc_0to30_noneutron->GetEntries());
  h_cc_30to60_noneutron ->Scale(1./h_cc_30to60_noneutron->GetEntries());
  h_es_0to30            ->Scale(1./h_es_0to30->GetEntries());
  h_es_30to60           ->Scale(1./h_es_30to60->GetEntries());

  h_ratio_cc_neutron    ->Scale(1./h_ratio_cc_neutron->GetEntries());
  h_ratio_cc_noneutron    ->Scale(1./h_ratio_cc_noneutron->GetEntries());

  h_cc_0to30_neutron  ->SetLineColor(kRed);
  h_cc_0to30_neutron  ->SetLineStyle(9);
  h_cc_0to30_noneutron  ->SetLineColor(kBlue);
  h_cc_0to30_noneutron  ->SetLineStyle(1);
  h_cc_30to60_neutron  ->SetLineColor(kRed);
  h_cc_30to60_neutron  ->SetLineStyle(9);
  h_cc_30to60_noneutron  ->SetLineColor(kBlue);
  h_cc_30to60_noneutron  ->SetLineStyle(1);
  h_es_0to30              ->SetLineColor(kGreen+2);
  h_es_0to30              ->SetLineStyle(2);
  h_es_30to60             ->SetLineColor(kGreen+2);
  h_es_30to60             ->SetLineStyle(2);
  
  h_cc_0to30_neutron  ->SetLineWidth(2);
  h_cc_0to30_noneutron  ->SetLineWidth(2);
  h_cc_30to60_neutron  ->SetLineWidth(2);
  h_cc_30to60_noneutron  ->SetLineWidth(2);
  h_es_0to30  ->SetLineWidth(2);
  h_es_0to30  ->SetLineWidth(2);

  gStyle->SetOptStat(0);

  TCanvas* c1 = new TCanvas("c1","c1",600,450);
//  h_cc_0to30_neutron->GetYaxis()->SetRangeUser(0,1000);
  h_cc_0to30_neutron->Draw("hist");
  h_cc_0to30_noneutron->Draw("same hist");
  h_es_0to30->Draw("same hist");
  
  TCanvas* c2 = new TCanvas("c2","c2",600,450);
  //h_cc_30to60_neutron->GetYaxis()->SetRangeUser(0,0.5);
  h_cc_30to60_neutron->Draw("hist");
  h_cc_30to60_noneutron->Draw("same hist");
  //h_es_30to60->Draw("same hist");

  
  TCanvas* c3 = new TCanvas("c3","c3",600,450);
  //h_cc_30to60_neutron->GetYaxis()->SetRangeUser(0,0.5);
  h_ratio_cc_neutron->Draw("hist");
  h_ratio_cc_noneutron->Draw("same hist");

}
