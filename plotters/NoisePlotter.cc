//////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////// 
#include "blip.h"
#include "tools.h"

void NoisePlotter(){
  
  TGraph* gr_75 = new TGraph();
  TGraph* gr_150 = new TGraph();

  gr_75->SetPoint(gr_75->GetN(), 10., 8.3);
  gr_75->SetPoint(gr_75->GetN(), 30., 9.3);
  gr_75->SetPoint(gr_75->GetN(), 50., 10.9);
  gr_75->SetPoint(gr_75->GetN(), 70., 13.2);
  gr_150->SetPoint(gr_150->GetN(), 10., 13.2);
  gr_150->SetPoint(gr_150->GetN(), 30., 13.9);
  gr_150->SetPoint(gr_150->GetN(), 50., 15.6);
  gr_150->SetPoint(gr_150->GetN(), 70., 16.3);

  
  TCanvas* cfrac = new TCanvas("cfrac","cfrac",700,450);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.05);
  auto mg = new TMultiGraph();
  mg->Add(gr_75,     "APL");
  mg->Add(gr_150,     "APL");
  mg->Draw("a");
  //mg->SetTitle("Total Energy");
  FormatAxes(mg,0.045,0.045,1,1);
  mg->GetXaxis()->SetTitle("Per-Blip Smearing [keV]");
  mg->GetYaxis()->SetTitle("FWHM Energy Resolution (#sigma = FWHM/2.355)");
  mg->GetXaxis()->CenterTitle();
  mg->GetYaxis()->CenterTitle();
  cfrac->Update();

  /////////////////////////////////////////////////////////
  TCanvas* c2 = new TCanvas("c2","c2",700,450);
  TGraph* gr = new TGraph();
  /*
  gr->SetPoint(gr->GetN(), 4.,  2.54);
  gr->SetPoint(gr->GetN(), 5.,  3.08);
  gr->SetPoint(gr->GetN(), 7.,  4.63);
  gr->SetPoint(gr->GetN(), 9.,  6.34);
  gr->SetPoint(gr->GetN(), 10., 7.16);
  gr->SetPoint(gr->GetN(), 12., 7.87);
  //gr->SetPoint(gr->GetN(), 14., 9.09);
  //gr->SetPoint(gr->GetN(), 16., 11.1);
  */
  gr->SetPoint(gr->GetN(), 4.,  2.2);
  gr->SetPoint(gr->GetN(), 5.,  2.85);
  gr->SetPoint(gr->GetN(), 7.,  4.34);
  gr->SetPoint(gr->GetN(), 9.,  6.10);
  gr->SetPoint(gr->GetN(), 10., 6.86);
  gr->SetPoint(gr->GetN(), 12., 7.71);
  gr->SetMarkerSize(1);
  gr->SetMarkerStyle(20);
  gr->Draw("AP");

}
