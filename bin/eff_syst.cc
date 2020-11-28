#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TEfficiency.h>
#include <TGraphErrors.h>
#include <iostream>

using namespace std;

#define particle 1 //0 = B+;   1 = Bs;

int eff_syst(){
  TFile* f_raw_yield = particle ?  new TFile("~/work2/BinQGP/results/Bs/Bpt/pT.root") : new TFile("~/work2/BinQGP/results/Bu/Bpt/pT.root");
  TGraphAsymmErrors* raw_yield = (TGraphAsymmErrors*)f_raw_yield->Get("Graph;1");
  
  TFile* f_eff0 = particle ? new TFile("~/work2/BinQGP/results/Bs/efficiency/root_files/efficiency0.root") : new TFile("~/work2/BinQGP/results/Bu/efficiency/root_files/efficiency0.root");
  TFile* f_eff1 = particle ? new TFile("~/work2/BinQGP/results/Bs/efficiency/root_files/efficiency1.root") : new TFile("~/work2/BinQGP/results/Bu/efficiency/root_files/efficiency1.root");

  int n_pt_bins = raw_yield->GetN();
  double *x_values = raw_yield->GetX();
 
  double pt_bins[] = {5, 10, 15, 20, 50};

  TEfficiency* efficiency0 = new TEfficiency("efficiency0", "efficiency0", n_pt_bins, pt_bins);
  TEfficiency* efficiency1 = new TEfficiency("efficiency1", "efficiency1", n_pt_bins, pt_bins);
  efficiency0 = (TEfficiency*)f_eff0->Get("hist_tot_noweights_clone");
  efficiency1 = (TEfficiency*)f_eff1->Get("hist_tot_weights_clone");

  double eff0;
  double eff1;
  double syst;

  double y_values[n_pt_bins];
  double y_errors[n_pt_bins];
  double x_error_low[n_pt_bins];
  double x_error_high[n_pt_bins];

  for(int i = 0; i < n_pt_bins; i++)
    {
      eff0 = efficiency0->GetEfficiency(i + 1);
      eff1 = efficiency1->GetEfficiency(i + 1);
      syst = (eff1 - eff0) / eff0;
      y_values[i] = syst;
      y_errors[i] = 0;
      x_error_low[i] = raw_yield->GetErrorXlow(i);
      x_error_high[i] = raw_yield->GetErrorXhigh(i);
    }

  double pt_zero[n_pt_bins];
  for (int i=0;i<n_pt_bins;i++) pt_zero[i]= 0.;

  TGraphAsymmErrors* systematic_errors = new TGraphAsymmErrors(n_pt_bins, x_values, y_values, x_error_low, x_error_high, y_errors, y_errors);
  TCanvas c;
  systematic_errors->SetMarkerColor(4);
  systematic_errors->SetMarkerStyle(5);
  systematic_errors->GetXaxis()->SetTitle("pT(GeV)");
  systematic_errors->Draw("AP");
  systematic_errors->SetTitle("Efficiency systematic error");

  if (particle == 0){
     c.SaveAs("~/work2/BinQGP/results/Bu/efficiency/plots/systematic_error.gif");
     c.SaveAs("~/work2/BinQGP/results/Bu/efficiency/plots/systematic_error.pdf");
  }
  else if (particle == 1){
     c.SaveAs("~/work2/BinQGP/results/Bs/efficiency/plots/systematic_error.gif");
     c.SaveAs("~/work2/BinQGP/results/Bs/efficiency/plots/systematic_error.pdf");
  }

  TFile* f1 = particle ? new TFile("~/work2/BinQGP/results/Bs/efficiency/root_files/efficiency_systematic_errors.root", "recreate") : new TFile("~/work2/BinQGP/results/Bu/efficiency/root_files/efficiency_systematic_errors.root", "recreate");
  f1->cd();
  systematic_errors->Write();
  f1->Write();
  f1->ls();
  f1->Close();
  
  return 0;  
}
