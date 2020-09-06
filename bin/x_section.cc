/////////////////////////////////////////////////////////////////////////
//
// X-section, ratio and systematics
//
//
// August 2019
//
/////////////////////////////////////////////////////////////////////////


#include <sstream>
#include <vector>
#include <TStyle.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TNtupleD.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <TSystem.h>
#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooConstVar.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooChebychev.h>
#include <RooBernstein.h>
#include <RooExponential.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooBifurGauss.h>
#include "TMath.h"
#include <RooGenericPdf.h>
#include "TRatioPlot.h"
#include <RooBifurGauss.h>
#include <RooProduct.h>
#include <RooHist.h>
#include "RooStats/SPlot.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooMCStudy.h"
#include <iostream>
#include <TF1.h>
#include <RooPolynomial.h>
#include <fstream>
#include <TGraphAsymmErrors.h>
#include <TEfficiency.h>
#include <TGraphErrors.h>

using namespace RooStats;
using namespace RooFit;
using namespace std;

//functions

//particle
// 0 = Bu
// 1 = Bs

//void plot_xsection(int bin_n,double* pt_m, double* pt_l, double* pt_h, double* x_sec, double* stat, double* syst);
//void error_syst_final(double s_errors_bin);

#define particle 1

int main(){

  //Raw yield files
  TFile* f_raw_yield = particle ?  new TFile("~/work2/BinQGP/results/Bs/Bpt/pT.root") : new TFile("~/work2/BinQGP/results/Bu/Bpt/pT.root");
  
  //Efficiency files
  TFile* f_efficiency = particle ?  new TFile("~/work2/BinQGP/results/Bs/efficiency/root_files/efficiency0.root") : new TFile("~/work2/BinQGP/results/Bu/efficiency/root_files/efficiency0.root");

  //Efficiency systematic error files
  TFile* f_eff_syst = particle ? new TFile("~/work2/BinQGP/results/Bs/efficiency/root_files/efficiency_systematic_errors.root") : new TFile("~/work2/BinQGP/results/Bu/efficiency/root_files/efficiency_systematic_errors.root");

  const double branching_fraction = particle ? 0.0000313 : 0.0000599;

  const double branching_fraction_error = particle ?  0.0000030 : 0.0000023;

  const double luminosity = 0.0000000000003023;  //wrong value
  const double luminosity_error = 0.0000000001;  //wrong value

  double pt_bins[] = {5, 10, 15, 20, 50};
  double n_pt_bins = 4;

  TGraphAsymmErrors* raw_yield = (TGraphAsymmErrors*)f_raw_yield->Get("Graph");

  TEfficiency* efficiency = particle ? new TEfficiency("efficiency_Bs", "efficiency_Bs", n_pt_bins, pt_bins) : new TEfficiency("efficiency_Bu", "efficiency_Bu", n_pt_bins, pt_bins);

  efficiency = (TEfficiency*)f_efficiency->Get("hist_tot_noweights_clone");
  TGraphErrors* eff_syst = (TGraphErrors*)f_eff_syst->Get("Graph");

  TH1F* x_section = particle ? new TH1F("x_section_Bs", "x_section_Bs", n_pt_bins, pt_bins) : new TH1F("x_section_Bu", "x_section_Bu", n_pt_bins, pt_bins);

  double x_sec[4];
  double x_sec0;

  double n;
  double eff;

  double* raw = raw_yield->GetY();

  double* eff_s = eff_syst->GetY();

 /*
  for(int i = 0; i < n_pt_bins; i++)
    {
      cout << "Bin " << i+1 << endl;
      if (particle == 1){
         cout << "Efficiency Bs = " << efficiency->GetBinContent(i+1) << endl;
         cout << "Yield Bs = " << raw[i] << endl;
      }
      else if (particle == 0){
         cout << "Efficiency Bu = " << efficiency->GetBinContent(i+1) << endl;
         cout << "Yield Bu = " << raw[i] << endl;
      }
      cout << endl;
    }
 */ 
  
  for(int i = 0; i < n_pt_bins; i++)
    {
      n = raw[i];
      cout << "n ="<< n << endl;
      eff = efficiency->GetEfficiency(i+1);
      cout << "eff =" << eff << endl;
      x_sec0 = n/(eff*branching_fraction*luminosity);
      cout << "x_sec0 =" << x_sec0 << endl;
      x_sec[i] = x_sec0;
      cout << "x_sec[i] =" << x_sec[i] << endl;
      x_section->SetBinContent(i, x_sec0);
      cout << "bin content =" << x_section->GetBinContent(i) << endl;
    }
  
 
  double syst_errors[5][4];

  //Eff-Acc systematic
  syst_errors[0][0] = eff_s[0];
  syst_errors[0][1] = eff_s[1];
  syst_errors[0][2] = eff_s[2];
  syst_errors[0][3] = eff_s[3];

  //Fit systematic ---------- Luminosity systematic --------- Branching fraction systematic
  for(int i = 0; i < 4; i++)
    {
      syst_errors[1][i] = 0.005;
      syst_errors[2][i] = luminosity_error;
      syst_errors[3][i] = 0.0000023;
    }
  
  TCanvas c;
  //x_section->SetMinimum(100000000000000000);
  //x_section->SetMaximum(2000000000000000000);
  x_section->Draw();
  if (particle == 0){
     c.SaveAs("~/work2/BinQGP/results/Bu/x_section/x_section.gif");
     c.SaveAs("~/work2/BinQGP/results/Bu/x_section/x_section.pdf");
  }
  else if (particle == 1){
     c.SaveAs("~/work2/BinQGP/results/Bs/x_section/x_section.gif");
     c.SaveAs("~/work2/BinQGP/results/Bs/x_section/x_section.pdf");
  }

  TFile* f = particle ? new TFile("~/work2/BinQGP/results/Bs/x_section/x_section.root", "recreate") : new TFile("~/work2/BinQGP/results/Bu/x_section/x_section.root", "recreate");
  f->cd();
  x_section->Write();
  f->Write();

  

  return 0;
  
}
 
/*
void plot_xsection(int bin_n,double* pt_m, double* pt_l, double* pt_h, double* x_sec, double* stat, double* syst){
  TCanvas c;
  TMultiGraph* mg = new TMultiGraph();
  
  TGraphAsymmErrors* g_stat = new TGraphAsymmErrors(bin_n,pt_m,x_sec,pt_l,pt_h,stat,stat);
  g_stat->SetTitle("");
  g_stat->SetMarkerColor(4);
  g_stat->SetMarkerStyle(1);
  g_stat->SetLineColor(1);
  
  double pt_zero[bin_n];
  for (int i=0;i<bin_n;i++) pt_zero[i]= 0.;
  
  TGraphAsymmErrors* g_syst= new TGraphAsymmErrors(bin_n,pt_m,x_sec,pt_zero,pt_zero,syst,syst);
  g_syst->SetTitle("");
  g_syst->SetMarkerColor(4);
  g_syst->SetMarkerStyle(1);
  g_syst->SetLineColor(2);  
  
  mg->Add(gr);
  mg->Add(grs);
  mg->Draw("AP");
  mg->GetXaxis()->SetTitle("p_{T}(B) [GeV]");
  mg->GetYaxis()->SetTitle("X-section [GeV^{-1}]");
}
*/
/*
//function that evaluates the final systematic error
//s_errors_bin = array com os valores dos 5 erros por bin
void error_syst_final(double s_errors_bin){ 
  const int n_s_errors = 5;
  //N, B, L, eff+A, eff MC
  //number of sources of systematics errors
  double s_errors_bin [n_pt_bins][n_s_errors];
  //value of the systematic error for each source per pT bin -->need update
  double final_syst[n_pt_bins];
  //value of the systematic error per bin
  double sum_pow_syst[n_pt_bins];
  //sum of the squares ofthe syst errors

   for(int s = 0; s<n_pt_bins; s++){
    sum_pow_syst[s]=0;
  } 
 
  //loop through the number of bins
  for(int i=0;i<n_pt_bins;i++){
    //loop through the array of sources of error
    for(int k=0;k<n_s_errors;k++){
      sum_pow_syst[i] += pow(s_errors_bin[i][k],2); 
    }
    final_syst[i] = sqrt(sum_pow_syst[i]);
  }
}
*/
//error_syst_final ends

