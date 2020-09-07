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

void plot_xsection(int bin_n,double* pt_m, double* pt_l, double* pt_h, double* x_sec, double* stat, double* syst);
//void error_syst_final(double s_errors_bin);

#define particle 1

void x_section(){

  //Raw yield files
  TFile* f_raw_yield = particle ?  new TFile("~/work2/BinQGP/results/Bs/Bpt/pT.root") : new TFile("~/work2/BinQGP/results/Bu/Bpt/pT.root");
  
  //Efficiency files
  TFile* f_efficiency = particle ?  new TFile("~/work2/BinQGP/results/Bs/efficiency/root_files/efficiency0.root") : new TFile("~/work2/BinQGP/results/Bu/efficiency/root_files/efficiency0.root");

  //Efficiency systematic error files
  TFile* f_eff_syst = particle ? new TFile("~/work2/BinQGP/results/Bs/efficiency/root_files/efficiency_systematic_errors.root") : new TFile("~/work2/BinQGP/results/Bu/efficiency/root_files/efficiency_systematic_errors.root");
	
  const double branching_fraction = particle ? 0.0000313 : 0.0000599;//check
  const double branching_fraction_error = particle ?  0.0000030 : 0.0000023; //absolute error

  const double luminosity = 302.3;  //pb^-1
  const double luminosity_error = 0.023;  //(relative error)

  //yields
  TGraphAsymmErrors* raw_yield = (TGraphAsymmErrors*)f_raw_yield->Get("Graph;1"); //has absolute stat errors (y error bars)
  TGraphAsymmErrors* raw_yield_s = (TGraphAsymmErrors*)f_raw_yield->Get("Graph;2"); //has absolute syst errors (y error bars)

  double* pt_bins = raw_yield->GetX();
  const int n_pt_bins = raw_yield->GetN();

  //efficiencies
  TEfficiency* efficiency = particle ? new TEfficiency("efficiency_Bs", "efficiency_Bs", n_pt_bins, pt_bins) : new TEfficiency("efficiency_Bu", "efficiency_Bu", n_pt_bins, pt_bins);
  efficiency = (TEfficiency*)f_efficiency->Get("hist_tot_noweights_clone");//has efficiencies stat errors (y error bars)

  //systematic error
  TGraphErrors* efficiency_syst = (TGraphErrors*)f_eff_syst->Get("Graph");// has efficiencies syst errors (y values)

  ///// VARIABLES /////

  double* raw = raw_yield->GetY(); // yield value
  double* raw_errX_low = raw_yield->GetEXlow(); //lower pt  error 
  double* raw_errX_high = raw_yield->GetEXhigh(); //upper pt  error 

  double* eff_s = efficiency_syst->GetY(); //efficiency syst error (relative)

  double n;
  double eff;
  double x_sec0;
  double norm;
  double raw_stat;
  double raw_syst;
  double eff_syst;
  double lumi_syst = luminosity_error; // (relative) 
  double branch_syst = branching_fraction_error/branching_fraction; // (relative)
  double tot_syst;

  double x_sec[n_pt_bins]; //cross section central value 
  double x_sec_stat[n_pt_bins]; //cross section stat error
  double x_sec_syst[n_pt_bins]; //cross section syst error
 

  for(int i = 0; i < n_pt_bins; i++)
    {
      /// Central Value ///
      n = raw[i]; //yield value
      eff = efficiency->GetEfficiency(i+1); //efficiency value
      norm = eff*branching_fraction*luminosity;
      x_sec0 = n/norm; //cross section (central value)
      x_sec[i] = x_sec0; 

      /// Stat Error ///
      raw_stat = raw_yield->GetErrorY(i);
      x_sec_stat[i] = raw_stat/norm;
      
      /// Syst Error /// (relative)
      raw_syst = (raw_yield_s->GetErrorY(i))/n;
      eff_syst = fabs(eff_s[i]);
  
      tot_syst = raw_syst*raw_syst + eff_syst*eff_syst + lumi_syst*lumi_syst + branch_syst*branch_syst;
      x_sec_syst[i] = sqrt(tot_syst)*x_sec0;
    }

  /*
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
  */
 
  plot_xsection(n_pt_bins, pt_bins, raw_errX_low, raw_errX_high, x_sec, x_sec_stat, x_sec_syst); 
}
 

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

  TFile* f = particle ? new TFile("~/work2/BinQGP/results/Bs/x_section/x_section.root", "recreate") : new TFile("~/work2/BinQGP/results/Bu/x_section/x_section.root", "recreate");
  f->cd();
 
  mg->Add(g_stat);
  mg->Add(g_syst);
  mg->Draw("AP");
  mg->GetXaxis()->SetTitle("p_{T}(B) [GeV]");
  mg->GetYaxis()->SetTitle("#frac{d#sigma}{dp_{T}} [pb/GeV]");
 
  f->Write();

  if (particle == 0){
     c.SaveAs("~/work2/BinQGP/results/Bu/x_section/x_section.gif");
     c.SaveAs("~/work2/BinQGP/results/Bu/x_section/x_section.pdf");
  }
  else if (particle == 1){
     c.SaveAs("~/work2/BinQGP/results/Bs/x_section/x_section.gif");
     c.SaveAs("~/work2/BinQGP/results/Bs/x_section/x_section.pdf");
  }
  f->Close();
}

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

