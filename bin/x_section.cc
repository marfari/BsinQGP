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
using namespace RooStats;
using namespace RooFit;
using namespace std;

//functions

//particle
// 0 = Bu
// 1 = Bs

//void plot_xsection(int bin_n,double* pt_m, double* pt_l, double* pt_h, double* x_sec, double* stat, double* syst);
//void error_syst_final(double s_errors_bin);


int main(){

  TFile* f_raw_yield_Bs = new TFile("./results/Bs/Bpt/pT.root");
  TFile* f_raw_yield_Bu = new TFile("/lstore/cms/ev19u032/pT_Bu_meson.root");
  
  TFile* f_efficiency_Bs = new TFile("/home/t3cms/julia/LSTORE/CMSSW_7_5_8_patch5/src/UserCode/Bs_analysis/for_students/MCstudiesPbPbPtBin.root");
  TFile* f_efficiency_Bu = new TFile("/home/t3cms/julia/LSTORE/CMSSW_7_5_8_patch5/src/UserCode/Bs_analysis/for_students/MCstudiesPbPb_Bsbin.root");

  const double branching_fraction_Bs = 0.0000313;
  const double branching_fraction_Bu = 0.0000599;

  //const double branching_fraction_error_Bs = 0.0000030;
  //const double branching_fraction_error_Bu = 0.0000023;

  const double luminosity = 0.0000000015;  //wrong value
  //const double luminosity_error = 0.0000000001;  //wrong value

  double pt_bins[] = {5, 10, 15, 20, 50};
  double n_pt_bins = 4;

  TGraphAsymmErrors* raw_yield_Bs = (TGraphAsymmErrors*)f_raw_yield_Bs->Get("Graph");
  TGraphAsymmErrors* raw_yield_Bu = (TGraphAsymmErrors*)f_raw_yield_Bu->Get("Graph");

  TH1D* efficiency_Bs = new TH1D("efficiency_Bs", "efficiency_Bs", n_pt_bins, pt_bins);
  efficiency_Bs = (TH1D*)f_efficiency_Bs->Get("hEff");
  TH1D* efficiency_Bu = new TH1D("efficiency_Bu", "efficiency_Bu", n_pt_bins, pt_bins);
  efficiency_Bu = (TH1D*)f_efficiency_Bu->Get("hEff");

  //TH1F* x_section_Bu = new TH1F("x_section_Bu", "x_section_Bu", n_pt_bins, pt_bins);
  //TH1F* x_section_Bs = new TH1F("x_section_Bs", "x_section_Bs", n_pt_bins, pt_bins);

  double x_sec_Bu[4];
  double x_sec_Bs[4];
  double x_sec0;

  double n;
  double eff;

  double* raw_Bs_y = raw_yield_Bs->GetY();
  double* raw_Bu_y = raw_yield_Bu->GetY();

  /*
  for(int i = 0; i < n_pt_bins; i++)
    {
      cout << "Bin " << i+1 << endl;
      cout << "Efficiency Bs = " << efficiency_Bs->GetBinContent(i+1) << endl;
      cout << "Efficiency Bu = " << efficiency_Bu->GetBinContent(i+1) << endl;
      cout << "Yield Bs = " << raw_Bs_y[i] << endl;
      cout << "Yield Bu = " << raw_Bu_y[i] << endl;
      cout << endl;
    }
  */

  for(int i = 0; i < n_pt_bins; i++)
    {
      //Bu
      n = raw_Bu_y[i];
      eff = efficiency_Bu->GetBinContent(i+1);
      x_sec0 = n/(eff*branching_fraction_Bu*luminosity);
      x_sec_Bu[i] = x_sec0;
      cout << x_sec_Bu[i] << endl;
      cout << endl;

      //Bs
      n = raw_Bs_y[i];
      eff = efficiency_Bs->GetBinContent(i+1);
      x_sec0 = n/(eff*branching_fraction_Bs*luminosity);
      x_sec_Bs[i] = x_sec0;
      cout << x_sec_Bs[i] << endl;
      cout << endl;
    }

  /*
  TCanvas Bu_c;
  x_section_Bu->SetMinimum(100000000000000);
  x_section_Bu->SetMaximum(200000000000000000);
  x_section_Bu->Draw();
  Bu_c.SaveAs("./results/Bu/x_section/x_section.gif");
  Bu_c.SaveAs("./results/Bu/x_section/x_section.pdf");
  TFile* Bu_f = new TFile("./results/Bu/x_section/x_section.root", "recreate");
  Bu_f->cd();
  x_section_Bu->Write();
  Bu_f->Write();

  TCanvas Bs_c;
  x_section_Bs->Draw();
  Bs_c.SaveAs("./results/Bs/x_section/x_section.gif");
  Bs_c.SaveAs("./results/Bs/x_section/x_section.pdf");
  TFile* Bs_f = new TFile("./results/Bs/x_section/x_section.root", "recreate");
  Bs_f->cd();
  x_section_Bs->Write();
  Bs_f->Write();
  */

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

//error_syst_final ends
*/
