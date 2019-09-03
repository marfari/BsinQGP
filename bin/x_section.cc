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
//  0 = Bu
// 1 = Bs

int main(){

  TFile* f_raw_yield_Bs = new TFile("./results/Bs/Bpt/pT.root");
  TFile* f_raw_yield_Bu = new TFile("/lstore/cms/ev19u032/pT_Bu_meson.root");
  
  TFile* f_efficiency_Bs = new TFile("/home/t3cms/julia/LSTORE/CMSSW_7_5_8_patch5/src/UserCode/Bs_analysis/for_students/MCstudiesPbPbPtBin.root");
  TFile* f_efficiency_Bu = new TFile("/home/t3cms/julia/LSTORE/CMSSW_7_5_8_patch5/src/UserCode/Bs_analysis/for_students/MCstudiesPbPb_Bsbin.root");

  const double branching_fraction_Bs = 0.0000313;
  const double branching_fraction_Bu = 0.0000599;

  const double branching_fraction_error_Bs = 0.0000030;
  const double branching_fraction_error_Bu = 0.0000023;

  const double luminosity = 0.0000000015;  //wrong value
  const double luminosity_error = 0.0000000001;  //wrong value

  double pt_bins[] = {5, 10, 15, 20, 50};
  double n_pt_bins = 4;

  TGraphAsymmErrors* raw_yield_Bs = (TGraphAsymmErrors*)f_raw_yield_Bs->Get("Graph");
  TGraphAsymmErrors* raw_yield_Bu = (TGraphAsymmErrors*)f_raw_yield_Bu->Get("Graph");

  TEfficiency* efficiency_Bs = new TEfficiency("efficiency_Bs", "efficiency_Bs", n_pt_bins, pt_bins);
  TEfficiency* efficiency_Bu = new TEfficiency("efficiency_Bu", "efficiency_Bu", n_pt_bins, pt_bins);
  efficiency_Bs = (TEfficiency*)f_efficiency_Bs->Get("hEff");
  efficiency_Bu = (TEfficiency*)f_efficiency_Bu->Get("hEff");

  TH1F* x_section_Bu = new TH1F("x_section_Bu", "x_section_Bu", n_pt_bins, pt_bins);
  TH1F* x_section_Bs = new TH1F("x_section_Bs", "x_section_Bs", n_pt_bins, pt_bins);

  double x_sec;

  double n;
  double eff;

  for(int i = 0; i < 4; i++)
    {
      n = raw_yield_Bs->
    }

  return 0;
}
  
