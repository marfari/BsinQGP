/////////////////////////////////////////////////////////////////////////
//
// Sideband subtraction and SPlot methods
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
using namespace RooStats;
using namespace RooFit;
using namespace std;

#define particle 0//0 = B+;    1 = Bs;

std::vector<TH1D*> sideband_subtraction(RooWorkspace* w, int* n, int n_var);
std::vector<TH1D*> splot_method(RooWorkspace& w, int* n, TString* label, int n_var);

void set_up_workspace_variables(RooWorkspace& w);
TH1D* create_histogram_mc(RooRealVar var, TTree* t, int n); //mc
TH1D* create_histogram(RooRealVar var,TString name, double factor, RooDataSet* reduced, RooDataSet* central, RooDataSet* total, int n); //data
void read_data(RooWorkspace& w, TString f_input);
void build_pdf (RooWorkspace& w);
void plot_complete_fit(RooWorkspace& w);
void do_splot(RooWorkspace& w);
TH1D* make_splot(RooWorkspace& w, int n, TString label);
void validate_fit(RooWorkspace* w);

int main(){
  
  int n_var;

  TString input_file_data = particle ? "/home/t3cms/julia/LSTORE/CMSSW_7_5_8_patch5/src/UserCode/Bs_analysis/prefiltered_trees/selected_data_ntphi_PbPb_2018_corrected_test_new.root" : "/home/t3cms/julia/LSTORE/CMSSW_7_5_8_patch5/src/UserCode/Bs_analysis/prefiltered_trees/selected_data_ntKp_PbPb_2018_corrected_test.root";
  TString input_file_mc = particle ? "/home/t3cms/julia/LSTORE/CMSSW_7_5_8_patch5/src/UserCode/Bs_analysis/prefiltered_trees/selected_mc_ntphi_PbPb_2018_corrected_test_new.root" : "/home/t3cms/julia/LSTORE/CMSSW_7_5_8_patch5/src/UserCode/Bs_analysis/prefiltered_trees/selected_mc_ntKp_PbPb_2018_corrected_test.root";

  std::vector<TH1D*> histos_data;
  std::vector<TH1D*> histos_mc;
  std::vector<TH1D*> histos_splot;

#if particle == 0
  int n_bins[] = {20, 20, 10, 10, 10, 10, 10, 10, 10, 10, 35, 10, 70, 10, 10, 10, 10, 10,10,10,10};
  TString variables[] = {"Bpt","By","Btrk1eta","Btrk1Y","Btrk1pt","Bmu1eta","Bmu2eta","Bmu1pt","Bmu2pt","Bchi2cl", "BsvpvDistance", "BsvpvDistance_Err","Balpha","Btrk1Dz1","BvtxX", "BvtxY", "Btrk1DzError1", "Btrk1Dxy1", "Btrk1DxyError1", "Bd0","Bd0err"};
#elif particle == 1
  int n_bins[] = {10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10};
  TString variables[] = {"Bpt","By","Btrk1eta","Btrk1Y", "Btrk2eta", "Btrk1pt", "Btrk2pt", "Bmu1eta","Bmu2eta","Bmu1pt","Bmu2pt","Bchi2cl", "Bmumumass", "Btrktrkmass", "BsvpvDistance", "BsvpvDistance_Err","Balpha","BvtxX", "BvtxY"};
#endif

  int n_n_bins = sizeof(n_bins)/sizeof(n_bins[0]);
  //std::cout << n_n_bins << std::endl;
  int n_variables = sizeof(variables)/sizeof(variables[0]);
  //std::cout << n_variables << std::endl;

  if(n_n_bins != n_variables){
    std::cout << "Error: number of bins does not correspond to number of variables." << std::endl;
    return 0;
  }

  n_var = n_variables;
  
  RooWorkspace* ws = new RooWorkspace("ws");
  set_up_workspace_variables(*ws);
  read_data(*ws,input_file_data);
  build_pdf(*ws);
  plot_complete_fit(*ws);

  //validate_fit(ws);

  //sideband_sub histograms
  histos_data = sideband_subtraction(ws, n_bins, n_var);
  
  TFile *fin_mc = new TFile(input_file_mc);
  TTree* t1_mc = particle ? (TTree*)fin_mc->Get("ntphi") : (TTree*)fin_mc->Get("ntKp");

  std::vector<TString> names;


  for(int i=0; i<(int)histos_data.size(); i++){
    std::cout<< "Var names: "<< histos_data[i]->GetName()<<std::endl;
  }

  ///////////HERE///////////
  for(int i=0; i<(int)histos_data.size(); i++){
    histos_mc.push_back(create_histogram_mc((*ws->var(histos_data[i]->GetName())), t1_mc, n_bins[i]));
    names.push_back(TString(histos_data[i]->GetName()));
  }
  

  //splot histograms
  do_splot(*ws);
  histos_splot = splot_method(*ws,n_bins,variables, n_var);
  
  //sideband subtraction method vs. monte carlo
  for(int i=0; i<(int)histos_data.size(); i++)
    {
      TCanvas c;
      histos_mc[i]->SetXTitle(TString(histos_data[i]->GetName()));
      histos_mc[i]->SetStats(0);
      histos_data[i]->SetStats(0);

      //normalization
      histos_mc[i]->Scale(1/histos_mc[i]->Integral());
      histos_data[i]->Scale(1/histos_data[i]->Integral());

      histos_mc[i]->GetYaxis()->SetRangeUser(2*histos_data[i]->GetMinimum(),2*histos_mc[i]->GetMaximum());
      histos_mc[i]->Draw();
      histos_data[i]->Draw("same");
      auto rp = new TRatioPlot(histos_data[i], histos_mc[i], "divsym");
      c.SetTicks(0, 1);
      rp->SetH1DrawOpt("E");
      rp->Draw("nogrid");
      rp->GetLowerRefYaxis()->SetTitle("Data(ss)/MC");
      rp->GetUpperRefYaxis()->SetTitle("normalized entries");
      c.Update();
	
      TLegend* leg;

      leg = new TLegend(0.7, 0.7, 0.9, 0.9);
      leg->AddEntry(histos_data[i]->GetName(), "S. Subtraction", "l");
      leg->AddEntry(histos_mc[i]->GetName(), "Monte Carlo", "l");
      leg->SetTextSize(0.03);
      leg->Draw("same");

      if(particle == 0){
      c.SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/teste2/mc_validation_plots/ss_mc/" + names[i]+"_mc_validationB+.pdf");
      }else if(particle == 1){
      c.SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/Bs/mc_validation_plots/ss_mc/" + names[i]+"_mc_validationBs.pdf");
      }

      leg->Delete();

    }

  //sideband subtraction method vs. monte carlo vs splot
  for(int i=0; i<(int)histos_data.size(); i++)
    {
      TCanvas a;
      histos_mc[i]->SetXTitle(TString(histos_data[i]->GetName()));
      histos_splot[i]->SetXTitle(TString(histos_data[i]->GetName()));
      histos_mc[i]->SetStats(0);
      histos_splot[i]->SetStats(0);
      histos_data[i]->SetStats(0);

      //normalization
      histos_mc[i]->Scale(1/histos_mc[i]->Integral());
      histos_splot[i]->Scale(1/histos_splot[i]->Integral());
      histos_data[i]->Scale(1/histos_data[i]->Integral());

      histos_mc[i]->GetYaxis()->SetRangeUser(2*histos_data[i]->GetMinimum(),2*histos_mc[i]->GetMaximum());
      histos_mc[i]->Draw();
      histos_splot[i]->Draw("same");
      histos_data[i]->Draw("same");

      //TRatio Plot Data(sideband_sub)/MC
      auto rp = new TRatioPlot(histos_data[i], histos_mc[i], "divsym");
      a.SetTicks(0, 1);
      rp->SetH1DrawOpt("E");
      rp->Draw("nogrid");
      rp->GetLowerRefYaxis()->SetTitle("Data(ss)/MC");
      rp->GetUpperRefYaxis()->SetTitle("normalized entries");
      a.Update();
	
      TLegend* leg;

      leg = new TLegend(0.7, 0.7, 0.9, 0.9);
      leg->AddEntry(histos_data[i]->GetName(), "S. Subtraction", "l");
      leg->AddEntry(histos_mc[i]->GetName(), "Monte Carlo", "l");
      leg->AddEntry(histos_splot[i]->GetName(), "SPlot", "l");
      leg->SetTextSize(0.03);
      leg->Draw("same");

      if(particle == 0){
	a.SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/teste2/mc_validation_plots/ss_mc_sp/"+names[i]+"_mc_validationB+.pdf");
      }else if(particle == 1){
	a.SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/Bs/mc_validation_plots/ss_mc_sp/"+names[i]+"_mc_validationBs.pdf");
      }
      leg->Delete();

    }
 
  //Splot vs MC comparison

  for(int i=0; i<(int)histos_data.size(); i++)
    {
      TCanvas a;
      histos_mc[i]->SetXTitle(TString(histos_data[i]->GetName()));
      histos_mc[i]->SetYTitle("normalized entries");
      histos_splot[i]->SetXTitle(TString(histos_data[i]->GetName()));
      histos_mc[i]->SetStats(0);
      histos_splot[i]->SetStats(0);
      histos_mc[i]->Scale(1/histos_mc[i]->Integral());
      histos_splot[i]->Scale(1/histos_splot[i]->Integral());
      histos_mc[i]->GetYaxis()->SetRangeUser(2*histos_data[i]->GetMinimum(),2*histos_mc[i]->GetMaximum());
      histos_mc[i]->Draw();
      histos_splot[i]->Draw("same");

      //--TRATIOS--//
      
      //TRatio Plot Data(sideband_sub)/MC
      auto rp = new TRatioPlot(histos_splot[i], histos_mc[i], "divsym");
      a.SetTicks(0, 1);
      rp->SetH1DrawOpt("E");
      rp->Draw("nogrid");
      rp->GetLowerRefYaxis()->SetTitle("Data(sp)/MC");
      rp->GetUpperRefYaxis()->SetTitle("normalized entries");
      a.Update();
     
      TLegend* leg;

      leg = new TLegend(0.7, 0.7, 0.9, 0.9);
      leg->AddEntry(histos_mc[i]->GetName(), "Monte Carlo", "l");
      leg->AddEntry(histos_splot[i]->GetName(), "SPlot", "l");
      leg->SetTextSize(0.03);
      leg->Draw("same");

      if(particle == 0){
	a.SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/teste2/mc_validation_plots/mc_sp/" + names[i]+"_mc_validationB+.pdf");
      }else if(particle == 1){
	a.SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/Bs/mc_validation_plots/mc_sp/"+names[i]+"_mc_validationBs.pdf");
      }

      leg->Delete();
      histos_mc[i]->Delete();
      histos_data[i]->Delete();
      histos_splot[i]->Delete();
    }     
  
}

void read_data(RooWorkspace& w, TString f_input){

  TFile* fin_data = new TFile(f_input);
  //TNtupleD* _nt = (TNtupleD*)fin_data->Get("ntKp");
  TTree* t1_data = particle ? (TTree*)fin_data->Get("ntphi") : (TTree*)fin_data->Get("ntKp"); //ntKp
  RooArgList arg_list ("arg_list");

  arg_list.add(*(w.var("Bmass")));
  arg_list.add(*(w.var("Bpt")));
  arg_list.add(*(w.var("By")));
  arg_list.add(*(w.var("Btrk1eta")));
  arg_list.add(*(w.var("Btrk1Y")));
  if(particle == 1){arg_list.add(*(w.var("Btrk2eta")));}
  arg_list.add(*(w.var("Btrk1pt"))); 
  if(particle == 1){arg_list.add(*(w.var("Btrk2pt")));}
  arg_list.add(*(w.var("Bmu1eta")));
  arg_list.add(*(w.var("Bmu2eta")));
  arg_list.add(*(w.var("Bmu1pt")));
  arg_list.add(*(w.var("Bmu2pt")));
  arg_list.add(*(w.var("Bchi2cl")));
  if(particle == 1){arg_list.add(*(w.var("Bmumumass")));}
  if(particle == 1){arg_list.add(*(w.var("Btrktrkmass")));}
  arg_list.add(*(w.var("BsvpvDistance")));
  arg_list.add(*(w.var("BsvpvDistance_Err")));
  arg_list.add(*(w.var("Balpha")));
  if(particle == 0){arg_list.add(*(w.var("Btrk1Dz1")));}
  arg_list.add(*(w.var("BvtxX")));
  arg_list.add(*(w.var("BvtxY")));
  if(particle == 0){
  arg_list.add(*(w.var("Btrk1DzError1")));
  arg_list.add(*(w.var("Btrk1Dxy1")));
  arg_list.add(*(w.var("Btrk1DxyError1")));
  arg_list.add(*(w.var("Bd0")));
  arg_list.add(*(w.var("Bd0err")));
  }

  RooDataSet* data = new RooDataSet("data","data",t1_data,arg_list);

  w.import(*data, Rename("data"));
}

void build_pdf(RooWorkspace& w) {

  RooRealVar Bmass = *(w.var("Bmass"));
  RooDataSet* data = (RooDataSet*) w.data("data");

  RooDataSet* reduceddata_central;

  double left = 5.15;
  double right = 5.4;
  double mass_peak = 5.265;

  reduceddata_central = (RooDataSet*) data->reduce(Form("Bmass>%lf",left));
  reduceddata_central = (RooDataSet*) reduceddata_central->reduce(Form("Bmass<%lf",right));

  //SIGNAL//

  RooRealVar mean("mean","mean",mass_peak,5.26,5.29);
  RooRealVar sigma1("sigma1","sigma1",0.021,0.020,0.030);
  RooGaussian signal1("signal1","signal_gauss1",Bmass,mean,sigma1);
  RooRealVar sigma2("sigma2","sigma2",0.011,0.010,0.020);
  RooGaussian signal2("signal2","signal_gauss2",Bmass,mean,sigma2);
  RooRealVar cofs("cofs", "cofs", 0.5, 0., 1.);
  RooAddPdf signal("signal", "signal", RooArgList(signal1,signal2),cofs);

  //BACKGROUND//

  //error function
  RooRealVar m_nonprompt_scale("m_nonprompt_scale", "m_nonprompt_scale", 4.74168e-02, 0, 1);
  RooRealVar m_nonprompt_shift("m_nonprompt_shift", "m_nonprompt_shift", 5.14425, 4.5, 6.);
  
  m_nonprompt_shift.setConstant(kTRUE);
  m_nonprompt_scale.setConstant(kTRUE);

  RooGenericPdf erf("erf","erf","TMath::Erfc((Bmass-m_nonprompt_shift)/m_nonprompt_scale)", RooArgList(Bmass,m_nonprompt_scale,m_nonprompt_shift));
 
  //exponential
  RooRealVar lambda("lambda","lambda",-2.,-5.,0.0);
  RooExponential fit_side("fit_side", "fit_side_exp", Bmass, lambda);

  //jpsi_pi component
  RooRealVar m_jpsipi_mean1("m_jpsipi_mean1","m_jpsipi_mean1",5.34693e+00,Bmass.getAsymErrorLo(),Bmass.getAsymErrorHi());
  RooRealVar m_jpsipi_mean2("m_jpsipi_mean2","m_jpsipi_mean2",5.46876e+00,Bmass.getAsymErrorLo(),Bmass.getAsymErrorHi());
  RooRealVar m_jpsipi_mean3("m_jpsipi_mean3","m_jpsipi_mean3",5.48073e+00,Bmass.getAsymErrorLo(),Bmass.getAsymErrorHi());
  RooRealVar m_jpsipi_sigma1l("m_jpsipi_sigma1l","m_jpsipi_sigma1l",2.90762e-02,0.010,0.150);
  RooRealVar m_jpsipi_sigma1r("m_jpsipi_sigma1r","m_jpsipi_sigma1r",6.52519e-02,0.010,0.150);
  RooRealVar m_jpsipi_sigma2("m_jpsipi_sigma2","m_jpsipi_sigma2",9.94712e-02,0.020,0.500);

  RooRealVar m_jpsipi_sigma3("m_jpsipi_sigma3","m_jpsipi_sigma3",3.30152e-01,0.020,0.500);
  RooRealVar m_jpsipi_fraction2("m_jpsipi_fraction2","m_jpsipi_fraction2",2.34646e-01,0.0,1.0);
  RooRealVar m_jpsipi_fraction3("m_jpsipi_fraction3","m_jpsipi_fraction3",1.14338e-01,0.0,1.0);

  m_jpsipi_mean1.setConstant(kTRUE);
  m_jpsipi_mean2.setConstant(kTRUE);
  m_jpsipi_mean3.setConstant(kTRUE);
  m_jpsipi_sigma1l.setConstant(kTRUE);
  m_jpsipi_sigma1r.setConstant(kTRUE);
  m_jpsipi_sigma2.setConstant(kTRUE);
  m_jpsipi_sigma3.setConstant(kTRUE);
  m_jpsipi_fraction2.setConstant(kTRUE);
  m_jpsipi_fraction3.setConstant(kTRUE);

  RooBifurGauss m_jpsipi_gaussian1("m_jpsipi_gaussian1","m_jpsipi_gaussian1",Bmass,m_jpsipi_mean1,m_jpsipi_sigma1l,m_jpsipi_sigma1r);
  RooGaussian m_jpsipi_gaussian2("m_jpsipi_gaussian2","m_jpsipi_gaussian2",Bmass,m_jpsipi_mean2,m_jpsipi_sigma2);
  RooGaussian m_jpsipi_gaussian3("m_jpsipi_gaussian3","m_jpsipi_gaussian3",Bmass,m_jpsipi_mean3,m_jpsipi_sigma3);

  RooAddPdf jpsipi("jpsipi","jpsipi",RooArgList(m_jpsipi_gaussian3,m_jpsipi_gaussian2,m_jpsipi_gaussian1),RooArgList(m_jpsipi_fraction3,m_jpsipi_fraction2));

  Bmass.setRange("all", Bmass.getMin(),Bmass.getMax());
  Bmass.setRange("right",right,Bmass.getMax());
  Bmass.setRange("left",Bmass.getMin(),left);
  Bmass.setRange("peak",left,right);

  //n values
  double n_signal_initial = data->sumEntries(TString::Format("abs(Bmass-%g)<0.05",mass_peak)) - data->sumEntries(TString::Format("abs(Bmass-%g)<0.10&&abs(Bmass-%g)>0.05",mass_peak,mass_peak));
  double n_combinatorial_initial = data->sumEntries() - n_signal_initial;
  RooRealVar n_signal("n_signal","n_signal",n_signal_initial,0.,data->sumEntries());
  RooRealVar n_combinatorial("n_combinatorial","n_combinatorial",n_combinatorial_initial,0.,data->sumEntries());

  RooRealVar f_erf("f_erf","f_erf",2.50259e-01,0,1);
  RooProduct n_erf("n_erf","n_erf",RooArgList(n_signal,f_erf));
  
  RooRealVar f_jpsipi("f_jpsipi","f_jpsipi",4.1E-5/1.026E-3,0.,0.1); 
  f_jpsipi.setConstant(kTRUE);
  RooProduct n_jpsipi("n_jpsipi","n_jpsipi",RooArgList(n_signal,f_jpsipi));

  if(particle == 0){
    RooAddPdf model("model", "model", RooArgList(signal,fit_side,erf,jpsipi),RooArgList(n_signal,n_combinatorial,n_erf,n_jpsipi));
    model.fitTo(*data,Range("all"));
    w.import(model);
  }else if(particle == 1){
    RooAddPdf model("model", "model", RooArgList(signal,fit_side), RooArgList(n_signal, n_combinatorial)); 
    model.fitTo(*data,Range("all"));
    w.import(model);
    {
  
      w.import(fit_side);
      w.import(signal);

    } 
  }
}

void plot_complete_fit(RooWorkspace& w){

  RooAbsPdf*  model = w.pdf("model");
  RooDataSet* data = (RooDataSet*) w.data("data");
  RooRealVar Bmass = *(w.var("Bmass"));
  RooRealVar* lambda   = w.var("lambda");

  RooPlot* massframe = Bmass.frame();

  data->plotOn(massframe, RooFit::Name("Data"));
  model->plotOn(massframe, RooFit::Name("Fit"),Range("all"),LineColor(kRed),LineStyle(1),LineWidth(2));
  model->plotOn(massframe, RooFit::Name("Combinatorial"),Components("fit_side"),Range("all"),LineColor(kBlue),LineStyle(kDashed));
  model->plotOn(massframe, RooFit::Name("Signal"),Components("signal"),Range("all"),LineColor(kOrange),LineStyle(kDashed));
  model->plotOn(massframe, RooFit::Name("B->J/psi X"),Components("erf"),Range("all"),LineColor(kGreen+3),LineStyle(kDashed));
  model->plotOn(massframe, RooFit::Name("B->J/psi pi"),Components("jpsipi"),Range("all"),LineColor(kPink+10),LineStyle(kDashed));
  model->paramOn(massframe,Layout(0.60,0.90,0.75));
  massframe->getAttText()->SetTextSize(0.028);
  massframe->GetYaxis()->SetTitleOffset(1.3);
  massframe->SetXTitle("Bmass (GeV)");

  TCanvas d;

  TPad *p1 = new TPad("p1","p1",0.0,0.27,0.82,0.99);
  p1->SetBorderMode(1); 
  p1->SetFrameBorderMode(0); 
  p1->SetBorderSize(2);

  p1->SetBottomMargin(0.10);

  p1->Draw(); 
     
  TPad *p2 = new TPad("p2","p2",0.0,0.065,0.82,0.24);
  p2->SetTopMargin(0.); 
  p2->SetBottomMargin(0.2);
   
  p2->SetBorderMode(1);
  p2->SetFrameBorderMode(0);
  p2->SetBorderSize(1); 
  
  p2->Draw();

  p1->cd();

  massframe->Draw();
  TLatex* tex11 = new TLatex(0.6,0.8,"1.5 nb^{-1} (PbPb) 5.02 TeV");
  tex11->SetNDC(kTRUE);
  tex11->SetLineWidth(2);
  tex11->SetTextSize(0.04);
  tex11->Draw();
  tex11 = new TLatex(0.6,0.85,"CMS Preliminary");
  tex11->SetNDC(kTRUE);
  tex11->SetTextFont(42);
  tex11->SetTextSize(0.04);
  tex11->SetLineWidth(2);
  tex11->Draw();

  double lambda_str = lambda->getVal();
  double lambda_err = lambda->getError();

  double chis = massframe->chiSquare();

  TLatex* tex12 = new TLatex(0.15, 0.85, Form("#lambda_{exp} = %.3lf #pm %.3lf",lambda_str,lambda_err));
  tex12->SetNDC(kTRUE);
  tex12->SetTextFont(42);
  tex12->SetTextSize(0.04);
 
  TLatex* tex13 = new TLatex(0.15, 0.8, Form("#chi/DOF = %.3lf",chis));
  tex13->SetNDC(kTRUE);
  tex13->SetTextFont(42);
  tex13->SetTextSize(0.04);

  TLegend *leg = new TLegend (0.4, 0.5, 0.6, 0.7);
  leg->SetTextSize(0.03);
  leg->AddEntry(massframe->findObject("Data"), "Data", "l");
  leg->AddEntry(massframe->findObject("B->J/psi X"), "B->J/psi X", "l");
  leg->AddEntry(massframe->findObject("Signal"), "Signal", "l");
  leg->AddEntry(massframe->findObject("Combinatorial"), "Combinatorial", "l");
  leg->AddEntry(massframe->findObject("B->J/psi pi"), "B->J/psi pi", "l");
  leg->AddEntry(massframe->findObject("Fit"),"Fit","l");
  leg->Draw("same");
  
  //pull dists

  RooHist* pull_hist = massframe->pullHist("Data","Fit");
  RooPlot *pull_plot = Bmass.frame();

  pull_plot->addPlotable(static_cast<RooPlotable*>(pull_hist),"P");
  pull_plot->SetTitle("");
  pull_plot->GetXaxis()->SetTitle("");
  pull_plot->GetXaxis()->SetTitleFont(42);
  pull_plot->GetXaxis()->SetTitleSize(0.17);
  pull_plot->GetXaxis()->SetTitleOffset(1.09);

  pull_plot->GetXaxis()->SetLabelFont(42);
  pull_plot->GetXaxis()->SetLabelSize(0.15);
  pull_plot->GetXaxis()->SetLabelOffset(0.01);
  pull_plot->GetXaxis()->SetTickLength(0.13);
  
  pull_plot->GetYaxis()->SetTitle("Pull hist");
  pull_plot->GetYaxis()->SetTitleFont(42);  
  pull_plot->GetYaxis()->SetTitleSize(0.10);
  pull_plot->GetYaxis()->SetTitleOffset(1.09);

  pull_plot->GetYaxis()->SetLabelFont(42);
  pull_plot->GetYaxis()->SetLabelSize(0.13);
  pull_plot->GetYaxis()->SetLabelOffset(0.005);
  
  pull_plot->GetYaxis()->SetNdivisions(305);

  p2->cd();
  pull_plot->Draw();
 
  if(particle == 0){
    d.SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/teste2/complete_fitB+.pdf"); 
  }else if(particle == 1){
    d.SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/Bs/complete_fitBs.pdf");
  }
}

//SIDEBAND SUBTRACTION//
std::vector<TH1D*> sideband_subtraction(RooWorkspace* w, int* n, int n_var){
  
  RooDataSet* data = (RooDataSet*) w->data("data");

  RooAbsPdf* fit_side = w->pdf("fit_side");

  vector<RooRealVar> variables;

  variables.push_back(*(w->var("Bmass")));
  variables.push_back(*(w->var("Bpt")));
  variables.push_back(*(w->var("By")));
  variables.push_back(*(w->var("Btrk1eta")));
  variables.push_back(*(w->var("Btrk1Y")));
  if(particle == 1){variables.push_back(*(w->var("Btrk2eta")));}
  variables.push_back(*(w->var("Btrk1pt")));
  if(particle == 1){variables.push_back(*(w->var("Btrk2pt")));}
  variables.push_back(*(w->var("Bmu1eta")));
  variables.push_back(*(w->var("Bmu2eta")));
  variables.push_back(*(w->var("Bmu1pt")));
  variables.push_back(*(w->var("Bmu2pt")));
  variables.push_back(*(w->var("Bchi2cl")));
  if(particle == 1){variables.push_back(*(w->var("Bmumumass")));}
  if(particle == 1){variables.push_back(*(w->var("Btrktrkmass")));}
  variables.push_back(*(w->var("BsvpvDistance")));
  variables.push_back(*(w->var("BsvpvDistance_Err")));
  variables.push_back(*(w->var("Balpha")));
  if(particle == 0){variables.push_back(*(w->var("Btrk1Dz1")));}
  variables.push_back(*(w->var("BvtxX")));
  variables.push_back(*(w->var("BvtxY")));
  if(particle == 0){
    variables.push_back(*(w->var("Btrk1DzError1")));
    variables.push_back(*(w->var("Btrk1Dxy1")));
    variables.push_back(*(w->var("Btrk1DxyError1")));
    variables.push_back(*(w->var("Bd0")));
    variables.push_back(*(w->var("Bd0err")));
  }

  RooDataSet* reduceddata_side;
  RooDataSet* reduceddata_central; 

  double left = 5.15;
  double right = 5.4;

  reduceddata_side = (RooDataSet*) data->reduce(Form("Bmass>%lf",right));
  reduceddata_central = (RooDataSet*) data->reduce(Form("Bmass>%lf",left));
  reduceddata_central = (RooDataSet*) reduceddata_central->reduce(Form("Bmass<%lf",right));
 
  //Integrating the background distribution

  RooAbsReal* int_fit_side_right = fit_side->createIntegral(variables[0], "right");
  RooAbsReal* int_fit_peak = fit_side->createIntegral(variables[0], "peak");

  std::cout<< std::endl << "Integral right band: " << int_fit_side_right->getVal() << std::endl;

  double factor = (int_fit_peak->getVal())/(int_fit_side_right->getVal());

  std::cout << std::endl << "Factor: " << factor << std::endl;
  for(int i=0; i<n_var; i++){
    std::cout << "bins: " << n[i] << std::endl;
  } 

  std::vector<TH1D*> histos;

  if(particle == 0){
    histos.push_back(create_histogram(variables[1],"Bpt", factor, reduceddata_side, reduceddata_central, data, n[0]));
    histos.push_back(create_histogram(variables[2], "By",factor, reduceddata_side, reduceddata_central, data, n[1]));
    histos.push_back(create_histogram(variables[3], "Btrk1eta",factor, reduceddata_side, reduceddata_central, data, n[2]));
    histos.push_back(create_histogram(variables[4], "Btrk1Y",factor, reduceddata_side, reduceddata_central, data, n[3]));
    histos.push_back(create_histogram(variables[5], "Btrk1pt",factor, reduceddata_side, reduceddata_central, data, n[4]));
    histos.push_back(create_histogram(variables[6], "Bmu1eta",factor, reduceddata_side, reduceddata_central, data, n[5]));
    histos.push_back(create_histogram(variables[7], "Bmu2eta",factor, reduceddata_side, reduceddata_central, data, n[6]));
    histos.push_back(create_histogram(variables[8], "Bmu1pt",factor, reduceddata_side, reduceddata_central, data, n[7]));
    histos.push_back(create_histogram(variables[9], "Bmu2pt",factor, reduceddata_side, reduceddata_central, data, n[8]));
    histos.push_back(create_histogram(variables[10], "Bchi2cl",factor, reduceddata_side, reduceddata_central, data, n[9]));
    histos.push_back(create_histogram(variables[11], "BsvpvDistance",factor, reduceddata_side, reduceddata_central, data, n[10]));
    histos.push_back(create_histogram(variables[12], "BsvpvDistance_Err",factor, reduceddata_side, reduceddata_central, data, n[11]));
    histos.push_back(create_histogram(variables[13], "Balpha",factor, reduceddata_side, reduceddata_central, data, n[12]));
    histos.push_back(create_histogram(variables[14], "Btrk1Dz1",factor, reduceddata_side, reduceddata_central, data, n[13]));
    histos.push_back(create_histogram(variables[15], "BvtxX",factor, reduceddata_side, reduceddata_central, data, n[14]));
    histos.push_back(create_histogram(variables[16], "BvtxY",factor, reduceddata_side, reduceddata_central, data, n[15]));
    histos.push_back(create_histogram(variables[17], "Btrk1DzError1",factor, reduceddata_side, reduceddata_central, data, n[16]));
    histos.push_back(create_histogram(variables[18], "Btrk1Dxy1",factor, reduceddata_side, reduceddata_central, data, n[17]));
    histos.push_back(create_histogram(variables[19], "Btrk1DxyError1",factor, reduceddata_side, reduceddata_central, data, n[18]));
    histos.push_back(create_histogram(variables[20], "Bd0",factor, reduceddata_side, reduceddata_central, data, n[19]));
    histos.push_back(create_histogram(variables[21], "Bd0err",factor, reduceddata_side, reduceddata_central, data, n[20]));
  }else if(particle == 1){
    histos.push_back(create_histogram(variables[1],"Bpt", factor, reduceddata_side, reduceddata_central, data, n[0]));
    histos.push_back(create_histogram(variables[2], "By",factor, reduceddata_side, reduceddata_central, data, n[1]));
    histos.push_back(create_histogram(variables[3], "Btrk1eta",factor, reduceddata_side, reduceddata_central, data, n[2]));
    histos.push_back(create_histogram(variables[4], "Btrk1Y",factor, reduceddata_side, reduceddata_central, data, n[3]));
    histos.push_back(create_histogram(variables[5], "Btrk2eta",factor, reduceddata_side, reduceddata_central, data, n[4]));
    histos.push_back(create_histogram(variables[5], "Btrk1pt",factor, reduceddata_side, reduceddata_central, data, n[4]));
    histos.push_back(create_histogram(variables[5], "Btrk2pt",factor, reduceddata_side, reduceddata_central, data, n[4]));
    histos.push_back(create_histogram(variables[6], "Bmu1eta",factor, reduceddata_side, reduceddata_central, data, n[5]));
    histos.push_back(create_histogram(variables[7], "Bmu2eta",factor, reduceddata_side, reduceddata_central, data, n[6]));
    histos.push_back(create_histogram(variables[8], "Bmu1pt",factor, reduceddata_side, reduceddata_central, data, n[7]));
    histos.push_back(create_histogram(variables[9], "Bmu2pt",factor, reduceddata_side, reduceddata_central, data, n[8]));
    histos.push_back(create_histogram(variables[10], "Bchi2cl",factor, reduceddata_side, reduceddata_central, data, n[9]));
    histos.push_back(create_histogram(variables[5], "Bmumumass",factor, reduceddata_side, reduceddata_central, data, n[4]));
    histos.push_back(create_histogram(variables[5], "Btrktrkmass",factor, reduceddata_side, reduceddata_central, data, n[4]));
    histos.push_back(create_histogram(variables[11], "BsvpvDistance",factor, reduceddata_side, reduceddata_central, data, n[10]));
    histos.push_back(create_histogram(variables[12], "BsvpvDistance_Err",factor, reduceddata_side, reduceddata_central, data, n[11]));
    histos.push_back(create_histogram(variables[13], "Balpha",factor, reduceddata_side, reduceddata_central, data, n[12]));
    histos.push_back(create_histogram(variables[15], "BvtxX",factor, reduceddata_side, reduceddata_central, data, n[14]));
    histos.push_back(create_histogram(variables[16], "BvtxY",factor, reduceddata_side, reduceddata_central, data, n[15]));
  }

  return histos;
  //data histograms
}
//sideband subtraction ends


TH1D* create_histogram_mc(RooRealVar var, TTree* t, int n){


  TH1D* h = new TH1D(var.GetName(), var.GetName(), n, var.getMin(), var.getMax());

  TString name_string = TString(var.GetName()) + ">>htemp(" + Form("%d",n) +"," + Form("%lf", var.getMin()) + "," + Form("%lf", var.getMax()) + ")";

  t->Draw(name_string, "Pthatweight");

  h = (TH1D*)gDirectory->Get("htemp")->Clone();
  h->SetTitle("");
  h->SetMarkerColor(kGreen);
  h->SetLineColor(kGreen);
  return h;
}
//create_histogram_mc ends

TH1D* create_histogram(RooRealVar var,TString name, double factor, RooDataSet* reduced, RooDataSet* central, RooDataSet* total, int n){


  std::cout<< "n in create_histogram = "<< n <<std::endl;
  TH1D* dist_side = (TH1D*)reduced->createHistogram("dist_side",var, Binning(n, var.getMin(), var.getMax()));
  dist_side->SetMarkerColor(kRed);
  dist_side->SetLineColor(kRed);
  dist_side->SetNameTitle("dist_side", "");

  TH1D* hist_dist_peak = (TH1D*)central->createHistogram(var.GetName(), var, Binning(n, var.getMin(), var.getMax()));
  TH1D* dist_peak = new TH1D(*hist_dist_peak);
  dist_peak->SetMarkerColor(kBlue);
  dist_peak->SetLineColor(kBlue);
  dist_peak->SetNameTitle(var.GetName(), "");

  hist_dist_peak->SetMarkerColor(kBlack);
  hist_dist_peak->SetLineColor(kBlack);
  hist_dist_peak->SetNameTitle("dist_total", "");

  dist_peak->Add(dist_side, -factor);
  dist_side->Scale(factor);

  dist_peak->SetStats(0);
  dist_side->SetStats(0);
  hist_dist_peak->SetStats(0);
  TCanvas c;

  hist_dist_peak->Draw();
  dist_side->Draw("same");
  dist_peak->Draw("same");
  
  dist_peak->SetXTitle(var.GetName());
  dist_side->SetXTitle(var.GetName());
  hist_dist_peak->SetXTitle(var.GetName());

  hist_dist_peak->GetYaxis()->SetRangeUser(0, 1.3*hist_dist_peak->GetMaximum());
  TLatex* tex = new TLatex(0.6,0.8,"1.5 nb^{-1} (PbPb) 5.02 TeV");
  tex->SetNDC(kTRUE);
  tex->SetLineWidth(2);
  tex->SetTextSize(0.04);
  tex->Draw();
  tex = new TLatex(0.68,0.85,"CMS Preliminary");
  tex->SetNDC(kTRUE);
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();

  TLegend *leg = new TLegend (0.7, 0.5, 0.85, 0.65);
  leg->AddEntry(var.GetName(), "Signal", "l");
  leg->AddEntry("dist_side", "Background", "l");
  leg->AddEntry("hist_dist_peak", "Sinal inicial", "l");
  leg->Draw("same");

  std::cout<<"name: "<<var.GetName()<<std::endl;
  std::cout<<"histo name: "<<dist_peak->GetName()<<std::endl;

  if(particle == 0){
    c.SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/teste2/sideband_sub/"+name + "sideband_subB+.pdf");
  }else if(particle == 1){
    c.SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/Bs/sideband_sub/"+name + "sideband_subBs.pdf");
  }
  return dist_peak;

}
//create_histogram ends


void do_splot(RooWorkspace& w){

   RooDataSet* data = (RooDataSet*) w.data("data");   
   RooAbsPdf* model = w.pdf("model");
   //we need the fit and the dataset previously saved in the woorkspace

   RooRealVar* BpYield = w.var("n_signal");
   RooRealVar* BgYield = w.var("n_combinatorial");
   //we need the n values previously saved in the woorkspace

   //fit the model to the data
   model->fitTo(*data,Extended());

   //sPlot technique requires model parameters (other than the yields) to be fixed

   RooRealVar* mean  = w.var("mean");
   RooRealVar* sigma1 = w.var("sigma1");
   RooRealVar* sigma2 = w.var("sigma2");
   RooRealVar* cofs = w.var("cofs");
   RooRealVar* lambda = w.var("lambda");

   mean->setConstant();
   sigma1->setConstant();
   sigma2->setConstant();
   cofs->setConstant();
   lambda->setConstant();

   RooMsgService::instance().setSilentMode(true);

  //add sWeights to dataset based on model and yield variables
  //sPlot class adds a new variable that has the name of the corresponding yield + "_sw".
   SPlot* sData = new SPlot("sData","An sPlot",*data, model, RooArgList(*BpYield,*BgYield));

  //An SPlot gives us the distribution of some variable, x in our data sample for a given species (eg. signal or background). The result is similar to a likelihood projection plot, but no cuts are made, so every event contributes to the distribution.


  cout << endl <<  "Yield of B+ is "
       << BpYield->getVal() << ".  From sWeights it is "
       << sData->GetYieldFromSWeight("n_signal") << endl;

  cout << "Yield of background is "
       << BgYield->getVal() << ".  From sWeights it is "
       << sData->GetYieldFromSWeight("n_combinatorial") << endl
       << endl;

  for(Int_t i=0; i < 10; i++) {
    if(0)
      cout << "y Weight   "     << sData->GetSWeight(i,"BpYield")
	   << "\tb Weight   "     << sData->GetSWeight(i,"BgYield")
	   << "\ttotal Weight   " << sData->GetSumOfEventSWeight(i)
	   << endl;
  }

  cout << endl;

  w.import(*data, Rename("dataWithSWeights"));
  //the reweighted data is saved in the woorkspace
 
 }


//do_splot ends

TH1D* make_splot(RooWorkspace& w, int n, TString label){

  //saves the plots of signal distributions, background distributions and signal+background distributions
  //in the end returns the histogram of signal
 

  TCanvas* cdata = new TCanvas("sPlot","sPlot", 800, 600);
  cdata->Divide(2,2);

  RooAbsPdf* model  = w.pdf("model");
  RooAbsPdf* BpModel = w.pdf("signal");
  RooAbsPdf* BgModel = w.pdf("fit_side");

  RooRealVar* Bmass  = w.var("Bmass");
  RooRealVar* variable = w.var(label);

  RooRealVar* BpYield = w.var("n_signal");
  RooRealVar* BgYield = w.var("n_combinatorial");

  double sigYield = BpYield->getVal();
  double bkgYield = BgYield->getVal();


  RooDataSet* data = (RooDataSet*) w.data("data");

  cdata->cd(1);
  RooPlot* mframe = Bmass->frame();
  if(particle == 0){
    mframe->GetXaxis()->SetTitle(TString::Format("mass of B+ [GeV]"));
  }else if(particle == 1){
    mframe->GetXaxis()->SetTitle(TString::Format("mass of Bs [GeV]"));
  }
  data->plotOn(mframe);
  model->plotOn(mframe,LineColor(kRed));
  model->plotOn(mframe,Components(*BpModel),LineStyle(kDashed),LineColor(kOrange));
  model->plotOn(mframe,Components(*BgModel),LineStyle(kDashed),LineColor(kGreen+3));
  mframe->SetTitle("Bmass");
  mframe->Draw();

  cdata->cd(2);
  RooPlot* ptframe = variable->frame();
  data->plotOn(ptframe);
  if(particle == 0){
    ptframe->SetTitle(label + " of B+: total sample");
  }else if(particle == 1){
    ptframe->SetTitle(label + " of Bs: total sample");
  }
  ptframe->Draw();

  //get the dataset with sWeights
  RooDataSet* dataW = (RooDataSet*) w.data("dataWithSWeights");
  RooDataSet* dataWBp = new RooDataSet(dataW->GetName(),dataW->GetTitle(),dataW,*dataW->get(),0,"n_signal_sw");
  RooDataSet* dataWBg = new RooDataSet(dataW->GetName(),dataW->GetTitle(),dataW,*dataW->get(),0,"n_combinatorial_sw");

  RooPlot* ptframe2Bp = variable->frame();
  RooPlot* ptframe2Bg = variable->frame();

  ptframe2Bp->GetYaxis()->SetTitle(TString::Format("Events / (%g)",(variable->getMax()-variable->getMin())/n));
  ptframe2Bg->GetYaxis()->SetTitle(TString::Format("Events / (%g)",(variable->getMax()-variable->getMin())/n));

  if(particle == 0){
    ptframe2Bp->GetXaxis()->SetTitle(label + " of B+");
    ptframe2Bg->GetXaxis()->SetTitle(label + " of B+");
  }else if(particle == 1){
    ptframe2Bp->GetXaxis()->SetTitle(label + " of Bs");
    ptframe2Bg->GetXaxis()->SetTitle(label + " of Bs");
  }

  dataWBp->plotOn(ptframe2Bp, DataError(RooAbsData::SumW2),Binning(n));
  dataWBg->plotOn(ptframe2Bg, DataError(RooAbsData::SumW2),Binning(n));

  if(particle == 0){
    ptframe2Bp->SetTitle(label+" distribution of B+ for signal (splot)");
    ptframe2Bg->SetTitle(label+" distribution of B+ for background (splot)");
  }else if(particle == 1){
    ptframe2Bp->SetTitle(label+" distribution of Bs for signal (splot)");
    ptframe2Bg->SetTitle(label+" distribution of Bs for background (splot)");
  }

  cdata->cd(3);  ptframe2Bp->Draw();
  cdata->cd(4);  ptframe2Bg->Draw();

  if(particle == 0){
    cdata->SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/teste2/splot/Bmass/"+label+"sPlotB+.gif");
    cdata->SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/teste2/splot/Bmass/"+label+"sPlotB+.pdf");
  }else if(particle == 1){
    cdata->SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/Bs/splot/Bmass/"+label+"sPlotBs.gif");
    cdata->SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/Bs/splot/Bmass/"+label+"sPlotBs.pdf");
  }

  TH1D* histo_Bp_sig = (TH1D*)dataWBp->createHistogram(label,n,0,0);
  TH1D* histo_Bp_bkg = (TH1D*)dataWBg->createHistogram(label,n,0,0);

    for (int i=1; i<=n; i++) {

      if (histo_Bp_sig->GetBinContent(i)==0) histo_Bp_sig->SetBinError(i,0.);
      if (histo_Bp_bkg->GetBinContent(i)==0) histo_Bp_bkg->SetBinError(i,0.);

       histo_Bp_sig->SetBinContent(i,histo_Bp_sig->GetBinContent(i)/sigYield);
       histo_Bp_sig->SetBinError(i,histo_Bp_sig->GetBinError(i)/sigYield);

       histo_Bp_bkg->SetBinContent(i,histo_Bp_bkg->GetBinContent(i)/bkgYield);
       histo_Bp_bkg->SetBinError(i,histo_Bp_bkg->GetBinError(i)/bkgYield);

  }

  TCanvas* prov = new TCanvas ("prov","c1",200,10,700,500);
  prov->cd();
  histo_Bp_sig->SetMarkerStyle(20);
  histo_Bp_sig->SetMarkerSize(0.);
  histo_Bp_sig->SetMarkerColor(kRed);
  histo_Bp_sig->SetLineColor(kRed);
  histo_Bp_sig->SetTitle("");
  histo_Bp_sig->GetYaxis()->SetTitle(TString::Format("Events /(%g)",(variable->getMax()-variable->getMin())/n));
  histo_Bp_sig->GetXaxis()->SetTitle(label );

  histo_Bp_sig->SetStats(0);
  histo_Bp_sig->Draw("E");

  if(particle == 0){
    prov->SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/teste2/splot/sig/"+label+"sPlotB+.gif");
    //prov->SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/teste2/splot/sig/"+label+"sPlotB+.pdl"); 
    prov->SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/teste2/splot/sig/"+label+"sPlotB+.pdf");
  }else if(particle == 1){
    prov->SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/Bs/splot/sig/"+label+"sPlotBs.gif");
    //prov->SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/teste2/splot/sig/"+label+"sPlotBs.pdl"); 
    prov->SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/Bs/splot/sig/"+label+"sPlotBs.pdf");
  }

  TCanvas* prov_bkg = new TCanvas ("prov_bkg","c2",200,10,700,500);
  prov_bkg->cd();
  histo_Bp_bkg->SetMarkerStyle(20);
  histo_Bp_bkg->SetMarkerSize(0.);
  histo_Bp_bkg->SetMarkerColor(kBlue);
  histo_Bp_bkg->SetLineColor(kBlue);
  histo_Bp_bkg->SetTitle("");
  histo_Bp_bkg->GetYaxis()->SetTitle(TString::Format("Events /(%g)",(variable->getMax()-variable->getMin())/n));
  histo_Bp_bkg->GetXaxis()->SetTitle(label);

  histo_Bp_bkg->SetStats(0);
  histo_Bp_bkg->Draw("E");

  if(particle == 0){
    prov_bkg->SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/teste2/splot/bkg/"+label+"sPlotB+.gif");
    prov_bkg->SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/teste2/splot/bkg/"+label+"sPlotB+.pdf");
  }else if(particle == 1){
    prov_bkg->SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/Bs/splot/bkg/"+label+"sPlotBs.gif");
    prov_bkg->SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/Bs/splot/bkg/"+label+"sPlotBs.pdf");
  }

  TCanvas* sig_bkg = new TCanvas ("sig_bkg","c3",200,10,700,500); 
  sig_bkg->cd();

  histo_Bp_sig->Draw();
  histo_Bp_bkg->Draw("SAME");

  TLegend* legend = new TLegend(0.7,0.9,0.9,0.8);
   legend->AddEntry(histo_Bp_sig,"Signal","lep");
   legend->AddEntry(histo_Bp_bkg,"Background","lep");
   legend->Draw();

   if(particle == 0){
     sig_bkg->SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/teste2/splot/sig_bkg/"+label+"sPlotB+.gif");
     sig_bkg->SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/teste2/splot/sig_bkg/"+label+"sPlotB+.pdf");
   }else if(particle == 1){
     sig_bkg->SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/Bs/splot/sig_bkg/"+label+"sPlotBs.gif");
     sig_bkg->SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/Bs/splot/sig_bkg/"+label+"sPlotBs.pdf");
   }

  //cleanup
  delete cdata;
  delete prov;
  delete prov_bkg;
  delete sig_bkg;

  
  return histo_Bp_sig;

} 
//make_splot ends

//SPLOT_METHOD//
std::vector<TH1D*> splot_method(RooWorkspace& w, int* n, TString* label, int n_var){

  std::vector<TH1D*> histos;

  for(int i = 0;i<n_var;i++){
    histos.push_back(make_splot(w,n[i],label[i]));
  }

  return histos;
}

void validate_fit(RooWorkspace* w)
{
  /*
  RooRealVar Bmass = *(w->var("Bmass"));
  RooRealVar n_signal = *(w->var("n_signal"));
  RooAbsPdf& model  = *(w->pdf("model"));

  vector<RooRealVar> params;
  params.push_back(n_signal);

  RooMCStudy* mcstudy = new RooMCStudy(model, Bmass, Binned(kTRUE), Silence(), Extended(), FitOptions(Save(kTRUE), PrintEvalErrors(0)));

  mcstudy->generateAndFit(1000);

  vector<RooPlot*> framesPull, framesParam;

  for(int i = 0; i < params.size(); ++i)
    {
      framesPull.push_back(mcstudy->plotPull(params.at(i),FrameBins(200),FrameRange(-5,5)));
      framesPull[i]-SetTitle("");
      framesParam.push_back(mcstudy->plotParam(params.at(i),FrameBins(50)));
      framesParam[i]->SetTitle("");
    }

  vector<TGraph*> h;

  for

  gStyle->SetOptFit(0111);

  TCanvas* c_pull = new TCanvas("pulls", "pulls", 900, 800);

  gPad->SetLeftMargin(0.15);

  c_pull->cd();
  framePULL->SetTitle("");
  framePULL->Draw();
  c_pull->Update();
  framePULL->Fit("gaus","","",-5,5);
  framePULL->GetFunction("gaus")->SetLineColor(4);
  framePULL->GetFunction("gaus")->SetLineWidth(5);
  framePULL->GetXaxis()->SetTitle("Pull");
  framePULL->Draw("same");
  */
}

void set_up_workspace_variables(RooWorkspace& w)
{

  if(particle == 0){

     double mass_min, mass_max;
     double pt_min, pt_max;
     double y_min, y_max;
     double trk1eta_min, trk1eta_max;
     double trk1Y_min, trk1Y_max;
     double trk1pt_min, trk1pt_max;
     double mu1eta_min, mu1eta_max;
     double mu2eta_min, mu2eta_max;
     double mu1pt_min, mu1pt_max;
     double mu2pt_min, mu2pt_max;
     double chi2cl_min, chi2cl_max;
     double svpvDistance_min, svpvDistance_max;
     double svpvDistanceErr_min, svpvDistanceErr_max;
     double alpha_min, alpha_max;
     double trk1Dz_min, trk1Dz_max;
     double vtxX_min, vtxX_max;
     double vtxY_min, vtxY_max;
     double trk1DzError1_min, trk1DzError1_max;
     double trk1Dxy1_min, trk1Dxy1_max;
     double trk1DxyErr1_min, trk1DxyErr1_max;
     double d0_min, d0_max;
     double d0Err_min, d0Err_max;
  
     mass_min=5.;
     mass_max=6.;

     pt_min=5.;
     pt_max=100.;

     y_min=-2.4;
     y_max=2.4;

     trk1eta_min=-2.5;
     trk1eta_max=2.5;

     trk1Y_min = -2.5;
     trk1Y_max = 2.5;

     trk1pt_min = 0.;
     trk1pt_max = 25.;

     mu1eta_min=-2.5;
     mu1eta_max=2.5;

     mu2eta_min = -2.6;
     mu2eta_max = 2.6;

     mu1pt_min=0.;
     mu1pt_max=52.;

     mu2pt_min = 0.;
     mu2pt_max = 54.;

     chi2cl_min = 0.;
     chi2cl_max = 1.05;

     svpvDistance_min=0.;
     svpvDistance_max=9.5;

     svpvDistanceErr_min=0.;
     svpvDistanceErr_max=0.08;

     alpha_min=0.;
     alpha_max=0.1;

     trk1Dz_min=-9;
     trk1Dz_max=5.;

     vtxX_min = -0.85;
     vtxX_max = 0.8;

     vtxY_min = -0.9;
     vtxY_max = 0.9;

     trk1DzError1_min = 0;
     trk1DzError1_max = 1.5;

     trk1Dxy1_min = -0.45;
     trk1Dxy1_max = 0.35;

     trk1DxyErr1_min = 0;
     trk1DxyErr1_max = 0.2;

     d0_min=0.;
     d0_max=0.9;

     d0Err_min=0.;
     d0Err_max=0.00033;

     RooRealVar Bmass("Bmass","Bmass",mass_min,mass_max);
     RooRealVar Bpt("Bpt","Bpt",pt_min,pt_max);
     RooRealVar By("By","By",y_min,y_max);
     RooRealVar Btrk1eta("Btrk1eta","Btrk1eta",trk1eta_min,trk1eta_max);
     RooRealVar Btrk1Y("Btrk1Y","Btrk1Y",trk1Y_min,trk1Y_max);
     RooRealVar Btrk1pt("Btrk1pt","Btrk1pt",trk1pt_min,trk1pt_max);
     RooRealVar Bmu1eta("Bmu1eta","Bmu1eta",mu1eta_min,mu1eta_max);
     RooRealVar Bmu2eta("Bmu2eta","Bmu2eta",mu2eta_min,mu2eta_max);
     RooRealVar Bmu1pt("Bmu1pt","Bmu1pt",mu1pt_min,mu1pt_max);
     RooRealVar Bmu2pt("Bmu2pt","Bmu2pt",mu2pt_min,mu2pt_max);
     RooRealVar Bchi2cl("Bchi2cl", "Bchi2cl", chi2cl_min, chi2cl_max);
     RooRealVar BsvpvDistance("BsvpvDistance", "BsvpvDistance", svpvDistance_min, svpvDistance_max);
     RooRealVar BsvpvDistance_Err("BsvpvDistance_Err", "BsvpvDistance_Err", svpvDistanceErr_min, svpvDistanceErr_max);
     RooRealVar Balpha("Balpha", "Balpha", alpha_min, alpha_max);
     RooRealVar Btrk1Dz1("Btrk1Dz1","Btrk1Dz1",trk1Dz_min,trk1Dz_max);
     RooRealVar BvtxX("BvtxX","BvtxX", vtxX_min,vtxX_max);
     RooRealVar BvtxY("BvtxY","BvtxY",vtxY_min,vtxY_max);
     RooRealVar Btrk1DzError1("Btrk1DzError1","Btrk1DzError1",trk1DzError1_min,trk1DzError1_max);
     RooRealVar Btrk1Dxy1("Btrk1Dxy1","Btrk1Dxy1",trk1Dxy1_min,trk1Dxy1_max);
     RooRealVar Btrk1DxyError1("Btrk1DxyError1","Btrk1DxyError1",trk1DxyErr1_min,trk1DxyErr1_max);
     RooRealVar Bd0("Bd0", "Bd0", d0_min, d0_max);
     RooRealVar Bd0err("Bd0err", "Bd0err", d0Err_min, d0Err_max);
 
     w.import(Bmass);
     w.import(Bpt);
     w.import(By);
     w.import(Btrk1eta);
     w.import(Btrk1Y);
     w.import(Btrk1pt);
     w.import(Bmu1eta);
     w.import(Bmu2eta);
     w.import(Bmu1pt);
     w.import(Bmu2pt);
     w.import(Bchi2cl);
     w.import(BsvpvDistance);
     w.import(BsvpvDistance_Err);
     w.import(Balpha);
     w.import(Btrk1Dz1);
     w.import(BvtxX);
     w.import(BvtxY);
     w.import(Btrk1DzError1);
     w.import(Btrk1Dxy1);
     w.import(Btrk1DxyError1);
     w.import(Bd0);
     w.import(Bd0err);
  }
   
  
  if(particle == 1){
    double mass_min, mass_max;
    double pt_min, pt_max;
    double y_min, y_max;
    double trk1eta_min, trk1eta_max;
    double trk1Y_min, trk1Y_max;
    double trk2eta_min, trk2eta_max; 
    double trk1pt_min, trk1pt_max;
    double trk2pt_min, trk2pt_max;
    double mu1eta_min, mu1eta_max;
    double mu2eta_min, mu2eta_max;
    double mu1pt_min, mu1pt_max;
    double mu2pt_min, mu2pt_max;
    double chi2cl_min, chi2cl_max;
    double mumumass_min, mumumass_max;
    double trktrkmass_min, trktrkmass_max;
    double svpvDistance_min, svpvDistance_max;
    double svpvDistanceErr_min, svpvDistanceErr_max;
    double alpha_min, alpha_max;
    double vtxX_min, vtxX_max;
    double vtxY_min, vtxY_max;

    ///ESTES LIMITES ESTÃO TODOS MAL
  
    mass_min=5.;
    mass_max=6.;

    pt_min=5.;
    pt_max=100.;

    y_min=-2.4;
    y_max=2.4;

    trk1eta_min=-2.5;
    trk1eta_max=2.5;

    trk1Y_min = -2.5;
    trk1Y_max = 2.5;

    trk2eta_min = 0.;
    trk2eta_max = 0.;

    trk1pt_min = 0.;
    trk1pt_max = 25.;

    trk2pt_min = 0.;
    trk2pt_max = 0.;

    mu1eta_min=-2.5;
    mu1eta_max=2.5;

    mu2eta_min = -2.6;
    mu2eta_max = 2.6;

    mu1pt_min=0.;
    mu1pt_max=52.;

    mu2pt_min = 0.;
    mu2pt_max = 54.;

    chi2cl_min = 0.;
    chi2cl_max = 1.05;

    mumumass_min = 0.;
    mumumass_max = 0.;
     
    trktrkmass_min = 0.;
    trktrkmass_max = 0.;

    svpvDistance_min=0.;
    svpvDistance_max=9.5;

    svpvDistanceErr_min=0.;
    svpvDistanceErr_max=0.08;

    alpha_min=0.;
    alpha_max=0.1;

    vtxX_min = -0.85;
    vtxX_max = 0.8;

    vtxY_min = -0.9;
    vtxY_max = 0.9;

    RooRealVar Bmass("Bmass","Bmass",mass_min,mass_max);
    RooRealVar Bpt("Bpt","Bpt",pt_min,pt_max);
    RooRealVar By("By","By",y_min,y_max);
    RooRealVar Btrk1eta("Btrk1eta","Btrk1eta",trk1eta_min,trk1eta_max);
    RooRealVar Btrk1Y("Btrk1Y","Btrk1Y",trk1Y_min,trk1Y_max);    
    RooRealVar Btrk2eta("Btrk2eta","Btrk2eta",trk2eta_min,trk2eta_max);
    RooRealVar Btrk1pt("Btrk1pt","Btrk1pt",trk1pt_min,trk1pt_max);
    RooRealVar Btrk2pt("Btrk2pt","Btrk2pt",trk2pt_min, trk2pt_max);
    RooRealVar Bmu1eta("Bmu1eta","Bmu1eta",mu1eta_min,mu1eta_max);
    RooRealVar Bmu2eta("Bmu2eta","Bmu2eta",mu2eta_min,mu2eta_max);
    RooRealVar Bmu1pt("Bmu1pt","Bmu1pt",mu1pt_min,mu1pt_max);
    RooRealVar Bmu2pt("Bmu2pt","Bmu2pt",mu2pt_min,mu2pt_max);
    RooRealVar Bchi2cl("Bchi2cl", "Bchi2cl", chi2cl_min, chi2cl_max);
    RooRealVar Bmumumass("Bmumumass","Bmumumass",mumumass_min,mumumass_max);
    RooRealVar Btrktrkmass("Btrktrkmass","Btrktrkmass",trktrkmass_min,trktrkmass_max);
    RooRealVar BsvpvDistance("BsvpvDistance", "BsvpvDistance", svpvDistance_min, svpvDistance_max);
    RooRealVar BsvpvDistance_Err("BsvpvDistance_Err", "BsvpvDistance_Err", svpvDistanceErr_min, svpvDistanceErr_max);
    RooRealVar Balpha("Balpha", "Balpha", alpha_min, alpha_max);
    RooRealVar BvtxX("BvtxX","BvtxX", vtxX_min,vtxX_max);
    RooRealVar BvtxY("BvtxY","BvtxY",vtxY_min,vtxY_max);
 
    w.import(Bmass);
    w.import(Bpt);
    w.import(By);
    w.import(Btrk1eta);
    w.import(Btrk1Y);
    w.import(Btrk2eta);
    w.import(Btrk1pt);
    w.import(Btrk2pt);
    w.import(Bmu1eta);
    w.import(Bmu2eta);
    w.import(Bmu1pt);
    w.import(Bmu2pt);
    w.import(Bchi2cl);
    w.import(Bmumumass);
    w.import(Btrktrkmass);
    w.import(BsvpvDistance);
    w.import(BsvpvDistance_Err);
    w.import(Balpha);
    w.import(BvtxX);
    w.import(BvtxY);
  }
}
