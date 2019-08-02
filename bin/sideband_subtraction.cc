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
//tens que acrescentar:
#include <RooBifurGauss.h>
#include <RooProduct.h>
#include <RooHist.h>

#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "TCanvas.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"
#include "TFile.h"
#include "TNtupleD.h"
#include "TH1D.h"
#include "TLegend.h"

using namespace RooStats;
using namespace RooFit;
using namespace std;

void DoSPlot(RooWorkspace&);
vector<TH1D*> GetSPlot(RooWorkspace&, int, TString);
//void MakePlots(RooWorkspace&, int, TString);
void read_data(RooWorkspace&, TString, TString);
void set_up_workspace_variables(RooWorkspace& w);
TH1D* create_histogram(RooRealVar var, TTree* t, int n);
TH1D* create_histogram(RooRealVar var,TString name, double factor, RooDataSet* reduced, RooDataSet* central, RooDataSet* total, int n);
std::vector<TH1D*> sideband_subtraction(RooWorkspace* w, TString input_file_data, int* n);

int main(){

  const int n_var = 16;
  TString input_file_data = "/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/prefiltered_trees/selected_data_ntKp_PbPb_2018.root";
  TString input_file_mc = "/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/prefiltered_trees/selected_mc_ntKp_PbPb_2018_pthatweight.root";

  std::vector<TH1D*> histos_data;
  std::vector<TH1D*> histos_mc;
  int n_bins[]= {10, 20, 10, 10, 10, 10, 10, 10, 15, 10, 10, 10, 10, 10, 10, 10, 10, 10};
  TString variables[]={"Bpt","By","Btrk1D0Err","Bmu1pt","Bmu1eta","Btrk1pt","Btrk1eta","Bchi2cl","BsvpvDistance","BsvpvDistance_Err","Balpha","Btrk1D0","Btrk1Dz","Bd0","Blxy","Bd0err"};

  RooWorkspace* ws = new RooWorkspace("ws");

  set_up_workspace_variables(*ws);
  histos_data = sideband_subtraction(ws, input_file_data, n_bins);

  TFile *fin_mc = new TFile(input_file_mc);
  TTree* t1_mc = (TTree*)fin_mc->Get("ntKp");

  std::vector<TString> names;

  for(int i=0; i<(int)histos_data.size(); i++){
    std::cout<< "Var names: "<< histos_data[i]->GetName()<<std::endl;
  }
  for(int i=0; i<(int)histos_data.size(); i++){histos_mc.push_back(create_histogram((*ws->var(histos_data[i]->GetName())), t1_mc, n_bins[i]));
    names.push_back(TString(histos_data[i]->GetName()));
  } 

  for(int i=0; i<(int)histos_data.size(); i++)
    {
      TCanvas c;
      histos_mc[i]->SetXTitle(TString(histos_data[i]->GetName()));
      histos_mc[i]->SetStats(0);
      histos_data[i]->SetStats(0);
      histos_mc[i]->Scale(1/histos_mc[i]->Integral());
      histos_data[i]->Scale(1/histos_data[i]->Integral());
      histos_mc[i]->GetYaxis()->SetRangeUser(2*histos_data[i]->GetMinimum(),2*histos_mc[i]->GetMaximum());
      histos_mc[i]->Draw();
      histos_data[i]->Draw("same");
      auto rp = new TRatioPlot(histos_data[i], histos_mc[i], "divsym");
      c.SetTicks(0, 1);
      rp->SetH1DrawOpt("E");
      rp->Draw("nogrid");
      rp->GetLowerRefYaxis()->SetTitle("Data/MC");
      rp->GetUpperRefYaxis()->SetTitle("normalized entries");
      c.Update();
	
      TLegend* leg;

      leg = new TLegend(0.7, 0.7, 0.9, 0.9);
      leg->AddEntry(histos_data[i]->GetName(), "S. Subtraction", "l");
      leg->AddEntry(histos_mc[i]->GetName(), "Monte Carlo", "l");
      leg->SetTextSize(0.03);
      leg->Draw("same");

      c.SaveAs("mc_validation_plots/"+names[i]+"_mc_validation.pdf");
      leg->Delete();
      //tex->Delete();
      histos_mc[i]->Delete();
      histos_data[i]->Delete();

    }

  //delete ws;
  
  vector<TH1D*> histos_weight;

  //SPlot
  for(int i = 0;i<n_var;i++)   {
    //RooWorkspace* WS = new RooWorkspace("WS");
    set_up_workspace_variables(*ws);
    read_data(*ws, input_file_data, variables[i]);
    DoSPlot(*ws);
    histos_weight = GetSPlot(*ws, n_bins[i], variables[i]);
    //delete WS;
  }
}

void read_data(RooWorkspace& w, TString filename, TString variable)
{

  TFile* f = new TFile(filename);
  TNtupleD* _nt = (TNtupleD*)f->Get("ntKp");
  
  RooArgList arg_list ("arg_list");
  arg_list.add(*(w.var("Bmass")));
  arg_list.add(*(w.var(variable)));

  RooDataSet* data = new RooDataSet("data","data",_nt,arg_list);

  w.import(*data, Rename("data"));
}

void DoSPlot(RooWorkspace& ws)
{
 
  RooAbsPdf* model = ws.pdf("model");
  RooDataSet* data = (RooDataSet*)ws.data("data");
  
  RooRealVar* BpYield = ws.var("n_signal");
  RooRealVar* BgYield = ws.var("n_combinatorial");

  

  double sigYield = BpYield->getVal();
  double bkgYield = BgYield->getVal();

  cout<< "BpYield (before SPlot) = " << sigYield <<endl;
  cout<< "BgYield (before Splot) = " << bkgYield <<endl;
  
  model->fitTo(*data, Extended() );

  RooRealVar* mean = ws.var("mean");
  RooRealVar* sigma1 = ws.var("sigma1");
  RooRealVar* sigma2 = ws.var("sigma2");
  RooRealVar* lambda = ws.var("lambda");
  RooRealVar* m_jpsipi_mean1 = ws.var("m_jpsipi_mean1");
  RooRealVar* m_jpsipi_mean2 = ws.var("m_jpsipi_mean2");
  RooRealVar* m_jpsipi_mean3 = ws.var("m_jpsipi_mean3");
  RooRealVar* m_jpsipi_sigma1l = ws.var("m_jpsipi_sigma1l");
  RooRealVar* m_jpsipi_sigma1r = ws.var("m_jpsipi_sigma1r");
  RooRealVar* m_jpsipi_sigma2 = ws.var("m_jpsipi_sigma2");
  RooRealVar* m_jpsipi_sigma3 = ws.var("m_jpsipi_sigma3");
  
  mean ->setConstant();
  sigma1->setConstant();
  sigma2->setConstant();
  lambda->setConstant();
  m_jpsipi_mean1->setConstant();
  m_jpsipi_mean2->setConstant();
  m_jpsipi_mean3->setConstant();
  m_jpsipi_sigma1l->setConstant();
  m_jpsipi_sigma1r->setConstant();
  m_jpsipi_sigma2->setConstant();
  m_jpsipi_sigma3->setConstant();

  RooMsgService::instance().setSilentMode(true);

  SPlot* sData = new SPlot("sData","An sPlot",*data, model, RooArgList(*BpYield,*BgYield));

  cout << "BpYield (after SPlot) = " << sData->GetYieldFromSWeight("n_signal") << endl;
  cout << "BgYield (after SPlot) = " << sData->GetYieldFromSWeight("n_combinatorial") << endl;

  ws.import(*data, Rename("dataWithSWeights"));

}

vector<TH1D*> GetSPlot(RooWorkspace& ws, int nob, TString label){

  RooAbsPdf* model = ws.pdf("model");
  RooRealVar* BpModel = ws.var("signal");
  RooRealVar* BgModel = ws.var("fit_side");

  RooRealVar* Bmass = ws.var("Bmass");
  RooRealVar* variable = ws.var(label);

  RooRealVar* BpYield = ws.var("n_signal");
  RooRealVar* BgYield = ws.var("n_combinatorial");

  double sigYield = BpYield->getVal();
  double bkgYield = BgYield->getVal();

  

  //Get the data with weights
  RooDataSet* dataW = (RooDataSet*)ws.data("dataWithSWeights");
  RooDataSet* dataWBp = new RooDataSet(dataW->GetName(),dataW->GetTitle(),dataW,*dataW->get(),0,"n_signal_sw");
  RooDataSet* dataWBg = new RooDataSet(dataW->GetName(),dataW->GetTitle(),dataW,*dataW->get(),0,"n_combinatorial_sw");

  //Create histograms for data with weights
  TH1D* histo_Bp_sig = (TH1D*)dataWBp->createHistogram(variable,nob,0,0);
  TH1D* histo_Bp_bkg = (TH1D*)dataWBg->createHistogram(variable,nob,0,0);

  for (int i=1; i<=nob; i++) {

    if (histo_Bp_sig->GetBinContent(i)==0) histo_Bp_sig->SetBinError(i,0.);
    if (histo_Bp_bkg->GetBinContent(i)==0) histo_Bp_bkg->SetBinError(i,0.);

    cout << "HERE!!!!!!!!!!!!" << endl;

    histo_Bp_sig->SetBinContent(i,histo_Bp_sig->GetBinContent(i)/sigYield);
    histo_Bp_sig->SetBinError(i,histo_Bp_sig->GetBinError(i)/sigYield);

    histo_Bp_bkg->SetBinContent(i,histo_Bp_bkg->GetBinContent(i)/bkgYield);
    histo_Bp_bkg->SetBinError(i,histo_Bp_bkg->GetBinError(i)/bkgYield);
  }

  vector<TH1D*> histos;
  histos.push_back(histo_Bp_sig);
  histos.push_back(histo_Bp_bkg);

  return histos;
  
 
}

std::vector<TH1D*> sideband_subtraction(RooWorkspace* w,TString f_input, int* n){

  TFile* fin_data = new TFile(f_input);
//tal como há bocado lê o ficheiro root 
TTree* t1_data = (TTree*)fin_data->Get("ntKp");
//a partir do ficheiro devolve o histograma

RooRealVar Bmass = *(w->var("Bmass"));
RooRealVar Bpt = *(w->var("Bpt"));
RooRealVar By = *(w->var("By"));
RooRealVar Btrk1D0Err = *(w->var("Btrk1D0Err"));
RooRealVar Bmu1pt = *(w->var("Bmu1pt"));
RooRealVar Bmu1eta = *(w->var("Bmu1eta"));
RooRealVar Btrk1pt = *(w->var("Btrk1pt"));
RooRealVar Btrk1eta = *(w->var("Btrk1eta"));
RooRealVar Bchi2cl = *(w->var("Bchi2cl"));
RooRealVar BsvpvDistance = *(w->var("BsvpvDistance"));
RooRealVar BsvpvDistance_Err = *(w->var("BsvpvDistance_Err"));
RooRealVar Balpha = *(w->var("Balpha"));
RooRealVar Btrk1D0 = *(w->var("Btrk1D0"));
RooRealVar Btrk1Dz = *(w->var("Btrk1Dz"));
RooRealVar Bd0 = *(w->var("Bd0"));
RooRealVar Bd0err = *(w->var("Bd0err"));
RooRealVar Blxy = *(w->var("Blxy"));


RooArgSet arg_list(Bmass,Bpt,By, Btrk1D0Err, Bmu1pt, Bmu1eta, Btrk1pt, Btrk1eta);
arg_list.add(Bchi2cl);
arg_list.add(BsvpvDistance);
arg_list.add(Balpha);
arg_list.add(Btrk1D0);
arg_list.add(Btrk1Dz);
arg_list.add(Bd0);
arg_list.add(Blxy);
arg_list.add(Bd0err);
arg_list.add(BsvpvDistance_Err);
RooDataSet* data = new RooDataSet("data","data",t1_data,arg_list);

//COMEÇA O QUE INTERESSA//

 RooDataSet* reduceddata_side; //direita
 // RooDataSet* reduceddata_aux; 
 RooDataSet* reduceddata_central;

 double left = 5.15;
 double right = 5.4;
 double mass_peak = 5.265;


 // reduceddata_aux = (RooDataSet*) data->reduce(Form("Bmass<%lf", left));
 reduceddata_side = (RooDataSet*) data->reduce(Form("Bmass>%lf",right));

//reduceddata_side->append(*reduceddata_aux);

reduceddata_central = (RooDataSet*) data->reduce(Form("Bmass>%lf",left));
reduceddata_central = (RooDataSet*) reduceddata_central->reduce(Form("Bmass<%lf",right));

//SINAL//

 RooRealVar mean("mean","mean",mass_peak,5.26,5.29);

 RooRealVar sigma1("sigma1","sigma1",0.019,0.017,0.024);

 RooGaussian signal1("signal1","signal_gauss1",Bmass,mean,sigma1);

 RooRealVar sigma2("sigma2","sigma2",0.004,0.0035,0.010);

 RooGaussian signal2("signal2","signal_gauss2",Bmass,mean,sigma2);

 RooRealVar cofs("cofs", "cofs", 0.5, 0., 1.);

 RooAddPdf signal("signal", "signal", RooArgList(signal1,signal2),cofs);



 //BACKGROUND//

//ERROR FUNCTION//
 RooRealVar m_nonprompt_scale("m_nonprompt_scale", "m_nonprompt_scale", 4.74168e-02, 0, 1);
 //1.93204e-02, 0.001, 0.3);
 RooRealVar m_nonprompt_shift("m_nonprompt_shift", "m_nonprompt_shift", 5.14425, 4.5, 6.);
//5.14357e+00,5.12,5.16);
  
 m_nonprompt_shift.setConstant(kTRUE);
 m_nonprompt_scale.setConstant(kTRUE);

 RooGenericPdf erf("erf","erf","TMath::Erfc((Bmass-m_nonprompt_shift)/m_nonprompt_scale)", RooArgList(Bmass,m_nonprompt_scale,m_nonprompt_shift));
 
//EXPONENCIAL//
 RooRealVar lambda("lambda","lambda",-2.,-5.,0.0);
 RooExponential fit_side("fit_side", "fit_side_exp", Bmass, lambda);


 ////////////////////////////////////////////////
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

 ///////////////////////////////////////////////
  Bmass.setRange("all", Bmass.getMin(),Bmass.getMax());
  Bmass.setRange("right",right,Bmass.getMax());
  Bmass.setRange("left",Bmass.getMin(),left);
  Bmass.setRange("peak",left,right);

  std::cout<<"mass minimum: "<<Bmass.getMin()<<std::endl;
  std::cout<<"mass maximum: "<<Bmass.getMax()<<std::endl;


  //JUNTANDO TUDO//

  //valores dos N

  double n_signal_initial = data->sumEntries(TString::Format("abs(Bmass-%g)<0.05",mass_peak)) - data->sumEntries(TString::Format("abs(Bmass-%g)<0.10&&abs(Bmass-%g)>0.05",mass_peak,mass_peak));
  //não percebi bem porque se calcula o sinal inicial desta forma

  double n_combinatorial_initial = data->sumEntries() - n_signal_initial;

  RooRealVar n_signal("n_signal","n_signal",n_signal_initial,0.,data->sumEntries());

  RooRealVar n_combinatorial("n_combinatorial","n_combinatorial",n_combinatorial_initial,0.,data->sumEntries());

  RooRealVar f_erf("f_erf","f_erf",2.50259e-01,0,1);
  RooProduct n_erf("n_erf","n_erf",RooArgList(n_signal,f_erf));
  //o N para a error function já foi calculado e depende do número de sinal
  
  RooRealVar f_jpsipi("f_jpsipi","f_jpsipi",4.1E-5/1.026E-3,0.,0.1); //BF(jpsi_pi) = (4.1+-0.4)*10^-5 / BF(jpsi K) = (1.026+-0.031)*10^-3
  f_jpsipi.setConstant(kTRUE);
  RooProduct n_jpsipi("n_jpsipi","n_jpsipi",RooArgList(n_signal,f_jpsipi));
  //o N para o jpsipi já foi calculado e depende do número de sinal


  RooAddPdf model("model", "model", RooArgList(signal,fit_side,erf,jpsipi),RooArgList(n_signal,n_combinatorial,n_erf,n_jpsipi));

  model.fitTo(*data,Range("all"));
  //TUDO


  RooPlot* massframe = Bmass.frame(); //removi o título
  //o gráfico é feito em função da massa

  data->plotOn(massframe, RooFit::Name("Data"));
  //os dados que vão permanecer a preto

  model.plotOn(massframe, RooFit::Name("Fit"),Range("all"),LineColor(kRed),LineStyle(1),LineWidth(2));
 //o fit vai ser uma linha encarnada e contínua.

  model.plotOn(massframe, RooFit::Name("Combinatorial"),Components("fit_side"),Range("all"),LineColor(kBlue),LineStyle(kDashed));
  //o background combinatório que vai ficar a azul e tracejado

  model.plotOn(massframe, RooFit::Name("Signal"),Components("signal"),Range("all"),LineColor(kOrange),LineStyle(kDashed));
 //o sinal que vai ficar a laranja e ser uma linha a tracejado 

  model.plotOn(massframe, RooFit::Name("B->J/psi X"),Components("erf"),Range("all"),LineColor(kGreen+3),LineStyle(kDashed));
 //o background específico (partial reconstructed decays) que fica a verde e tracejado

  model.plotOn(massframe, RooFit::Name("B-> J/psi pi"),Components("jpsipi"),Range("all"),LineColor(kPink+10),LineStyle(kDashed));
  //parte correspondente ao Jpsipi que está a cor de rosa


 //parâmetros//

model.paramOn(massframe,Layout(0.60,0.90,0.75));
 //o tamanho da caixa de parâmetros vai de 60 a 90% do eixo dos xx, a parte de cima da caixa está a 75% do eixo dos yy.
massframe->getAttText()->SetTextSize(0.028);

  massframe->GetYaxis()->SetTitleOffset(1.3);
  massframe->SetXTitle("Bmass (GeV)");

  TCanvas d;
  TPad *p1 = new TPad("p1","p1",0.0,0.27,0.82,0.99);
  p1->SetBorderMode(1); 
  p1->SetFrameBorderMode(0); 
  p1->SetBorderSize(2);
  p1->SetBottomMargin(0.0);
  p1->Draw(); 
     
  TPad *p2 = new TPad("p2","p2",0.0,0.065,0.82,0.24);
  p2->SetTopMargin(0.);    
  p2->SetBorderMode(0);
  p2->SetBorderSize(2); 
  p2->SetFrameBorderMode(0); 
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

  double lambda_str = lambda.getVal();
  double lambda_err = lambda.getError();
//  double n_str = n.getVal();
 // double n_err = n.getError();
  double chis = massframe->chiSquare(); //chisquare

//  TLatex* tex12 = new TLatex(0.15, 0.25, Form("N = %.3lf #pm %.3lf",n_str,n_err));
  TLatex* tex12 = new TLatex(0.15, 0.85, Form("#lambda_{exp} = %.3lf #pm %.3lf",lambda_str,lambda_err));
  tex12->SetNDC(kTRUE);
  tex12->SetTextFont(42);
  tex12->SetTextSize(0.04);
  //tex12->Draw();
  TLatex* tex13 = new TLatex(0.15, 0.8, Form("#chi/DOF = %.3lf",chis));
  tex13->SetNDC(kTRUE);
  tex13->SetTextFont(42);
  tex13->SetTextSize(0.04);
  // tex13->Draw();


  //LEGENDA//

 
  TLegend *leg = new TLegend (0.4, 0.5, 0.6, 0.7);
  //ordem: x1,y1,x2,y2 (estão em percentagens)
  leg->SetTextSize(0.03);
  leg->AddEntry(massframe->findObject("Data"), "Data", "l");
  leg->AddEntry(massframe->findObject("B->J/psi X"), "B->J/psi X", "l");
  leg->AddEntry(massframe->findObject("Signal"), "Signal", "l");
  leg->AddEntry(massframe->findObject("Combinatorial"), "Combinatorial", "l");
  leg->AddEntry(massframe->findObject("B->J/psi pi"), "B->J/psi pi", "l");
  leg->AddEntry(massframe->findObject("Fit"),"Fit","l");
  leg->Draw("same");

  ///////////////////////////////////////////////////////////////////////
  //pull dists

  //Construct a histogram with the pulls  of the data w.r.t. the curve
  RooHist* pull_hist = massframe->pullHist("Data","Fit");
  // Create a new frame to draw the pull distribution and add the distribution to the frame
  RooPlot *pull_plot = Bmass.frame();

  pull_plot->addPlotable(static_cast<RooPlotable*>(pull_hist),"P");
  pull_plot->SetTitle("");
  pull_plot->GetXaxis()->SetTitle("");
  pull_plot->GetXaxis()->SetTitleFont(42);
  pull_plot->GetXaxis()->SetTitleSize(0.17);
  pull_plot->GetXaxis()->SetTitleOffset(1.09);

  pull_plot->GetXaxis()->SetLabelFont(40);
  pull_plot->GetXaxis()->SetLabelSize(0.17);
  pull_plot->GetXaxis()->SetLabelOffset(0.01);
  pull_plot->GetXaxis()->SetTickLength(0.15);
  
  pull_plot->GetYaxis()->SetTitle("Pull hist");
  pull_plot->GetYaxis()->SetTitleFont(42);  
  pull_plot->GetYaxis()->SetTitleSize(0.20);
  pull_plot->GetYaxis()->SetTitleOffset(0.15);

  pull_plot->GetYaxis()->SetLabelFont(40);
  pull_plot->GetYaxis()->SetLabelSize(0.14);
  pull_plot->GetYaxis()->SetLabelOffset(0.005);
  
  pull_plot->GetYaxis()->SetNdivisions(305);

  p2->cd();
  pull_plot->Draw();

  ///////////////////////////////////////////////////////////////////////

  w->import(model);
  //w->import(n_signal);
  //w->import(n_combinatorial);
 
  d.SaveAs("mc_validation_plots/fit_side.pdf"); //ACABA O QUE INTERESSA//
									  

  std::cout << std::endl << "chisquare: " << massframe->chiSquare() << std::endl;
  //  std::cout << "LogLikelihood: " << nll->getVal() << std::endl;

  //Integrating the background distribution

  //RooAbsReal* int_fit_side_left = fit_side.createIntegral(Bmass, "left");
  RooAbsReal* int_fit_side_right = fit_side.createIntegral(Bmass, "right");
  RooAbsReal* int_fit_peak = fit_side.createIntegral(Bmass, "peak");

  //std::cout<< std::endl << "Integral left band: " << int_fit_side_left->getVal() << std::endl;
  std::cout<< std::endl << "Integral right band: " << int_fit_side_right->getVal() << std::endl;

//  double factor = (int_fit_peak->getVal())/(int_fit_side_left->getVal()+int_fit_side_right->getVal());
  double factor = (int_fit_peak->getVal())/(int_fit_side_right->getVal());

  std::cout << std::endl << "Factor: " << factor << std::endl;
  for(int i=0; i<16; i++){
   std::cout << "bins: " << n[i] << std::endl;
  } 

  std::vector<TH1D*> histos;

  histos.push_back(create_histogram(Bpt,"Bpt", factor, reduceddata_side, reduceddata_central, data, n[0]));
  histos.push_back(create_histogram(By, "By",factor, reduceddata_side, reduceddata_central, data, n[1]));
  histos.push_back(create_histogram(Btrk1D0Err, "Btrk1D0Err",factor, reduceddata_side, reduceddata_central, data, n[2]));
  histos.push_back(create_histogram(Bmu1pt, "Bmu1pt",factor, reduceddata_side, reduceddata_central, data, n[3]));
  histos.push_back(create_histogram(Bmu1eta, "Bmu1eta",factor, reduceddata_side, reduceddata_central, data, n[4]));
  histos.push_back(create_histogram(Btrk1pt, "Btrk1pt",factor, reduceddata_side, reduceddata_central, data, n[5]));
  histos.push_back(create_histogram(Btrk1eta, "Btrk1eta",factor, reduceddata_side, reduceddata_central, data, n[6]));
  histos.push_back(create_histogram(Bchi2cl, "Bchi2cl",factor, reduceddata_side, reduceddata_central, data, n[7]));
  histos.push_back(create_histogram(BsvpvDistance, "BsvpvDistance",factor, reduceddata_side, reduceddata_central, data, n[8]));
  histos.push_back(create_histogram(BsvpvDistance_Err, "BsvpvDistance_Err",factor, reduceddata_side, reduceddata_central, data, n[9]));
  histos.push_back(create_histogram(Balpha, "Balpha",factor, reduceddata_side, reduceddata_central, data, n[10]));
  histos.push_back(create_histogram(Btrk1D0, "Btrk1D0",factor, reduceddata_side, reduceddata_central, data, n[11]));
  histos.push_back(create_histogram(Btrk1Dz, "Btrk1Dz",factor, reduceddata_side, reduceddata_central, data, n[12]));
  histos.push_back(create_histogram(Bd0, "Bd0",factor, reduceddata_side, reduceddata_central, data, n[13]));
  histos.push_back(create_histogram(Blxy, "Blxy",factor, reduceddata_side, reduceddata_central, data, n[14]));
  histos.push_back(create_histogram(Bd0err, "Bd0err",factor, reduceddata_side, reduceddata_central, data, n[15]));
return histos;
//o que a sideband_subtraction me devolve para eu usar na função main para obter os histogramas dos dados
}


TH1D* create_histogram(RooRealVar var, TTree* t, int n){

  TH1D* h = new TH1D(var.GetName(), var.GetName(), n, var.getMin(), var.getMax());

  TString name_string = TString(var.GetName()) + ">>htemp(" + Form("%d",n) +"," + Form("%lf", var.getMin()) + "," + Form("%lf", var.getMax()) + ")";

  t->Draw(name_string, "Pthatweight");

  h = (TH1D*)gDirectory->Get("htemp")->Clone();
  h->SetTitle("");
  h->SetMarkerColor(kBlack);
  h->SetLineColor(kBlack);
  return h;
}


//def da função
TH1D* create_histogram(RooRealVar var,TString name, double factor, RooDataSet* reduced, RooDataSet* central, RooDataSet* total, int n){


  std::cout<< "n in create_histogram = "<< n <<std::endl;
  TH1D* dist_side = (TH1D*)reduced->createHistogram("dist_side",var, Binning(n, var.getMin(), var.getMax()));
  dist_side->SetMarkerColor(kBlue);
  dist_side->SetLineColor(kBlue);
  dist_side->SetNameTitle("dist_side", "");

  TH1D* hist_dist_peak = (TH1D*)central->createHistogram(var.GetName(), var, Binning(n, var.getMin(), var.getMax()));
  TH1D* dist_peak = new TH1D(*hist_dist_peak);
  dist_peak->SetMarkerColor(kRed);
  dist_peak->SetLineColor(kRed);
  dist_peak->SetNameTitle(var.GetName(), "");

  hist_dist_peak->SetMarkerColor(kBlack);
  hist_dist_peak->SetLineColor(kBlack);
  hist_dist_peak->SetNameTitle("dist_total", "");

  dist_peak->Add(dist_side, -factor);
  //dist_side->Add(dist_side, factor);
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

  //pt_dist_total->GetYaxis()->SetRangeUser(-1500.,70000.);
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
  leg->AddEntry("dist_total", "Total", "l");
  leg->AddEntry(var.GetName(), "Signal", "l");
  leg->AddEntry("dist_side", "Background", "l");
  leg->AddEntry("hist_dist_peak", "Sinal inicial", "l");
  leg->Draw("same");

  std::cout<<"name: "<<var.GetName()<<std::endl;
  std::cout<<"histo name: "<<dist_peak->GetName()<<std::endl;

  c.SaveAs("mc_validation_plots/sideband_sub/"+name + "sideband_sub.pdf");

  return dist_peak;

}

void set_up_workspace_variables(RooWorkspace& w)
{
  double mass_min, mass_max;
  double pt_min, pt_max;
  double y_min, y_max;
  double trkd0err_min, trkd0err_max;
  double mu1pt_min, mu1pt_max;
  double mu1eta_min, mu1eta_max;
  double trk1pt_min, trk1pt_max;
  double trk1eta_min, trk1eta_max;
  double chi2cl_min, chi2cl_max;
  double svpvDistance_min, svpvDistance_max;
  double alpha_min, alpha_max;
  double trk1D0_min, trk1D0_max;
  double trk1Dz_min, trk1Dz_max;
  double d0_min, d0_max;
  double lxy_min, lxy_max;
  double d0Err_min, d0Err_max;
  double svpvDistanceErr_min, svpvDistanceErr_max;

  mass_min=5.;
  mass_max=6.;

  pt_min=5.;
  pt_max=100.;

  y_min=-2.4;
  y_max=2.4;

  trkd0err_min = 0.;
  trkd0err_max = 0.01;

  mu1pt_min=0.;
  mu1pt_max=20.;


  mu1eta_min=-2.;
  mu1eta_max=2.;


  trk1pt_min=0.;
  trk1pt_max=8.;


  trk1eta_min=-3.;
  trk1eta_max=3.;

  chi2cl_min = 0.;
  chi2cl_max = 1.;

  svpvDistance_min=0.;
  svpvDistance_max=2.;
  alpha_min=0.;
  alpha_max=0.1;
  trk1D0_min=-0.3;
  trk1D0_max=0.3;
  trk1Dz_min=-0.1;
  trk1Dz_max=0.1;
  d0_min=0.; 
  d0_max=0.4;
  lxy_min=0.;
  lxy_max=3.;
  d0Err_min=0.;
  d0Err_max=0.0001;
  svpvDistanceErr_min=0.;
  svpvDistanceErr_max=0.05;

  //Depois de dar valores às variáveis de cima vou usá-los para o intervalo das variáveis que se seguem
  //RooRealVar=estou a construir/definir variáveis no roofit
  RooRealVar Bmass("Bmass","Bmass",mass_min,mass_max);
  RooRealVar Bpt("Bpt","Bpt",pt_min,pt_max);
  RooRealVar By("By","By",y_min,y_max);
  RooRealVar Btrk1D0Err("Btrk1D0Err","Btrk1D0Err",trkd0err_min,trkd0err_max);
  RooRealVar Bmu1pt("Bmu1pt","Bmu1pt",mu1pt_min,mu1pt_max);
  RooRealVar Bmu1eta("Bmu1eta","Bmu1eta",mu1eta_min,mu1eta_max);
  RooRealVar Btrk1pt("Btrk1pt","Btrk1pt",trk1pt_min,trk1pt_max);
  RooRealVar Btrk1eta("Btrk1eta","Btrk1eta",trk1eta_min,trk1eta_max);
  RooRealVar Bchi2cl("Bchi2cl", "Bchi2cl", chi2cl_min, chi2cl_max);
  RooRealVar BsvpvDistance("BsvpvDistance", "BsvpvDistance", svpvDistance_min, svpvDistance_max);
  RooRealVar Balpha("Balpha", "Balpha", alpha_min, alpha_max);
  RooRealVar Btrk1D0("Btrk1D0","Btrk1D0",trk1D0_min,trk1D0_max);
  RooRealVar Btrk1Dz("Btrk1Dz","Btrk1Dz",trk1Dz_min,trk1Dz_max);
  RooRealVar Bd0("Bd0", "Bd0", d0_min, d0_max); 
  RooRealVar Blxy("Blxy", "Blxy", lxy_min, lxy_max);
  RooRealVar Bd0err("Bd0err", "Bd0err", d0Err_min, d0Err_max);
  RooRealVar BsvpvDistance_Err("BsvpvDistance_Err", "BsvpvDistance_Err", svpvDistanceErr_min, svpvDistanceErr_max);

  w.import(Bmass);
  w.import(Bpt);
  w.import(By);
  w.import(Btrk1D0Err);
  w.import(Bmu1pt);
  w.import(Btrk1pt);
  w.import(Bmu1eta);
  w.import(Btrk1eta);
  w.import(Bchi2cl);
  w.import(BsvpvDistance);
  w.import(Balpha);
  w.import(Btrk1D0);
  w.import(Btrk1Dz);
  w.import(Bd0);
  w.import(Blxy);
  w.import(Bd0err);
  w.import(BsvpvDistance_Err);
}
