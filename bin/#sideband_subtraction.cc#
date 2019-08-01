

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
#include "TMath.h"
#include <RooGenericPdf.h>
#include "TRatioPlot.h"


using namespace RooFit;
using namespace std;
//using namespace RooStats;

void set_up_workspace_variables(RooWorkspace& w);
TH1D* create_histogram(RooRealVar var, TTree* t, int n);
TH1D* create_histogram(RooRealVar var,TString name, double factor, RooDataSet* reduced, RooDataSet* central, RooDataSet* total, int n);
std::vector<TH1D*> sideband_subtraction(RooWorkspace* w,TString f_input, int* n);


int main(){

  TString input_file_data = "/home/t3cms/ev19u032/test/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/prefiltered_trees/selected_data_ntKp_PbPb_2018.root"; 
//input - dados

TString input_file_mc = "/home/t3cms/ev19u032/test/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/prefiltered_trees/selected_mc_ntKp_PbPb_2018_pthatweight.root";
//input - monte carlo

std::vector<TH1D*> histos_data;
//vetor por enquanto vazio que vai ter os histogramas dos dados
std::vector<TH1D*> histos_mc;
//vetor por enquanto vazio que vai ter os histogramas de MC
int n_bins[]= {10, 20, 10, 10, 10, 10, 10, 10, 15, 10, 10, 10, 10, 10, 10, 10, 10, 10};

RooWorkspace* ws = new RooWorkspace("ws","Bmass");
//The RooWorkspace is a persistable container for RooFit projects.A workspace can contain and own variables, p.d.f.s, functions and datasets. All objects that live in the workspace are owned by the workspace.

set_up_workspace_variables(*ws);
//vou guardar variáveis no workspace (ver abaixo junto à definição de sideband_subtraction)

//HISTOGRAMA DOS DADOS
histos_data = sideband_subtraction(ws, input_file_data, n_bins);
//usar sideband_subtraction para obter o histograma dos dados(ver abaixo como esta está definida junto às variáveis no workspace)
//PARA O HISTROGAMA MC
TFile *fin_mc = new TFile(input_file_mc);
//lê o ficheiro root
TTree* t1_mc = (TTree*)fin_mc->Get("ntKp");
//a partir do ficheiro devolve o histograma

std::vector<TString> names;
//crio o vetor que vai conter os nomes

//percorro o vetor dos dados e dou os nomes às variáveis
for(int i=0; i<(int)histos_data.size(); i++){
std::cout<< "Var names: "<< histos_data[i]->GetName()<<std::endl;
}
//percorrro o vetor dos dados e adiciono ao vetor de monte carlo os nomes das variáveis(.push_back=preencho o vetor)
for(int i=0; i<(int)histos_data.size(); i++){histos_mc.push_back(create_histogram((*ws->var(histos_data[i]->GetName())), t1_mc, n_bins[i]));
names.push_back(TString(histos_data[i]->GetName()));

//com isto fui buscar as variáveis às trees e coloquei-as no roofit

 } 


      for(int i=0; i<(int)histos_data.size(); i++)
	{
	  TCanvas c;
	      histos_mc[i]->SetXTitle(TString(histos_data[i]->GetName()));
              histos_mc[i]->SetStats(0);
              histos_data[i]->SetStats(0);
	      histos_mc[i]->Scale(1/histos_mc[i]->Integral());
	      // histos_mc[i]->GetYaxis()->SetRangeUser(-0.1,2*histos_mc[i]->GetMaximum());
	      // histos_mc[i]->GetYaxis()->SetRangeUser(histos_data[i]->GetMinimum(),2*histos_mc[i]->GetMaximum());
	      // histos_mc[i]->GetYaxis()->SetRangeUser(-0.1,0.5*histos_data[i]->GetMaximum());
	      histos_data[i]->Scale(1/histos_data[i]->Integral());
	      // histos_mc[i]->GetYaxis()->SetRangeUser(histos_data[i]->GetMinimum(),2*histos_mc[i]->GetMaximum());
	      histos_mc[i]->GetYaxis()->SetRangeUser(2*histos_data[i]->GetMinimum(),2*histos_mc[i]->GetMaximum());
	      histos_mc[i]->Draw();
	      histos_data[i]->Draw("same");
	      //Fazendo o Tratio para ter o racio entre os dados e mc

	      auto rp = new TRatioPlot(histos_mc[i], histos_data[i],"divsym");
              c.SetTicks(0, 1);
	      rp->SetH1DrawOpt("E");
              rp->Draw("nogrid");
	      rp->GetLowerRefYaxis()->SetTitle("Data/MC");
	      rp->GetUpperRefYaxis()->SetTitle("normalized entries");
	
	  TLegend* leg;

	  leg = new TLegend(0.7, 0.7, 0.9, 0.9);

	  leg->AddEntry(histos_data[i]->GetName(), "S. Subtraction", "l");
	  leg->AddEntry(histos_mc[i]->GetName(), "Monte Carlo", "l");
	  leg->SetTextSize(0.03);
	  //  std::cout<<"NOME DO DATA: "<< histos_data[i]->GetName()<< std::endl;
	  // std::cout<<"NOME DO MC: "<< histos_mc[i]->GetName()<< std::endl;
	  leg->Draw("same");
	  //"DESENHOU"

    /*TLatex* tex = new TLatex(0.6,0.8,"1.5 nb^{-1} (PbPb) 5.02 TeV");
    tex->SetNDC(kTRUE);
    tex->SetLineWidth(2);
    tex->SetTextSize(0.04);
    tex->Draw();
    tex = new TLatex(0.68,0.85,"CMS Preliminary");
    tex->SetNDC(kTRUE);
    tex->SetTextFont(42);
    tex->SetTextSize(0.04);
    tex->SetLineWidth(2);
    tex->Draw();*/


	  //SALVA
    c.SaveAs("mc_validation_plots/"+names[i]+"_mc_validation.pdf");
    leg->Delete();
    //tex->Delete();
    histos_mc[i]->Delete();
    histos_data[i]->Delete();

	}//AQUI TERMINA O CICLO FOR


}//AQUI TERMINA A MAIN FUNCTION
//defino sideband_subtraction que vai ser usada para o histograma dos dados

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
 double mass_peak = 5.26;


 // reduceddata_aux = (RooDataSet*) data->reduce(Form("Bmass<%lf", left));
 reduceddata_side = (RooDataSet*) data->reduce(Form("Bmass>%lf",right));

//reduceddata_side->append(*reduceddata_aux);

reduceddata_central = (RooDataSet*) data->reduce(Form("Bmass>%lf",left));
reduceddata_central = (RooDataSet*) reduceddata_central->reduce(Form("Bmass<%lf",right));

//SINAL//

RooRealVar mean("mean","mean",mass_peak,5.15,5.4);

 RooRealVar sigma1("sigma1","sigma1",0.020,0.018,0.024);

RooGaussian signal1("signal1","signal_gauss1",Bmass,mean,sigma1);

 RooRealVar sigma2("sigma2","sigma2",0.015,0.010,0.018);

RooGaussian signal2("signal2","signal_gauss2",Bmass,mean,sigma2);

 RooRealVar sigma3("sigma3","sigma3",0.025,0.020,0.027);

RooGaussian signal3("signal3","signal_gauss3",Bmass,mean,sigma3);

RooRealVar cofs("cofs", "cofs", 0.5, 0., 1.);

 RooAddPdf signal("signal", "signal", RooArgList(signal1,signal2),RooArgList(cofs,cofs));



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

 Bmass.setRange("all", Bmass.getMin(),Bmass.getMax());
 Bmass.setRange("right",right,Bmass.getMax());
 Bmass.setRange("left",Bmass.getMin(),left);
 Bmass.setRange("peak",left,right);

 std::cout<<"mass minimum: "<<Bmass.getMin()<<std::endl;
 std::cout<<"mass maximum: "<<Bmass.getMax()<<std::endl;




  //JUNTANDO TUDO//

  RooRealVar cof1("cof", "cof", 0.5, 0., 1.);
  RooRealVar cof2("cof", "cof", 0.5, 0., 1.);

  RooAddPdf model("model", "model", RooArgList(signal,fit_side,erf),RooArgList(cof1,cof2));



  fit_side.fitTo(*reduceddata_side,Range("right"));
  //SIDEBAND RIGHT

  //erf.fitTo(*reduceddata_aux, Range("left"));
  erf.fitTo(*data, Range("left"));

  model.fitTo(*data,Range("all"));
  //TUDO


  RooPlot* massframe = Bmass.frame(Title("")); //removi o título
  //o gráfico é feito em função da massa

  data->plotOn(massframe, RooFit::Name("data"));
  //os dados que vão permanecer a preto

  // signal.plotOn(massframe,RooFit::Name("signal"),Range(5.15,5.4),LineColor(kOrange),LineStyle(kDashed));
 //o sinal que vai ficar a laranja e ser uma linha a tracejado 

 fit_side.plotOn(massframe, RooFit::Name("fit_side"),Range("all"),LineColor(kBlue),LineStyle(kDashed));
  //o background combinatório que vai ficar a azul e tracejado

 // erf.plotOn(massframe, RooFit::Name("erf"),Range("left"),LineColor(kGreen),LineStyle(kDashed));
 //mal feito --> corrigir
 //o background específico (partial reconstructed decays) que fica a verde e tracejado


 model.plotOn(massframe, RooFit::Name("model"),Range("all"),LineColor(kRed),LineStyle(kDashed));
 //o fit vai ser uma linha encarnada e tracejada

  massframe->GetYaxis()->SetTitleOffset(1.3);
  massframe->SetXTitle("Bmass (GeV)");
//  massframe->SetNameTitle("sideband_fit", "Exponential Fit - Sideband Subtraction");

  TCanvas d;
  massframe->Draw();
  TLatex* tex11 = new TLatex(0.6,0.8,"1.5 nb^{-1} (PbPb) 5.02 TeV");
  tex11->SetNDC(kTRUE);
  tex11->SetLineWidth(2);
  tex11->SetTextSize(0.04);
  tex11->Draw();
  tex11 = new TLatex(0.68,0.8,"CMS Preliminary");
  tex11->SetNDC(kTRUE);
  tex11->SetTextFont(42);
  tex11->SetTextSize(0.04);
  tex11->SetLineWidth(2);
  //tex11->Draw();

  double lambda_str = lambda.getVal();
  double lambda_err = lambda.getError();
//  double n_str = n.getVal();
 // double n_err = n.getError();
  double chis = massframe->chiSquare();

//  TLatex* tex12 = new TLatex(0.15, 0.25, Form("N = %.3lf #pm %.3lf",n_str,n_err));
  TLatex* tex12 = new TLatex(0.15, 0.85, Form("#lambda_{exp} = %.3lf #pm %.3lf",lambda_str,lambda_err));
  tex12->SetNDC(kTRUE);
  tex12->SetTextFont(42);
  tex12->SetTextSize(0.04);
  tex12->Draw();
  TLatex* tex13 = new TLatex(0.15, 0.8, Form("#chi/DOF = %.3lf",chis));
  tex13->SetNDC(kTRUE);
  tex13->SetTextFont(42);
  tex13->SetTextSize(0.04);
  //tex13->Draw();


  //LEGENDA//

 
  TLegend *leg = new TLegend (0.7, 0.5, 0.85, 0.65);
  leg->AddEntry(massframe->findObject("data"), "Data", "l");
  leg->AddEntry(massframe->findObject("erf"), "B->J/psi X", "l");
  leg->AddEntry(massframe->findObject("signal"), "Signal", "l");
  leg->AddEntry(massframe->findObject("fit_side"), "Background fit", "l");
  leg->AddEntry(massframe->findObject("model"),"Global Fit","l");
  leg->Draw("same");
 
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
