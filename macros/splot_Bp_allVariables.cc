/* example code illustrating the likelihood-based procedure for determining the angular distribution in data for signal,
without sideband-like subtraction the roostats splot implementaiton is utilized
 - code uses Bmass as (signal vs background) discriminating variable
 - polarization sample is generating assigning distinct polarization values for signal and background (lambda theta = +/-1)
 - a likelihood fit is performed to the data, using the Bmass dimension
 - the obtained per-event likelihood is used to weight the dataset to extract the signal angular distribution
 * comments on macro to nuno.leonardo@cern.ch 12/5/2013
 */

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


using namespace RooFit;
using namespace RooStats;
using namespace std;

void DoSPlot(RooWorkspace&,int);
void MakePlots(RooWorkspace&,int,TString);
void build_pdf(RooWorkspace&,int);
void read_data(RooWorkspace&,TString,TString,int);
TString channel_to_ntuple_name(int);
void set_up_workspace_variables(RooWorkspace&);


void splot_Bp_allVariables() {

const int n_var = 16;
int channel = 1; // B+
int n_bins[n_var]= {10, 20, 10, 10, 10, 10, 10, 10, 15, 10, 10, 10, 10, 10, 10, 10};
TString variables [n_var] = {"Bpt","By","Btrk1D0Err","Bmu1pt","Bmu1eta","Btrk1pt","Btrk1eta","Bchi2cl","BsvpvDistance","BsvpvDistance_Err","Balpha","Btrk1D0","Btrk1Dz","Bd0","Blxy","Bd0err"};
TString input_file_Bp = "/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/prefiltered_trees/selected_data_ntKp_PbPb_2018.root";

    for(int i = 0;i<n_var;i++)   { // loop over all the variables

      cout<<"Considering "<<variables[i]<<" variable...\n";
  RooWorkspace* ws = new RooWorkspace("WS");
  set_up_workspace_variables(*ws);
  read_data(*ws, input_file_Bp, variables[i], channel);


  cout<<"\nWorkSpace CONTENTS\n";
  ws->Print();   // print workspace contents
  DoSPlot(*ws,channel);   // fit data with 1D model, add sWeights
  MakePlots(*ws,n_bins[i],variables[i]); // project data for signal and background
  delete ws;     // cleanup
  }

  cout<<"\n...loop over variables done.\n";
}


//____________________________________
void DoSPlot(RooWorkspace& ws,int channel){

  cout << "\ndoing sPlot...\n" << endl;

  build_pdf(ws,channel);
  RooAbsPdf*  model = ws.pdf("model");
  RooDataSet* data = (RooDataSet*) ws.data("data");

  RooRealVar* BpYield = ws.var("n_signal");
  RooRealVar* BgYield = ws.var("n_combinatorial");

  cout<<"BpYield = "<<BpYield->getVal()<<endl;
  cout<<"BgYield = "<<BgYield->getVal()<<endl;

  // fit the model to the data
  model->fitTo(*data,Extended());

  // sPlot technique requires model parameters (other than the yields) to be fixed
  RooRealVar* mean  = ws.var("m_mean");
  RooRealVar* sigma = ws.var("m_sigma1");
  RooRealVar* dec   = ws.var("m_exp");

  mean ->setConstant();
  sigma->setConstant();
  dec  ->setConstant();


  RooMsgService::instance().setSilentMode(true);

  // add sWeights to dataset based on model and yield variables
  // sPlot class adds a new variable that has the name of the corresponding yield + "_sw".

  SPlot* sData = new SPlot("sData","An sPlot",*data, model, RooArgList(*BpYield,*BgYield));

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

  // import this new dataset with sWeights
  ws.import(*data, Rename("dataWithSWeights"));

  cout << "\n...done doSplot.\n" << endl;
} // "DO SPLOT" ENDS


void MakePlots(RooWorkspace& ws, int nob,TString label){

  cout << "\nplotting...\n" << endl;

  // Here we make plots of the discriminating variable (Bmass) after the fit
  // and of the control variable (angle) after unfolding with sPlot.

  TCanvas* cdata = new TCanvas("sPlot","sPlot", 800, 600);
  cdata->Divide(2,2);

  RooAbsPdf* model  = ws.pdf("model");
  RooAbsPdf* BpModel = ws.pdf("pdf_m_signal");
  RooAbsPdf* BgModel = ws.pdf("pdf_m_combinatorial");

  RooRealVar* Bmass  = ws.var("Bmass");
  RooRealVar* variable = ws.var(label);

  RooRealVar* BpYield = ws.var("n_signal");
  RooRealVar* BgYield = ws.var("n_combinatorial");


  double sigYield = BpYield->getVal();
  double bkgYield = BgYield->getVal();

  cout<<"sigYield = "<<sigYield<<endl;
  cout<<"bkgYield = "<<bkgYield<<endl;

  RooDataSet* data = (RooDataSet*) ws.data("data");

  //this shouldn't be necessary (set pars to their fit values)
  //model->fitTo(*data, Extended());

  cdata->cd(1);
  RooPlot* mframe = Bmass->frame();
  mframe->GetXaxis()->SetTitle(TString::Format("mass of B+ [GeV]"));
  data->plotOn(mframe);
  model->plotOn(mframe);
  model->plotOn(mframe,Components(*BpModel),LineStyle(kDashed),LineColor(kRed));
  model->plotOn(mframe,Components(*BgModel),LineStyle(kDashed),LineColor(kGreen));
  mframe->SetTitle("Bmass");
  mframe->Draw();

  cdata->cd(2);
  RooPlot* ptframe = variable->frame();
  data->plotOn(ptframe);
  ptframe->SetTitle(label+" of B+: total sample");
  ptframe->Draw();

  //get the dataset with sWeights
  RooDataSet* dataW = (RooDataSet*) ws.data("dataWithSWeights");
  RooDataSet* dataWBp = new RooDataSet(dataW->GetName(),dataW->GetTitle(),dataW,*dataW->get(),0,"n_signal_sw");
  RooDataSet* dataWBg = new RooDataSet(dataW->GetName(),dataW->GetTitle(),dataW,*dataW->get(),0,"n_combinatorial_sw");

  RooPlot* ptframe2Bp = variable->frame();
  RooPlot* ptframe2Bg = variable->frame();

  ptframe2Bp->GetYaxis()->SetTitle(TString::Format("Events / (%g)",(variable->getMax()-variable->getMin())/nob));
  ptframe2Bg->GetYaxis()->SetTitle(TString::Format("Events / (%g)",(variable->getMax()-variable->getMin())/nob));

  ptframe2Bp->GetXaxis()->SetTitle(label + " of B+");
  ptframe2Bg->GetXaxis()->SetTitle(label + " of B+");

  dataWBp->plotOn(ptframe2Bp, DataError(RooAbsData::SumW2),Binning(nob));
  dataWBg->plotOn(ptframe2Bg, DataError(RooAbsData::SumW2),Binning(nob));

  ptframe2Bp->SetTitle(label+" distribution of B+ for signal (splot)");
  ptframe2Bg->SetTitle(label+" distribution of B+ for background (splot)");

  cdata->cd(3);  ptframe2Bp->Draw();
  cdata->cd(4);  ptframe2Bg->Draw();

  cdata->SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/splot/sPlot_Bplus_Bmass_"+label+".gif");
  cdata->SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/splot/sPlot_Bplus_Bmass_"+label+".pdf");


  TH1D* histo_Bp_sig = (TH1D*)dataWBp->createHistogram(label,nob,0,0);
  TH1D* histo_Bp_bkg = (TH1D*)dataWBg->createHistogram(label,nob,0,0);

    for (int i=1; i<=nob; i++) {

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
  histo_Bp_sig->SetMarkerColor(kBlue);
  histo_Bp_sig->SetLineColor(kBlue);
  histo_Bp_sig->SetTitle("");
  histo_Bp_sig->GetYaxis()->SetTitle(TString::Format("Events /(%g)",(variable->getMax()-variable->getMin())/nob));
  histo_Bp_sig->GetXaxis()->SetTitle(label );

  histo_Bp_sig->SetStats(0);
  histo_Bp_sig->Draw("E");

  prov->SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/splot/sPlot_Bplus_"+label+".gif");
  prov->SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/splot/sPlot_Bplus_"+label+".pdl");

  TCanvas* prov_bkg = new TCanvas ("prov_bkg","c2",200,10,700,500);
  prov_bkg->cd();
  histo_Bp_bkg->SetMarkerStyle(20);
  histo_Bp_bkg->SetMarkerSize(0.);
  histo_Bp_bkg->SetMarkerColor(kRed);
  histo_Bp_bkg->SetLineColor(kRed);
  histo_Bp_bkg->SetTitle("");
  histo_Bp_bkg->GetYaxis()->SetTitle(TString::Format("Events /(%g)",(variable->getMax()-variable->getMin())/nob));
  histo_Bp_bkg->GetXaxis()->SetTitle(label);

  histo_Bp_bkg->SetStats(0);
  histo_Bp_bkg->Draw("E");

  prov_bkg->SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/splot/sPlot_Bplus_"+label+".gif");
  prov_bkg->SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/splot/sPlot_Bplus_"+label+".pdf");


  TCanvas* sig_bkg = new TCanvas ("sig_bkg","c3",200,10,700,500);
  sig_bkg->cd();

  histo_Bp_sig->Draw();
  histo_Bp_bkg->Draw("SAME");

  TLegend* legend = new TLegend(0.7,0.9,0.9,0.8);
   legend->AddEntry(histo_Bp_sig,"Signal","lep");
   legend->AddEntry(histo_Bp_bkg,"Background","lep");
   legend->Draw();

  sig_bkg->SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/splot/sPlot_SigBkg_"+label+".gif");
  sig_bkg->SaveAs("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/splot/sPlot_SigBkg_"+label+".pdf");


  cout << "\n...done plotting.\n" << endl;

  //cleanup
  delete cdata;
  delete prov;
  delete prov_bkg;
  delete sig_bkg;
} // "MAKE SPLOT" ENDS

TString channel_to_ntuple_name(int channel)
{
  //returns a TString with the ntuple name corresponding to the channel. It can be used to find the data on each channel saved in a file. or to write the name of a directory
  TString ntuple_name = "";

  switch(channel){
  default:
  case 1:
    ntuple_name="ntKp";
    break;
  case 2:
    ntuple_name="ntKstar";
    break;
  case 3:
   ntuple_name="ntKs";
   break;
  case 4:
    ntuple_name="ntphi";
    break;
  case 5:
   ntuple_name="ntmix";
   break;
  case 6:
    ntuple_name="ntlambda";
    break;
 }
  return ntuple_name;
} // "CHANNEL TO NTUPLE NAME" ENDS


void read_data(RooWorkspace& w, TString filename,TString var_label ,int channel) //SELECT DATA
{

  TFile* f = new TFile(filename);
  std::cout<<"new file"<<std::endl;
  TNtupleD* _nt = (TNtupleD*)f->Get(channel_to_ntuple_name(channel));
  std::cout<<"new ntuple"<<std::endl;

  RooArgList arg_list ("arg_list");
  arg_list.add(*(w.var("Bmass")));
  arg_list.add(*(w.var(var_label)));

  //RooRealVar* Bmass = w.var("Bmass");
  //RooRealVar* variable = w.var(var_label);

  RooDataSet* data = new RooDataSet("data","data",_nt,arg_list);

  w.import(*data, Rename("data"));

} // "READ DATA" ENDS

void build_pdf(RooWorkspace& w,int channel) {

 double mass_peak;
 RooRealVar Bmass = *(w.var("Bmass"));
 RooRealVar Bpt = *(w.var("Bpt"));
 RooAbsData* data = w.data("data");

switch (channel) {
 default:
 case 1:
   mass_peak = 5.27926; //BP_MASS
   break;
case 2:
  mass_peak = 5.27958;  //B0_MASS
  break;
 case 4:
   mass_peak = 5.36677; //BS_MASS
   break;
 }

 double n_signal_initial = data->sumEntries(TString::Format("abs(Bmass-%g)<0.015",mass_peak)) - data->sumEntries(TString::Format("abs(Bmass-%g)<0.030&&abs(Bmass-%g)>0.015",mass_peak,mass_peak));
 if(n_signal_initial<0)  n_signal_initial=1;

 double n_combinatorial_initial = data->sumEntries() - n_signal_initial;

 cout<<"Signal initial = "<<n_signal_initial<<endl<<"Bg initial = "<<n_combinatorial_initial<<endl<<"Total entries = "<<data->sumEntries()<<endl;

  RooRealVar m_mean("m_mean","m_mean",mass_peak,mass_peak-0.09,mass_peak+0.09);
  RooRealVar m_sigma1("m_sigma1","m_sigma1",0.010,0.009,0.200);
  RooRealVar m_sigma2("m_sigma2","m_sigma2",0.005,0.004,0.100);
  RooRealVar m_fraction("m_fraction","m_fraction", 0.5, 0, 1);
  RooGaussian m_gaussian1("m_gaussian1","m_gaussian1",Bmass,m_mean,m_sigma1);
  RooGaussian m_gaussian2("m_gaussian2","m_gaussian2",Bmass,m_mean,m_sigma2);

  RooAddPdf* pdf_m_signal;
  RooAddPdf* pdf_m_combinatorial;

  // use single Gaussian for low statistics
  if(n_signal_initial < 500)
  {
    m_sigma2.setConstant(kTRUE);
    m_fraction.setVal(1.);
    m_fraction.setConstant(kTRUE);
    std::endl(std::cout);
    std::cout<<"The initial signal was indeed < 500"<<std::endl;
    std::endl(std::cout);

  }

  //One Exponential
   RooRealVar m_exp("m_exp","m_exp",-0.3,-4.,0.);
   RooExponential pdf_m_combinatorial_exp("pdf_m_combinatorial_exp","pdf_m_combinatorial_exp",Bmass,m_exp);
   pdf_m_signal = new RooAddPdf("pdf_m_signal","pdf_m_signal",RooArgList(m_gaussian1,m_gaussian2),RooArgList(m_fraction));

   //Two Exponentials
   RooRealVar m_exp2("m_exp2","m_exp2",-0.3,-4.,0.);
   RooExponential pdf_m_combinatorial_exp2("pdf_m_combinatorial_exp2","pdf_m_combinatorial_exp2",Bmass,m_exp2);
   RooRealVar m_fraction_exp("m_fraction_exp","m_fraction_exp", 0.5);

   pdf_m_combinatorial=new RooAddPdf("pdf_m_combinatorial","pdf_m_combinatorial",RooArgList(pdf_m_combinatorial_exp,pdf_m_combinatorial_exp2),RooArgList(m_fraction_exp));
   m_exp2.setConstant(kTRUE);
   m_fraction_exp.setVal(1.);

   RooRealVar n_signal("n_signal","n_signal",n_signal_initial,0.,data->sumEntries());
   RooRealVar n_combinatorial("n_combinatorial","n_combinatorial",n_combinatorial_initial,0.,data->sumEntries());

   RooAddPdf* model;
   model = new RooAddPdf("model","model", RooArgList(*pdf_m_signal, *pdf_m_combinatorial), RooArgList(n_signal, n_combinatorial));

   w.import(*model);

} // "BUILD PDF" ENDS

void set_up_workspace_variables(RooWorkspace& w)
{
  double mass_min, mass_max;
  double pt_min, pt_max;
  double y_min, y_max;
  double svpvDistance_min, svpvDistance_max;
  double svpvDistance2D_min, svpvDistance2D_max;
  double svpvDistancerr_min, svpvDistancerr_max;
  double svpvDistanc2Derr_min, svpvDistanc2Derr_max;
  double mu1pt_min, mu1pt_max;
  double mu1eta_min, mu1eta_max;
  double mu2pt_min, mu2pt_max;
  double mu2eta_min, mu2eta_max;
  double trk1D0_min, trk1D0_max;
  double trk2D0_min, trk2D0_max;
  double trk1D0err_min, trk1D0err_max;
  double trk2D0err_min, trk2D0err_max;
  double trk1Dz_min, trk1Dz_max;
  double trk2Dz_min, trk2Dz_max;
  double d0_min, d0_max;
  double d0err_min, d0err_max;
  double lxy_min, lxy_max;
  double pthatweight_min, pthatweight_max;
  double hibin_min, hibin_max;
  double alpha_min, alpha_max;
  double trk1pt_min, trk1pt_max;
  double trk2pt_min, trk2pt_max;
  double trk1eta_min, trk1eta_max;
  double trk2eta_min, trk2eta_max;
  double dtheta_min, dtheta_max;
  double chi2cl_min, chi2cl_max;

  svpvDistance2D_min = 0.0;
  svpvDistance2D_max = 1.0;

  svpvDistanc2Derr_min = 0;
  svpvDistanc2Derr_max = 0.035;

  mu1pt_min = 0.;
  mu1pt_max = 20.;

  mu2pt_min = 0.;
  mu2pt_max = 20.;

  mu1eta_min = -2.;
  mu1eta_max = 2.;

  mu2eta_min = -2.;
  mu2eta_max = 2.;

  trk1D0_min = -0.3;
  trk1D0_max = 0.3;

  trk1D0err_min = 0.;
  trk1D0err_max = 0.01;

  trk2D0err_min = 0.;
  trk2D0err_max = 0.01;

  trk2D0_min = -0.3;
  trk2D0_max = 0.3;

  trk1Dz_min = -0.1;
  trk1Dz_max = 0.1;

  trk2Dz_min = -0.1;
  trk2Dz_max = 0.1;

  d0_min = 0.;
  d0_max = 0.4;

  d0err_min = 0.;
  d0err_max = 0.0001;

  lxy_min = 0.;
  lxy_max = 3.;

  pthatweight_min = -0.5;
  pthatweight_max = 0.5;

  hibin_min = 0;
  hibin_max = 200;

  alpha_min = 0.;
  alpha_max = 0.1;

  pt_min= 5;
  pt_max= 100;

  y_min=-2.4;
  y_max=2.4;

  mass_min=5.14;
  mass_max=5.6;

  svpvDistance_min = 0.;
  svpvDistance_max= 2.;

  svpvDistancerr_min = 0.;
  svpvDistancerr_max = 0.05;

  dtheta_min = 0;
  dtheta_max = 3.1415;

  trk1pt_min = 0.;
  trk1pt_max = 8.;

  trk2pt_min = 0.;
  trk2pt_max = 8.;

  trk1eta_min = -3.;
  trk1eta_max = 3.;

  trk2eta_min = -3.;
  trk2eta_max = 3.;

  chi2cl_min = 0.;
  chi2cl_max = 1.;


   RooRealVar Bmass("Bmass","Bmass",mass_min,mass_max);
   RooRealVar Bpt("Bpt","Bpt",pt_min,pt_max);
   RooRealVar By("By","By",y_min,y_max);
   RooRealVar Bdtheta ("Bdtheta","Bdtheta",dtheta_min,dtheta_max);
   RooRealVar BsvpvDistance ("BsvpvDistance","BsvpvDistance",svpvDistance_min,svpvDistance_max);
   RooRealVar BsvpvDistance_Err ("BsvpvDistance_Err","BsvpvDistance_Err",svpvDistancerr_min,svpvDistancerr_max);
   RooRealVar Btrk1pt ("Btrk1pt","Btrk1pt",trk1pt_min,trk1pt_max);
   RooRealVar Btrk2pt ("Btrk2pt","Btrk2pt",trk2pt_min,trk2pt_max);
   RooRealVar Bchi2cl ("Bchi2cl","Bchi2cl",chi2cl_min,chi2cl_max);
   RooRealVar Btrk1eta ("Btrk1eta","Btrk1eta",trk1eta_min,trk1eta_max);
   RooRealVar Btrk2eta ("Btrk2eta","Btrk2eta",trk2eta_min,trk2eta_max);
   RooRealVar BsvpvDistance_2D ("BsvpvDistance_2D","BsvpvDistance_2D",svpvDistance2D_min,svpvDistance2D_max);
   RooRealVar BsvpvDistance_2D_Err ("BsvpvDistance_2D_Err","BsvpvDistance_2D_Err",svpvDistanc2Derr_min,svpvDistanc2Derr_max);
   RooRealVar Balpha ("Balpha","Balpha",alpha_min,alpha_max);
   RooRealVar Bmu1pt ("Bmu1pt","Bmu1pt",mu1pt_min,mu1pt_max);
   RooRealVar Bmu2pt ("Bmu2pt","Bmu2pt",mu2pt_min,mu2pt_max);
   RooRealVar Bmu1eta ("Bmu1eta","Bmu1eta",mu1eta_min,mu1eta_max);
   RooRealVar Bmu2eta ("Bmu2eta","Bmu2eta",mu2eta_min,mu2eta_max);
   RooRealVar Bd0 ("Bd0","Bd0",d0_min,d0_max);
   RooRealVar Bd0err ("Bd0err","Bd0err",d0err_min,d0err_max);
   RooRealVar Btrk1D0 ("Btrk1D0","Btrk1D0",trk1D0_min,trk1D0_max);
   RooRealVar Btrk2D0 ("Btrk2D0","Btrk2D0",trk2D0_min,trk2D0_max);
   RooRealVar Btrk1D0Err ("Btrk1D0Err","Btrk1D0Err",trk1D0err_min,trk1D0err_max);
   RooRealVar Btrk2D0Err ("Btrk2D0Err","Btrk2D0Err",trk2D0err_min,trk2D0err_max);
   RooRealVar Btrk1Dz ("Btrk1Dz","Btrk1Dz",trk1Dz_min,trk1Dz_max);
   RooRealVar Btrk2Dz ("Btrk2Dz","Btrk2Dz",trk2Dz_min,trk2Dz_max);
   RooRealVar Blxy ("Blxy","Blxy",lxy_min,lxy_max);
   RooRealVar Pthatweight ("Pthatweight","Pthatweight",pthatweight_min,pthatweight_max);
   RooRealVar HiBin ("HiBin","HiBin",hibin_min,hibin_max);

    w.import(Bmass);
    w.import(Bdtheta);
    w.import(BsvpvDistance);
    w.import(BsvpvDistance_Err);
    w.import(Btrk1pt);
    w.import(Btrk2pt);
    w.import(Bpt);
    w.import(Bchi2cl);
    w.import(Btrk1eta);
    w.import(Btrk2eta);
    w.import(By);
    w.import(BsvpvDistance_2D);
    w.import(BsvpvDistance_2D_Err);
    w.import(Bmu1pt);
    w.import(Bmu2pt);
    w.import(Bmu1eta);
    w.import(Bmu2eta);
    w.import(Btrk1D0);
    w.import(Btrk2D0);
    w.import(Btrk1D0Err);
    w.import(Btrk2D0Err);
    w.import(Btrk1Dz);
    w.import(Btrk2Dz);
    w.import(Bd0);
    w.import(Bd0err);
    w.import(Blxy);
    w.import(Balpha);
    w.import(Pthatweight);
    w.import(HiBin);


  } //"SET UP WS VARIABLES" ENDS
