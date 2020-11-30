#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TEfficiency.h>
#include <TGraphErrors.h>
#include <iostream>
#include <TMultiGraph.h>

void frag_ratio(){
 TFile* f_syst_Bs = new TFile("~/work2/BinQGP/results/Bs/x_section/x_section_syst.root");
 TFile* f_stat_Bs = new TFile("~/work2/BinQGP/results/Bs/x_section/x_section_stat.root");

 TFile* f_syst_Bu = new TFile("~/work2/BinQGP/results/Bu/x_section/x_section_syst.root");
 TFile* f_stat_Bu = new TFile("~/work2/BinQGP/results/Bu/x_section/x_section_stat.root");


 TGraphAsymmErrors* syst_Bs = (TGraphAsymmErrors*)f_syst_Bs->Get("Graph;1");
 TGraphAsymmErrors* stat_Bs = (TGraphAsymmErrors*)f_stat_Bs->Get("Graph;1");

 TGraphAsymmErrors* syst_Bu = (TGraphAsymmErrors*)f_syst_Bu->Get("Graph;1");
 TGraphAsymmErrors* stat_Bu = (TGraphAsymmErrors*)f_stat_Bu->Get("Graph;1");

 const int n_pt_bins = syst_Bs->GetN();
 double pt_bins[] = {5,15,20,30,50};

 //Cross-section values
 double *xsec_Bs = syst_Bs->GetY();
 double *xsec_Bu = syst_Bu->GetY();

 //Pt error
 double *pt_error_low = stat_Bs->GetEXlow();
 double *pt_error_high = stat_Bs->GetEXhigh();
 double pt_zero[n_pt_bins]; 
 double *pt_values = stat_Bs->GetX();

 //Syst error values
 double syst_Bs_low;
 double syst_Bs_high;
 double syst_Bu_low;
 double syst_Bu_high;

 //Stat error values
 double stat_Bs_low;
 double stat_Bs_high;
 double stat_Bu_low;
 double stat_Bu_high;

 double ratio[n_pt_bins];
 double ratio_syst_low[n_pt_bins];
 double ratio_syst_high[n_pt_bins];
 double ratio_stat_low[n_pt_bins]; 
 double ratio_stat_high[n_pt_bins];

 double Dx;
 double Dy;

 for(int i=0; i<n_pt_bins; i++){
   stat_Bs_low = stat_Bs->GetErrorYlow(i);
   stat_Bs_high = stat_Bs->GetErrorYhigh(i);
   syst_Bs_low = syst_Bs->GetErrorYlow(i);
   syst_Bs_high = syst_Bs->GetErrorYhigh(i);

   stat_Bu_low = stat_Bu->GetErrorYlow(i);
   stat_Bu_high = stat_Bu->GetErrorYhigh(i);
   syst_Bu_low = syst_Bu->GetErrorYlow(i);
   syst_Bu_high = syst_Bu->GetErrorYhigh(i);

   pt_zero[i] = 0.;

   ratio[i] = xsec_Bs[i]/xsec_Bu[i];

   Dx = 1/xsec_Bu[i];
   Dy = xsec_Bs[i]/(xsec_Bu[i]*xsec_Bu[i]);
 
   
   ratio_syst_low[i] = sqrt(pow(Dx,2)*pow(syst_Bs_low,2)+pow(Dy,2)*pow(syst_Bu_low,2));   
   ratio_syst_high[i] = sqrt(pow(Dx,2)*pow(syst_Bs_high,2)+pow(Dy,2)*pow(syst_Bu_high,2));
   ratio_stat_low[i] = sqrt(pow(Dx,2)*pow(stat_Bs_low,2)+pow(Dy,2)*pow(stat_Bu_low,2));
   ratio_stat_high[i] = sqrt(pow(Dx,2)*pow(stat_Bs_high,2)+pow(Dy,2)*pow(stat_Bu_high,2));
 }


  TCanvas c;
  TMultiGraph* mg = new TMultiGraph();

  TGraphAsymmErrors* g_stat = new TGraphAsymmErrors(n_pt_bins,pt_values,ratio,pt_error_low,pt_error_high,ratio_stat_low,ratio_stat_high);
  g_stat->SetTitle("");
  g_stat->SetMarkerColor(4);
  g_stat->SetMarkerStyle(2);
  g_stat->SetLineColor(1);

  TGraphAsymmErrors* g_syst= new TGraphAsymmErrors(n_pt_bins,pt_values,ratio,pt_zero,pt_zero,ratio_syst_low,ratio_syst_high);
  g_syst->SetTitle("");
  g_syst->SetMarkerColor(4);
  g_syst->SetMarkerStyle(2);
  g_syst->SetLineColor(2);

  mg->Add(g_stat);
  mg->Add(g_syst);

  TLegend *leg = new TLegend(0.7, 0.7, 0.9,0.9);
  leg->SetFillColor(0);
  leg->AddEntry(g_stat, "Statistical Uncertainty", "lp");
  leg->AddEntry(g_syst, "Systematic Uncertainty", "lp");

  mg->Draw("AP");
  leg->Draw();
  mg->GetXaxis()->SetTitle("p_{T}(B) [GeV]");
  mg->GetYaxis()->SetTitle("Fragmentation Ratio f_{s}/f_{u}");
  mg->GetYaxis()->SetRangeUser(0.,0.4);
 
  c.SaveAs("~/work2/BinQGP/results/frag_ratio.gif");
  c.SaveAs("~/work2/BinQGP/results/frag_ratio.pdf");
}
