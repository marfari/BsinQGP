#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TEfficiency.h>
#include <TGraphErrors.h>
#include <iostream>

using namespace std;

#define particle 1; //0 = B+;      1 = Bs;

int main(){
  
  TFile* f_eff0 = new TFile("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/results/Bu/efficiency/root_files/efficiency0.root");
  TFile* f_eff1 = new TFile("/home/t3cms/ev19u033/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/results/Bu/efficiency/root_files/efficiency1.root");

  double pt_bins[] = {5, 7, 10, 15, 20, 30, 50, 100};

  TEfficiency* efficiency0 = new TEfficiency("efficiency0", "efficiency0", 7, pt_bins);
  TEfficiency* efficiency1 = new TEfficiency("efficiency1", "efficiency1", 7, pt_bins);
  efficiency0 = (TEfficiency*)f_eff0->Get("hist_tot_noweights_clone");
  efficiency1 = (TEfficiency*)f_eff1->Get("hist_tot_weights_clone");

  double eff0;
  double eff1;
  double syst;

  double y_values[7];
  double x_values[] = {6, 8.5, 12.5, 17.5, 25, 40, 75};

  double x_errors[] = {1, 1.5, 2.5, 2.5, 5, 10, 25};
  double y_errors[] = {0, 0, 0, 0, 0, 0, 0};

  for(int i = 0; i < 7; i++)
    {
      eff0 = efficiency0->GetEfficiency(i + 1);
      eff1 = efficiency1->GetEfficiency(i + 1);
      syst = (eff1 - eff0) / eff0;
      cout << syst << endl;
      y_values[i] = syst;
    }

  cout << endl;
  cout << "Efficiency0, bin 7: " << efficiency0->GetEfficiency(7) << endl;
  cout << "Efficiency1, bin 7: " << efficiency1->GetEfficiency(7) << endl;

  TGraphErrors* systematic_errors = new TGraphErrors(7, x_values, y_values, x_errors, y_errors);
  TCanvas c;
  systematic_errors->SetMarkerColor(4);
  systematic_errors->SetMarkerStyle(5);
  systematic_errors->Draw("AP");
  systematic_errors->SetTitle("Bpt efficiency systematic error");
  c.SaveAs("./results/Bu/efficiency/plots_Bpt/systematic_error.gif");
  c.SaveAs("./results/Bu/efficiency/plots_Bpt/systematic_error.pdf");
  

  TFile* f1 = new TFile("./results/Bu/efficiency/root_files_Bpt/efficiency_systematic_errors.root", "recreate");
  f1->cd();
  systematic_errors->Write();
  f1->Write();
  f1->ls();
  f1->Close();
  
}
