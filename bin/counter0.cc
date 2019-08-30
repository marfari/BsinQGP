#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TEfficiency.h>
#include <TString.h>
#include <iostream>


using namespace std;

double read_weights(TString var, double var_value);
double getWeight(double var_value, TH1D* h_weight);

#define particle 0; //0 = B+;   1 = Bs;

int main(){
  //int counter(){

  //TString input_nocuts = particle ? "./prefiltered_trees_2/selected_mc_ntphi_PbPb_2018_corrected_nocuts_BDT.root" : "./prefiltered_trees_2/selected_mc_ntKp_PbPb_2018_corrected_nocuts_BDT.root";
  //TString input_cuts = particle ? "/home/t3cms/julia/LSTORE/CMSSW_7_5_8_patch5/src/UserCode/Bs_analysis/prefiltered_trees/selected_mc_ntphi_PbPb_2018_corrected_test_new.root" : "/home/t3cms/julia/LSTORE/CMSSW_7_5_8_patch5/src/UserCode/Bs_analysis/prefiltered_trees/selected_mc_ntKp_PbPb_2018_corrected_test.root";

  //TString input_t = particle ? "ntphi" : "ntKp";

  TFile* f_mc_nocuts = new TFile("./prefiltered_trees_2/selected_mc_ntKp_PbPb_2018_corrected_nocuts_BDT.root");
  TTree* t_nocuts = (TTree*)f_mc_nocuts->Get("ntKp");

  TFile* f_mc_cuts = new TFile("/home/t3cms/julia/LSTORE/CMSSW_7_5_8_patch5/src/UserCode/Bs_analysis/prefiltered_trees/selected_mc_ntKp_PbPb_2018_corrected_test.root");
  TTree* t_cuts = (TTree*)f_mc_cuts->Get("ntKp");

  double pt_bins[] = {5, 7, 10, 15, 20, 30, 50, 100};

  TH1F* hist_passed_noweights = new TH1F("hist_passed_noweights", "hist_passed_noweights", 7, pt_bins);
  TH1F* hist_tot_noweights = new TH1F("hist_tot_noweights", "hist_tot_noweights", 7, pt_bins);

  TH1F* hist_passed_weights = new TH1F("hist_passed_weights", "hist_passed_weights", 7, pt_bins);
  TH1F* hist_tot_weights = new TH1F("hist_tot_weights", "hist_tot_weights", 7, pt_bins);

  float bpt1;
  t_cuts->SetBranchAddress("Bpt", &bpt1);
  double weight;

  for(int evt = 0; evt < t_cuts->GetEntries(); evt++)
    {
      t_cuts->GetEntry(evt);
      hist_passed_noweights->Fill(bpt1);
      weight = read_weights("Bpt", bpt1);
      hist_passed_weights->Fill(bpt1, weight);
    }

  TCanvas passed_noweights;
  hist_passed_noweights->Draw();
  passed_noweights.SaveAs("./bin/results/efficiency/plots/passed_noweights.pdf");
  TCanvas passed_weights;
  hist_passed_weights->Draw();
  passed_weights.SaveAs("./bin/results/efficiency/plots/passed_weights.pdf");
  
  float bpt2;
  t_nocuts->SetBranchAddress("Bpt", &bpt2);

  for(int evt = 0; evt < t_nocuts->GetEntries(); evt++)
    {
      t_nocuts->GetEntry(evt);
      hist_tot_noweights->Fill(bpt2);
      weight = read_weights("Bpt", bpt2);
      hist_tot_weights->Fill(bpt2, weight);
    }

  TCanvas total_noweights;
  hist_tot_noweights->Draw();
  total_noweights.SaveAs("./bin/results/efficiency/plots/total_noweights.pdf");
  TCanvas total_weights;
  hist_tot_weights->Draw();
  total_weights.SaveAs("./bin/results/efficiency/plots/total_weights.pdf");

  TEfficiency* efficiency0 = new TEfficiency(*hist_passed_noweights, *hist_tot_noweights);
  TEfficiency* efficiency1 = new TEfficiency(*hist_passed_weights, *hist_tot_weights);

  TCanvas c0;
  efficiency0->Draw("AP");
  c0.SaveAs("./bin/results/efficiency/plots/efficiency0.pdf");

  TFile* f0 = new TFile("./bin/results/efficiency/root_files/efficiency0.root" , "recreate");
  f0->cd();
  efficiency0->Write();
  f0->Write();

  TCanvas c1;
  efficiency1->Draw("AP");
  c1.SaveAs("./bin/results/efficiency/plots/efficiency1.pdf");

  TFile* f1 = new TFile("./bin/results/efficiency/root_files/efficiency1.root" , "recreate");
  f1->cd();
  efficiency1->Write();
  f1->Write();

  for(int i = 1; i < 8; i++)
    {
      cout << efficiency0->GetEfficiency(i) << endl;
    }

  cout << endl;

  for(int i = 1; i < 8; i++)
    {
      cout << efficiency1->GetEfficiency(i) << endl;
    }

  f_mc_nocuts->Close();
  f_mc_cuts->Close();
  delete f_mc_nocuts;
  delete f_mc_cuts;

  f0->ls();
  f0->Close();
  f1->ls();
  f1->Close();
  
  return 0;
  
}

double read_weights(TString variable, double var_value){
  
  TString input_file = "./bin/results/B+/mc_validation_plots/weights/weights.root";

  TFile* f_wei = new TFile(input_file, "read");

  TH1D* histo_variable = (TH1D*)f_wei->Get(Form("weights_"+variable));

  double weight;
  double variable_min;
  double variable_max;

  variable_min = histo_variable->GetXaxis()->GetXmin();
  variable_max = histo_variable->GetXaxis()->GetXmax();  
  
  //if the event is not in the range its weight is 1.
  if(var_value>=variable_min && var_value<=variable_max){  
    weight = getWeight(var_value,histo_variable);
  }
  else{
    weight = 1;
  }

  f_wei->Close();
  delete f_wei;

  return weight;
}

//definição da fc auxiliar

double getWeight(double var_value, TH1D* h_weight){
  int bin = h_weight->FindBin(var_value);
  return h_weight->GetBinContent(bin);
}
