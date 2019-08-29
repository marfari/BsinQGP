#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TEfficiency.h>
#include <iostream>

using namespace std;

double read_weights(TString var, double var_value);
double getWeight(double var_value, TH1D* h_weight);
TEfficiency* getEfficiency(double* pt_bins, int n_pt_bins, TString input_cuts, TString input_nocuts, bool weights);

int main(){
  //int counter(){

  TString input_cuts = "/home/t3cms/julia/LSTORE/CMSSW_7_5_8_patch5/src/UserCode/Bs_analysis/prefiltered_trees/selected_mc_ntKp_PbPb_2018_corrected_test.root";
  TString input_nocuts = "/home/t3cms/julia/LSTORE/CMSSW_7_5_8_patch5/src/UserCode/Bs_analysis/prefiltered_trees/selected_mc_ntKp_PbPb_2018_corrected_test_train_nocuts.root";
  
  double pt_bins[] = {5, 7, 10, 15, 20, 30, 50, 100};
  double n_pt_bins = 7;
  
  TEfficiency* efficiency0 = getEfficiency(pt_bins, n_pt_bins, input_cuts, input_nocuts, 0);
  TEfficiency* efficiency1 = getEfficiency(pt_bins, n_pt_bins, input_cuts, input_nocuts, 1);

  
  for(int i = 1; i < 8; i++)
    {
      cout << efficiency0->GetEfficiency(i) << endl;
    }

  cout << endl;

  for(int i = 1; i < 8; i++)
    {
      cout << efficiency1->GetEfficiency(i) << endl;
    }
  
  return 0;
  
}

TEfficiency* getEfficiency(double* pt_bins, int n_pt_bins, TString input_cuts, TString input_nocuts, bool weights){

  TFile* f_mc_cuts = new TFile(input_cuts);
  TTree* t_cuts = (TTree*)f_mc_cuts->Get("ntKp");

  TFile* f_mc_nocuts = new TFile(input_nocuts);
  TTree* t_nocuts = (TTree*)f_mc_nocuts->Get("ntKp");

  TH1F* hist_passed = new TH1F("hist_passed", "hist_passed", n_pt_bins, pt_bins);
  TH1F* hist_tot = new TH1F("hist_tot", "hist_tot", n_pt_bins, pt_bins);

  float bpt1;
  t_cuts->SetBranchAddress("Bpt", &bpt1);
  double weight;

  for(int evt = 0; evt < t_cuts->GetEntries(); evt++)
    {
      t_cuts->GetEntry(evt);
      if(weights){
	weight = read_weights("Bpt", bpt1);
	hist_passed->Fill(bpt1, weight);
      }else{
	hist_passed->Fill(bpt1);
      }
    }

  if(weights){
    TCanvas passed_weights;
    hist_passed->Draw();
    passed_weights.SaveAs("./bin/results/efficiency_test/passed_weights.pdf");
  }else{
    TCanvas passed_noweights;
    hist_passed->Draw();
    passed_noweights.SaveAs("./bin/results/efficiency_test/passed_noweights.pdf");
  }

  float bpt2;
  t_nocuts->SetBranchAddress("Bpt", &bpt2);

  for(int evt = 0; evt < t_nocuts->GetEntries(); evt++)
    {
      t_nocuts->GetEntry(evt);
      if(weights){
	weight = read_weights("Bpt", bpt2);
	hist_tot->Fill(bpt2, weight);
      }else{
	hist_tot->Fill(bpt2);
      }
    }

  if(weights){
    TCanvas total_weights;
    hist_tot->Draw();
    total_weights.SaveAs("./bin/results/efficiency_test/total_weights.pdf");
  }else{
    TCanvas total_noweights;
    hist_tot->Draw();
    total_noweights.SaveAs("./bin/results/efficiency_test/total_noweights.pdf");
  }

  TEfficiency* efficiency = new TEfficiency(*hist_passed, *hist_tot);

  TCanvas c1;
  efficiency->Draw("AP");
  if(weights){
    c1.SaveAs("./bin/results/efficiency_test/efficiency1.pdf");
  }else{
    c1.SaveAs("./bin/results/efficiency_test/efficiency0.pdf");
  }

  f_mc_nocuts->Close();
  f_mc_cuts->Close();
  delete f_mc_nocuts;
  delete f_mc_cuts;
  
  return efficiency;
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

  //testing
  //cout<<"min:"<<variable_min<<endl;
  //cout<<"max:"<<variable_max<<endl;
  
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
