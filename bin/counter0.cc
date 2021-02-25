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

#define particle 1 //0 = B+;   1 = Bs;

int counter0(){
  TFile* f_raw_yield = particle ?  new TFile("~/public/BinQGP/results/Bs/Bpt/pT.root") : new TFile("~/public/BinQGP/results/Bu/Bpt/pT.root");
  TGraphAsymmErrors* raw_yield = (TGraphAsymmErrors*)f_raw_yield->Get("Graph;1");

  TString input_f_mc_cuts = particle ? "/eos/cms/store/group/phys_heavyions/zshi/SamplesForMaria/NewBDTCut/Bs/BsMC.root" : "/eos/cms/store/group/phys_heavyions/zshi/SamplesForMaria/NewBDTCut/BP/BPMC.root";
  TFile* f_mc_cuts = new TFile(input_f_mc_cuts);
  
  TString input_f_mc_nocuts = particle ? "/eos/cms/store/group/phys_heavyions/zshi/SamplesForMaria/NewBDTCut/Bs/BsMCNoCut.root" : "/eos/cms/store/group/phys_heavyions/zshi/SamplesForMaria/NewBDTCut/BP/BPMCNoCut.root";
  TFile* f_mc_nocuts = new TFile(input_f_mc_nocuts);

  TString input_t = particle ? "ntphi" : "ntKp";

  TTree* t_cuts = (TTree*)f_mc_cuts->Get(input_t);
  TTree* t_nocuts = (TTree*)f_mc_nocuts->Get("Bfinder/ntphi");
  TTree* t_gen = (TTree*)f_mc_nocuts->Get("Bfinder/ntGen");

  //BDT tree (no cuts)  
  TTree* t_bdt_pt_3_5 = (TTree*)f_mc_nocuts->Get("BDT_pt_3_5");
  TTree* t_bdt_pt_5_7 = (TTree*)f_mc_nocuts->Get("BDT_pt_5_7");
  TTree* t_bdt_pt_7_10 = (TTree*)f_mc_nocuts->Get("BDT_pt_7_10");
  TTree* t_bdt_pt_10_15 = (TTree*)f_mc_nocuts->Get("BDT_pt_10_15");
  TTree* t_bdt_pt_15_20 = (TTree*)f_mc_nocuts->Get("BDT_pt_15_20");
  TTree* t_bdt_pt_20_50 = (TTree*)f_mc_nocuts->Get("BDT_pt_20_50");
  TTree* t_bdt_pt_50_100 = (TTree*)f_mc_nocuts->Get("BDT_pt_50_100");


  double pt_bins[] = {3,5,7,10,15,20,50,100};
  const int n_pt_bins = raw_yield->GetN();

  TH1F* hist_tot_noweights = new TH1F("hist_tot_noweights", "hist_tot_noweights", n_pt_bins, pt_bins);
  TH1F* hist_passed_noweights = new TH1F("hist_passed_noweights", "hist_passed_noweights", n_pt_bins, pt_bins);

  TH1F* hist_tot_weights = new TH1F("hist_tot_weights", "hist_tot_weights", n_pt_bins, pt_bins);
  TH1F* hist_passed_weights = new TH1F("hist_passed_weights", "hist_passed_weights", n_pt_bins, pt_bins);

  TH1F* hist_tot_gen = new TH1F("hist_tot_gen", "hist_tot_gen", n_pt_bins, pt_bins);

  //NOCUTS 
  float bpt1;
  //float pthat_nocuts;
  //float weight_nocuts;

  t_nocuts->SetBranchAddress("Bpt", &bpt1);
  //t_nocuts->SetBranchAddress("pthat", &pthat_nocuts);
  //t_nocuts->SetBranchAddress("weight", &weight_nocuts);

  double bdt_pt_3_5;
  double bdt_pt_5_7;
  double bdt_pt_7_10;
  double bdt_pt_10_15;
  double bdt_pt_15_20;
  double bdt_pt_20_50;
  double bdt_pt_50_100;

  t_bdt_pt_3_5->SetBranchAddress("BDT_pt_3_5", &bdt_pt_3_5);
  t_bdt_pt_5_7->SetBranchAddress("BDT_pt_5_7", &bdt_pt_5_7);
  t_bdt_pt_7_10->SetBranchAddress("BDT_pt_7_10", &bdt_pt_7_10);
  t_bdt_pt_10_15->SetBranchAddress("BDT_pt_10_15", &bdt_pt_10_15);
  t_bdt_pt_15_20->SetBranchAddress("BDT_pt_15_20", &bdt_pt_15_20);
  t_bdt_pt_20_50->SetBranchAddress("BDT_pt_20_50", &bdt_pt_20_50);
  t_bdt_pt_50_100->SetBranchAddress("BDT_pt_50_100", &bdt_pt_50_100);

  double weight = 1;

  TString variable;
  double bdt1_total = 0;

  for(int evt = 0; evt < t_nocuts->GetEntries(); evt++)
    {
      t_nocuts->GetEntry(evt);

      for(int kk=0; kk<n_pt_bins; kk++){
	if ( (bpt1 < pt_bins[kk]) || (bpt1 > pt_bins[kk+1]) )
	  continue;
        cout << "pt_bins" << pt_bins[kk] << "pt_bins + 1 = " << pt_bins[kk+1] << endl;
	variable.Form("BDT_pt_%g_%g", pt_bins[kk], pt_bins[kk+1]);
	if ((3<bpt1) && (bpt1<5))
	  {bdt1_total = bdt_pt_3_5;}
	else if ((5<bpt1) && (bpt1<7))
	  {bdt1_total = bdt_pt_5_7;}
        else if ((7<bpt1) && (bpt1<10))
          {bdt1_total = bdt_pt_7_10;}
        else if ((10<bpt1) && (bpt1<15))
          {bdt1_total = bdt_pt_10_15;}
	else if ((15<bpt1) && (bpt1<20))
	  {bdt1_total = bdt_pt_15_20;}
	else if ((20<bpt1) && (bpt1<50))
	  {bdt1_total = bdt_pt_20_50;}
        else if ((50<bpt1) && (bpt1<100))
          {bdt1_total = bdt_pt_50_100;}
	weight = read_weights(variable, bdt1_total);
        //weight *= pthat_nocuts*weight_nocuts;
      }      
      hist_tot_weights->Fill(bpt1, weight);
      hist_tot_noweights->Fill(bpt1);
    }
 
  
   
  if(particle == 0){
    TCanvas tot_noweights;
    hist_tot_noweights->Draw();
    tot_noweights.SaveAs("~/public/BinQGP/results/Bu/efficiency/plots/tot_noweights.pdf");
    TCanvas tot_weights;
    hist_tot_weights->Draw();
    tot_weights.SaveAs("~/public/BinQGP/results/Bu/efficiency/plots/totweights.pdf");
  }
  else if(particle == 1){
    TCanvas tot_noweights;
    hist_tot_noweights->Draw();
    tot_noweights.SaveAs("~/public/BinQGP/results/Bs/efficiency/plots/tot_noweights.pdf");
    TCanvas tot_weights;
    hist_tot_weights->Draw();
    tot_weights.SaveAs("~/public/BinQGP/results/Bs/efficiency/plots/totweights.pdf");
  }


  //CUTS
  float bpt2;
  //float pthat_cuts;
  //float weight_cuts;

  t_cuts->SetBranchAddress("Bpt", &bpt2);
  //t_cuts->SetBranchAddress("pthat", &pthat_cuts);
  //t_cuts->SetBranchAddress("weight", &weight_cuts);

  double bdt2_pt_3_5;
  double bdt2_pt_5_7;
  double bdt2_pt_7_10;
  double bdt2_pt_10_15;
  double bdt2_pt_15_20;
  double bdt2_pt_20_50;
  double bdt2_pt_50_100;
 
  t_cuts->SetBranchAddress("BDT_pt_3_5", &bdt2_pt_3_5);
  t_cuts->SetBranchAddress("BDT_pt_5_7", &bdt2_pt_5_7);
  t_cuts->SetBranchAddress("BDT_pt_7_10", &bdt2_pt_7_10);
  t_cuts->SetBranchAddress("BDT_pt_10_15", &bdt2_pt_10_15);
  t_cuts->SetBranchAddress("BDT_pt_15_20", &bdt2_pt_15_20);
  t_cuts->SetBranchAddress("BDT_pt_20_50", &bdt2_pt_20_50);
  t_cuts->SetBranchAddress("BDT_pt_50_100", &bdt2_pt_50_100);
 
  double weight2 = 1;
  double bdt2_total = 0;

  TString variable2;

  for(int evt = 0; evt < t_cuts->GetEntries(); evt++)
  {
    t_cuts->GetEntry(evt);
       
    for(int kk=0; kk<n_pt_bins; kk++){
       if ( (bpt2 < pt_bins[kk]) || (bpt2 > pt_bins[kk+1]) )
	  continue;
	variable2.Form("BDT_pt_%g_%g", pt_bins[kk], pt_bins[kk+1]);
	if ((3<bpt2) && (bpt2<5))
	  {bdt2_total = bdt2_pt_3_5;}
	else if ((5<bpt2) && (bpt2<7))
	  {bdt2_total = bdt2_pt_5_7;}
        else if ((7<bpt2) && (bpt2<10))
          {bdt2_total = bdt2_pt_7_10;}
        else if ((10<bpt2) && (bpt2<15))
          {bdt2_total = bdt2_pt_10_15;}
	else if ((15<bpt2) && (bpt2<20))
	  {bdt2_total = bdt2_pt_15_20;}
	else if ((20<bpt2) && (bpt2<50))
	  {bdt2_total = bdt2_pt_20_50;}
        else if ((50<bpt2) && (bpt2<100))
          {bdt2_total = bdt2_pt_50_100;}
	weight2 = read_weights(variable2, bdt2_total);
        //weight2 *= pthat_cuts*weight_cuts;
      }
      
      hist_passed_weights->Fill(bpt2, weight2);
      hist_passed_noweights->Fill(bpt2);
    }


  if(particle == 0){
    TCanvas passed_noweights;
    hist_passed_noweights->Draw();
    passed_noweights.SaveAs("~/public/BinQGP/results/Bu/efficiency/plots/passed_noweights.pdf");
    TCanvas passed_weights;
    hist_passed_weights->Draw();
    passed_weights.SaveAs("~/public/BinQGP/results/Bu/efficiency/plots/passed_weights.pdf");
  }else if(particle == 1){
    TCanvas passed_noweights;
    hist_passed_noweights->Draw();
    passed_noweights.SaveAs("~/public/BinQGP/results/Bs/efficiency/plots/passed_noweights.pdf");
    TCanvas passed_weights;
    hist_passed_weights->Draw();
    passed_weights.SaveAs("~/public/BinQGP/results/Bs/efficiency/plots/passed_weights.pdf");
  }

  //GEN
  
  float bpt3;
  
  t_gen->SetBranchAddress("Gpt", &bpt3);
  
  for(int evt = 0; evt < t_gen->GetEntries(); evt++){
    t_gen->GetEntry(evt);
  
    hist_tot_gen->Fill(bpt3);
  }
  
  
  if(particle == 0){
    TCanvas tot_gen;
    hist_tot_gen->Draw();
    tot_gen.SaveAs("~/public/BinQGP/results/Bu/efficiency/plots/totgen.pdf");
  }
  else if(particle == 1){
    TCanvas tot_gen;
    hist_tot_gen->Draw();
    tot_gen.SaveAs("~/public/BinQGP/results/Bs/efficiency/plots/totgen.pdf");
  }

  //efficiency0 = (gen + cuts) / gen (without weights)
  TEfficiency* efficiency0 = new TEfficiency(*hist_passed_noweights, *hist_tot_noweights);
  if(particle == 0){
    TCanvas c0;
    efficiency0->SetTitle("Nominal Efficiency #epsilon^{0};p_{T} (GeV);#epsilon^{0}");
    efficiency0->Draw("AP");
    c0.SaveAs("~/public/BinQGP/results/Bu/efficiency/plots/efficiency0.pdf");
  }else if(particle == 1){
    TCanvas c0;
    efficiency0->SetTitle("Nominal Efficiency #epsilon^{0};p_{T} (GeV);#epsilon^{0}");
    efficiency0->Draw("AP");
    c0.SaveAs("~/public/BinQGP/results/Bs/efficiency/plots/efficiency0.pdf");
  }
  
  //efficiency1 = (gen + cuts) /gen (with weights)
  TEfficiency* efficiency1 = new TEfficiency(*hist_passed_weights, *hist_tot_weights);
  if(particle == 0){
    TCanvas c1;
    efficiency1->SetTitle("#epsilon^{1};p_{T} (GeV);#epsilon^{1}");
    efficiency1->Draw("AP");
    c1.SaveAs("~/public/BinQGP/results/Bu/efficiency/plots/efficiency1.pdf");
  }else if(particle == 1){
    TCanvas c1;
    efficiency1->SetTitle("#epsilon^{1};p_{T} (GeV);#epsilon^{1}");
    efficiency1->Draw("AP");
    c1.SaveAs("~/public/BinQGP/results/Bs/efficiency/plots/efficiency1.pdf");
  }
  
  // efficiency x acceptance = (gen + cuts) / gen (without weights)
  TEfficiency* eff_x_acc = new TEfficiency(*hist_passed_noweights, *hist_tot_gen);
  if(particle == 0){
    TCanvas c3;
    eff_x_acc->SetTitle("Efficiency x Acceptance");
    eff_x_acc->Draw("AP");
    c3.SaveAs("~/public/BinQGP/results/Bu/efficiency/plots/eff_x_acc.pdf");
  }else if(particle == 1){
    TCanvas c3;
    eff_x_acc->SetTitle("Efficiency x Acceptance");
    eff_x_acc->Draw("AP");
    c3.SaveAs("~/public/BinQGP/results/Bs/efficiency/plots/eff_x_acc.pdf");
  }

  if(particle == 0){
    TFile* f0 = new TFile("~/public/BinQGP/results/Bu/efficiency/root_files/efficiency0.root" , "recreate");
    f0->cd();
    efficiency0->Write();
    f0->Write();
    f0->ls();
    f0->Close();
    
    TFile* f1 = new TFile("~/public/BinQGP/results/Bu/efficiency/root_files/efficiency1.root" , "recreate");
    f1->cd();
    efficiency1->Write();
    f1->Write();
    f1->ls();
    f1->Close();
 
    TFile* f2 = new TFile("~/public/BinQGP/results/Bu/efficiency/root_files/eff_x_acc.root", "recreate");
    f2->cd();
    eff_x_acc->Write();
    f2->Write();
    f2->ls();
    f2->Close();

 
  }else if(particle == 1){
    TFile* f0 = new TFile("~/public/BinQGP/results/Bs/efficiency/root_files/efficiency0.root" , "recreate");
    f0->cd();
    efficiency0->Write();
    f0->Write();
    f0->ls();
    f0->Close();
     
    TFile* f1 = new TFile("~/public/BinQGP/results/Bs/efficiency/root_files/efficiency1.root" , "recreate");
    f1->cd();
    efficiency1->Write();
    f1->Write();
    f1->ls();
    f1->Close();
   
    TFile* f2 = new TFile("~public/BinQGP/results/Bs/efficiency/root_files/eff_x_acc.root", "recreate");
    f2->cd();
    eff_x_acc->Write();
    f2->Write();
    f2->ls();
    f2->Close();
  }

  delete hist_tot_noweights;
  delete hist_passed_noweights;
  delete hist_tot_weights;
  delete hist_passed_weights;
  delete hist_tot_gen;

  f_mc_cuts->Close();
  delete f_mc_cuts;
  f_mc_nocuts->Close();
  delete f_mc_nocuts;
  
  return 0;
  
}

double read_weights(TString variable, double var_value){
  
  TString input_file = particle ? "~/public/BinQGP/results/Bs/mc_validation_plots/weights/weights.root" :"~/public/BinQGP/results/Bu/mc_validation_plots/weights/weights.root";

  TFile* f_wei = new TFile(input_file, "read");

  TH1D* histo_variable = (TH1D*)f_wei->Get("weights_"+variable);

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
