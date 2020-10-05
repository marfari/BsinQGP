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

  TString input_f_mc_cuts = particle ? "/lstore/cms/nuno/ppdata2017/000814/BsMC.root" : "/lstore/cms/nuno/ppdata2017/000814/BPMC.root";
  TFile* f_mc_cuts = new TFile(input_f_mc_cuts);
  
  TString input_t = particle ? "ntphi" : "ntKp";
  TTree* t_cuts = (TTree*)f_mc_cuts->Get(input_t); 

  TString input_f_mc_nocuts = particle ? "/lstore/cms/nuno/ppdata2017/000926_nocut_flat/BsMC.root" : "/lstore/cms/nuno/ppdata2017/000906_nocut/BPMC_nobdt.root";
  TFile* f_mc_nocuts = new TFile(input_f_mc_nocuts);

  TTree* t_nocuts = (TTree*)f_mc_nocuts->Get(input_t);

  TString input_f_accept = particle ? "/lstore/cms/nuno/ppdata2017/000923_nocut/BsMC.root" : "/lstore/cms/nuno/ppdata2017/000923_nocut/BsMC.root";
  TFile* f_accept = new TFile(input_f_accept);

  TTree* t_accept = (TTree*)f_accept->Get("Bfinder/ntGen");  


  double pt_bins[] = {5, 10, 15, 20, 50};
  //double pt_bins[] = {5, 7, 10, 15, 20, 30, 50, 100};
  //double pt_binsh[] = {5, 7, 10, 15, 20, 30, 40, 50, 60};
  
  int n_pt_bins = 4;
  //int n_pt_bins = 7;
  //int n_pt_binsh = 8;

  TH1F* hist_tot_noweights = new TH1F("hist_tot_noweights", "hist_tot_noweights", n_pt_bins, pt_bins);
  TH1F* hist_passed_noweights = new TH1F("hist_passed_noweights", "hist_passed_noweights", n_pt_bins, pt_bins);

  TH1F* hist_tot_weights = new TH1F("hist_tot_weights", "hist_tot_weights", n_pt_bins, pt_bins);
  TH1F* hist_passed_weights = new TH1F("hist_passed_weights", "hist_passed_weights", n_pt_bins, pt_bins);

  TH1F* hist_tot_gen = new TH1F("hist_tot_gen", "hist_tot_gen", n_pt_bins, pt_bins);

  float bpt1;
  /*
  double bdt_pt_5_7;
  double bdt_pt_7_10;
  double bdt_pt_10_15;
  double bdt_pt_15_20;
  double bdt_pt_20_30;
  double bdt_pt_30_50;
  double bdt_pt_50_100;
  */
  
  float bdt_pt_5_10;
  float bdt_pt_10_15;
  float bdt_pt_15_20;
  float bdt_pt_20_50;

  t_nocuts->SetBranchAddress("Bpt", &bpt1);
  /*
  t_nocuts->SetBranchAddress("BDT_pt_5_7", &bdt_pt_5_7);
  t_nocuts->SetBranchAddress("BDT_pt_7_10", &bdt_pt_7_10);
  t_nocuts->SetBranchAddress("BDT_pt_10_15", &bdt_pt_10_15);
  t_nocuts->SetBranchAddress("BDT_pt_15_20", &bdt_pt_15_20);
  t_nocuts->SetBranchAddress("BDT_pt_20_30", &bdt_pt_20_30);
  t_nocuts->SetBranchAddress("BDT_pt_30_50", &bdt_pt_30_50);
  t_nocuts->SetBranchAddress("BDT_pt_50_100", &bdt_pt_50_100);
  */
  
  t_nocuts->SetBranchAddress("BDT_pt_5_10", &bdt_pt_5_10);
  t_nocuts->SetBranchAddress("BDT_pt_10_15", &bdt_pt_10_15);
  t_nocuts->SetBranchAddress("BDT_pt_15_20", &bdt_pt_15_20);
  t_nocuts->SetBranchAddress("BDT_pt_20_50", &bdt_pt_20_50);
  

  double weight = 1;
  double bdt1_total = 0;

  /*
  for(int evt = 0; evt < t_nocuts->GetEntries(); evt++)
    {
      t_nocuts->GetEntry(evt);
      for(int kk=0; kk<n_pt_bins; kk++){
	if ( (bpt1 < pt_bins[kk]) || (bpt1 > pt_bins[kk+1]) )
	  continue;
	variable.Form("BDT_pt_%g_%g", pt_bins[kk], pt_bins[kk+1]);
	if ((5<bpt1) && (bpt1<7))
	  {bdt_total = bdt_pt_5_7;}

	else if ((7<bpt1) && (bpt1<10))
	  {bdt_total = bdt_pt_7_10;}

	else if ((10<bpt1) && (bpt1<15))
	  {bdt_total = bdt_pt_10_15;}

	else if ((15<bpt1) && (bpt1<20))
	  {bdt_total = bdt_pt_15_20;}

	else if ((20<bpt1) && (bpt1<30))
	  {bdt_total = bdt_pt_20_30;}

	else if ((30<bpt1) && (bpt1<50))
	  {bdt_total = bdt_pt_30_50;}

	else if ((50<bpt1) && (bpt1<100))
	  {bdt_total = bdt_pt_50_100;}
       
	weight = read_weights("Bmu1pt", Bmu1pt);
      }
      hist_tot_weights->Fill(bpt1, weight);
      hist_tot_noweights->Fill(bpt1);
    }
  cout << "for_loop" << endl;
  */

  //Bin by bin analysis of the BDT for Bs
  TString variable;

  for(int evt = 0; evt < t_nocuts->GetEntries(); evt++)
    {
      t_nocuts->GetEntry(evt);

      for(int kk=0; kk<n_pt_bins; kk++){
	if ( (bpt1 < pt_bins[kk]) || (bpt1 > pt_bins[kk+1]) )
	  continue;
	variable.Form("BDT_pt_%g_%g", pt_bins[kk], pt_bins[kk+1]);
	if ((5<bpt1) && (bpt1<10))
	  {bdt1_total = bdt_pt_5_10;}

	else if ((10<bpt1) && (bpt1<15))
	  {bdt1_total = bdt_pt_10_15;}

	else if ((15<bpt1) && (bpt1<20))
	  {bdt1_total = bdt_pt_15_20;}

	else if ((20<bpt1) && (bpt1<50))
	  {bdt1_total = bdt_pt_20_50;}
	weight = read_weights(variable, bdt1_total);
      }      
      hist_tot_weights->Fill(bpt1, weight);
      hist_tot_noweights->Fill(bpt1);
    }
  /*    
  for(int evt = 0; evt < t_nocuts->GetEntries(); evt++)
      {
      t_nocuts->GetEntry(evt);
      hist_tot_noweights->Fill(bpt1);
      weight = read_weights("Bpt", bpt1);
      hist_tot_weights->Fill(bpt1, weight);
      }
    */
    
  if(particle == 0)
    {
      TCanvas tot_noweights;
      hist_tot_noweights->Draw();
      tot_noweights.SaveAs("~/work2/BinQGP/results/Bu/efficiency/plots/tot_noweights.pdf");
      TCanvas tot_weights;
      hist_tot_weights->Draw();
      tot_weights.SaveAs("~/work2/BinQGP/results/Bu/efficiency/plots/totweights.pdf");
    }
  else if(particle == 1){
    TCanvas tot_noweights;
    hist_tot_noweights->Draw();
    tot_noweights.SaveAs("~/work2/BinQGP/results/Bs/efficiency/plots/tot_noweights.pdf");
    TCanvas tot_weights;
    hist_tot_weights->Draw();
    tot_weights.SaveAs("~/work2/BinQGP/results/Bs/efficiency/plots/totweights.pdf");
  }

  float bpt2;
  /*
  Double_t bdt2_pt_5_7;
  Double_t bdt2_pt_7_10;
  Double_t bdt2_pt_10_15;
  Double_t bdt2_pt_15_20;
  Double_t bdt2_pt_20_30;
  Double_t bdt2_pt_30_50;
  Double_t bdt2_pt_50_100;
  */
  
  float bdt2_pt_5_10;
  float bdt2_pt_10_15;
  float bdt2_pt_15_20;
  float bdt2_pt_20_50;
  

  t_cuts->SetBranchAddress("Bpt", &bpt2);
  /*
  t_cuts->SetBranchAddress("BDT_pt_5_7", &bdt2_pt_5_7);
  t_cuts->SetBranchAddress("BDT_pt_7_10", &bdt2_pt_7_10);
  t_cuts->SetBranchAddress("BDT_pt_10_15", &bdt2_pt_10_15);
  t_cuts->SetBranchAddress("BDT_pt_15_20", &bdt2_pt_15_20);
  t_cuts->SetBranchAddress("BDT_pt_20_30", &bdt2_pt_20_30);
  t_cuts->SetBranchAddress("BDT_pt_30_50", &bdt2_pt_30_50);
  t_cuts->SetBranchAddress("BDT_pt_50_100", &bdt2_pt_50_100);
  */
  
  t_cuts->SetBranchAddress("BDT_pt_5_10", &bdt2_pt_5_10);
  t_cuts->SetBranchAddress("BDT_pt_10_15", &bdt2_pt_10_15);
  t_cuts->SetBranchAddress("BDT_pt_15_20", &bdt2_pt_15_20);
  t_cuts->SetBranchAddress("BDT_pt_20_50", &bdt2_pt_20_50);
  

  double weight2 = 1;
  double bdt2_total = 0;
  /*
  for(int evt = 0; evt < t_cuts->GetEntries(); evt++)
    {
      t_cuts->GetEntry(evt);

      for(int kk=0; kk<n_pt_bins; kk++){
	if ( (bpt2 < pt_bins[kk]) || (bpt2 > pt_bins[kk+1]) )
	  continue;
	variable2.Form("BDT_pt_%g_%g", pt_bins[kk], pt_bins[kk+1]);
	if ((5<bpt2) && (bpt2<7))
	  {bdt2_total = bdt2_pt_5_7;}

	else if ((7<bpt2) && (bpt2<10))
	  {bdt2_total = bdt2_pt_7_10;}

	else if ((10<bpt2) && (bpt2<15))
	  {bdt2_total = bdt2_pt_10_15;}

	else if ((15<bpt2) && (bpt2<20))
	  {bdt2_total = bdt2_pt_15_20;}

	else if ((20<bpt2) && (bpt2<30))
	  {bdt2_total = bdt2_pt_20_30;}

	else if ((30<bpt2) && (bpt2<50))
	  {bdt2_total = bdt2_pt_30_50;}

	else if ((50<bpt2) && (bpt2<100))
	  {bdt2_total = bdt2_pt_50_100;}
	weight2 = read_weights("Bmu1pt", Bmu2pt);
      }
      hist_passed_weights->Fill(bpt2, weight2);
      hist_passed_noweights->Fill(bpt2);
    }
  */

  //Bin by bin analysis of the BDT for Bs
    TString variable2;

    for(int evt = 0; evt < t_cuts->GetEntries(); evt++)
    {
      t_cuts->GetEntry(evt);
        
      for(int kk=0; kk<n_pt_bins; kk++){
	if ( (bpt2 < pt_bins[kk]) || (bpt2 > pt_bins[kk+1]) )
	  continue;
	variable2.Form("BDT_pt_%g_%g", pt_bins[kk], pt_bins[kk+1]);
	if ((5<bpt2) && (bpt2<10))
	  {bdt2_total = bdt2_pt_5_10;}

	else if ((10<bpt2) && (bpt2<15))
	  {bdt2_total = bdt2_pt_10_15;}

	else if ((15<bpt2) && (bpt2<20))
	  {bdt2_total = bdt2_pt_15_20;}

	else if ((20<bpt2) && (bpt2<50))
	  {bdt2_total = bdt2_pt_20_50;}
	weight2 = read_weights(variable2, bdt2_total);
      }
      
      hist_passed_weights->Fill(bpt2, weight2);
      hist_passed_noweights->Fill(bpt2);
    }
  /*
  //Analysis of Bpt
  for(int evt = 0; evt < t_cuts->GetEntries(); evt++)
    {
      t_cuts->GetEntry(evt);
      hist_passed_noweights->Fill(bpt2);
      weight2 = read_weights("Bpt", bpt2);
      hist_passed_weights->Fill(bpt2, weight2);
      }
  */

  float bpt3;

  t_accept->SetBranchAddress("Gpt", &bpt3);

  for(int evt = 0; evt < t_accept->GetEntries(); evt++){
    t_accept->GetEntry(evt);

    hist_tot_gen->Fill(bpt3);
  }

  if(particle == 0)
    {
      TCanvas tot_gen;
      hist_tot_gen->Draw();
      tot_gen.SaveAs("~/work2/BinQGP/results/Bu/efficiency/plots/totgen.pdf");
    }
  else if(particle == 1){
    TCanvas tot_gen;
    hist_tot_gen->Draw();
    tot_gen.SaveAs("~/work2/BinQGP/results/Bs/efficiency/plots/totgen.pdf");
  }


  if(particle == 0){
    TCanvas passed_noweights;
    hist_passed_noweights->Draw();
    passed_noweights.SaveAs("~/work2/BinQGP/results/Bu/efficiency/plots/passed_noweights.pdf");
    TCanvas passed_weights;
    hist_passed_weights->Draw();
    passed_weights.SaveAs("~/work2/BinQGP/results/Bu/efficiency/plots/passed_weights.pdf");
  }else if(particle == 1){
    TCanvas passed_noweights;
    hist_passed_noweights->Draw();
    passed_noweights.SaveAs("~/work2/BinQGP/results/Bs/efficiency/plots/passed_noweights.pdf");
    TCanvas passed_weights;
    hist_passed_weights->Draw();
    passed_weights.SaveAs("~/work2/BinQGP/results/Bs/efficiency/plots/passed_weights.pdf");
  }

  TEfficiency* efficiency0 = new TEfficiency(*hist_passed_noweights, *hist_tot_noweights);
  if(particle == 0){
    TCanvas c0;
    efficiency0->Draw("AP");
    c0.SaveAs("~/work2/BinQGP/results/Bu/efficiency/plots/efficiency0.pdf");
  }else if(particle == 1){
    TCanvas c0;
    efficiency0->Draw("AP");
    c0.SaveAs("~/work2/BinQGP/results/Bs/efficiency/plots/efficiency0.pdf");
  }
  
  TEfficiency* efficiency1 = new TEfficiency(*hist_passed_weights, *hist_tot_weights);
  if(particle == 0){
    TCanvas c1;
    efficiency1->Draw("AP");
    c1.SaveAs("~/work2/BinQGP/results/Bu/efficiency/plots/efficiency1.pdf");
  }else if(particle == 1){
    TCanvas c1;
    efficiency1->Draw("AP");
    c1.SaveAs("~/work2/BinQGP/results/Bs/efficiency/plots/efficiency1.pdf");
  }
  
  TEfficiency* acceptance = new TEfficiency(*hist_tot_noweights, *hist_tot_gen);
  if(particle == 0){
    TCanvas c2;
    acceptance->Draw("AP");
    c2.SaveAs("~/work2/BinQGP/results/Bu/efficiency/plots/acceptance.pdf");
  }else if(particle == 1){
    TCanvas c2;
    acceptance->Draw("AP");
    c2.SaveAs("~/work2/BinQGP/results/Bs/efficiency/plots/acceptance.pdf");
  }

  if(particle == 0){
    TFile* f0 = new TFile("~/work2/BinQGP/results/Bu/efficiency/root_files/efficiency0.root" , "recreate");
    f0->cd();
    efficiency0->Write();
    f0->Write();
    f0->ls();
    f0->Close();
    
    TFile* f1 = new TFile("~/work2/BinQGP/results/Bu/efficiency/root_files/efficiency1.root" , "recreate");
    f1->cd();
    efficiency1->Write();
    f1->Write();
    f1->ls();
    f1->Close();

    TFile* f2 = new TFile("~/work2/BinQGP/results/Bu/efficiency/root_files/acceptance.root", "recreate");
    f2->cd();
    acceptance->Write();
    f2->Write();
    f2->ls();
    f2->Close();
    
  }else if(particle == 1){
    TFile* f0 = new TFile("~/work2/BinQGP/results/Bs/efficiency/root_files/efficiency0.root" , "recreate");
    f0->cd();
    efficiency0->Write();
    f0->Write();
    f0->ls();
    f0->Close();
     
    TFile* f1 = new TFile("~/work2/BinQGP/results/Bs/efficiency/root_files/efficiency1.root" , "recreate");
    f1->cd();
    efficiency1->Write();
    f1->Write();
    f1->ls();
    f1->Close();

    TFile* f2 = new TFile("~/work2/BinQGP/results/Bs/efficiency/root_files/acceptance.root", "recreate");
    f2->cd();
    acceptance->Write();
    f2->Write();
    f2->ls();
    f2->Close();
    
  }
  /*
  for(int i = 1; i < n_pt_bins + 1; i++)
    {
      cout << efficiency0->GetEfficiency(i) << endl;
    }

  cout << endl;
  
  for(int i = 1; i < n_pt_bins + 1; i++)
    {
      cout << efficiency1->GetEfficiency(i) << endl;
    }
  */

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
  
  TString input_file = particle ? "~/work2/BinQGP/results/Bs/mc_validation_plots/weights/weights.root" :"~/work2/BinQGP/results/Bu/mc_validation_plots/weights/weights.root";

  TFile* f_wei = new TFile(input_file, "read");

  TH1D* histo_variable = (TH1D*)f_wei->Get(TString("weights_"+variable));

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
