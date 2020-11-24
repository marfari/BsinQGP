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

#define particle 0 //0 = B+;   1 = Bs;

int counter0(){
  TFile* f_raw_yield = particle ?  new TFile("~/work2/BinQGP/results/Bs/Bpt/pT.root") : new TFile("~/work2/BinQGP/results/Bu/Bpt/pT.root");
  TGraphAsymmErrors* raw_yield = (TGraphAsymmErrors*)f_raw_yield->Get("Graph;1");

  TString input_f_mc_cuts = particle ? "/lstore/cms/nuno/ppdata2017/001120_nocut/TrkQualCut/BsMC.root" : "/lstore/cms/nuno/ppdata2017/001120_nocut/AnaCut/BPMC.root";
  TFile* f_mc_cuts = new TFile(input_f_mc_cuts);
  
  TString input_t = particle ? "ntphi" : "ntKp";
  TTree* t_cuts = (TTree*)f_mc_cuts->Get(input_t); 

  TString input_f_mc_nocuts = particle ? "/lstore/cms/nuno/ppdata2017/001120_nocut/NoCut/BsMC.root" : "/lstore/cms/nuno/ppdata2017/001120_nocut/NoCut/BPMC.root";
  TFile* f_mc_nocuts = new TFile(input_f_mc_nocuts);

  TTree* t_nocuts = (TTree*)f_mc_nocuts->Get("ntKp"); 
  TTree* t_gen = (TTree*)f_mc_nocuts->Get("ntGen");  

  double pt_bins[] = {7,10,15,20,50};
  const int n_pt_bins = raw_yield->GetN();

  TH1F* hist_tot_noweights = new TH1F("hist_tot_noweights", "hist_tot_noweights", n_pt_bins, pt_bins);
  TH1F* hist_passed_noweights = new TH1F("hist_passed_noweights", "hist_passed_noweights", n_pt_bins, pt_bins);

  TH1F* hist_tot_weights = new TH1F("hist_tot_weights", "hist_tot_weights", n_pt_bins, pt_bins);
  TH1F* hist_passed_weights = new TH1F("hist_passed_weights", "hist_passed_weights", n_pt_bins, pt_bins);

  TH1F* hist_tot_gen = new TH1F("hist_tot_gen", "hist_tot_gen", n_pt_bins, pt_bins);

  //NOCUTS 
  float bpt1;
  float pthat_nocut;
  float weight_nocut;

  t_nocuts->SetBranchAddress("Bpt", &bpt1);
  t_nocuts->SetBranchAddress("pthat", &pthat_nocuts);
  t_nocuts->SetBranchAddress("weight", &weight_nocuts);


  if(particle == 0){

  float bsvpvDisErr_2D; 

  t_nocuts->SetBranchAddress("BsvpvDisErr_2D", &bsvpvDisErr_2D);
  
  double weight = 1;

  for(int evt = 0; evt < t_nocuts->GetEntries(); evt++){

    t_nocuts->GetEntry(evt);
    
    weight = read_weights("BsvpvDisErr_2D", bsvpvDisErr_2D);
    weight *= pthat_nocuts*weight_nocuts;

    hist_tot_weights->Fill(bpt1, weight);
    hist_tot_noweights->Fill(bpt1);
   }
  
  }

  else if(particle == 1){

  float bdt_pt_5_10;
  float bdt_pt_10_15;
  float bdt_pt_15_20;
  float bdt_pt_20_50;
 
  t_nocuts->SetBranchAddress("BDT_pt_5_10", &bdt_pt_5_10);
  t_nocuts->SetBranchAddress("BDT_pt_10_15", &bdt_pt_10_15);
  t_nocuts->SetBranchAddress("BDT_pt_15_20", &bdt_pt_15_20);
  t_nocuts->SetBranchAddress("BDT_pt_20_50", &bdt_pt_20_50);

  double weight = 1;

  //Bin by bin analysis of the BDT for Bs
  TString variable;
  double bdt1_total = 0;


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
        weight *= pthat_nocuts*weight_nocuts;
      }      
      hist_tot_weights->Fill(bpt1, weight);
      hist_tot_noweights->Fill(bpt1);
    }
 }
  

    
  if(particle == 0){
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


  //CUTS
  float bpt2;
  float pthat_cuts;
  float weight_cuts;

  t_cuts->SetBranchAddress("Bpt", &bpt2);
  t_cuts->SetBranchAddress("pthat", &pthat_cuts);
  t_cuts->SetBranchAddress("weight", &weight_cuts);

  if(particle == 0){

  float bsvpvDisErr_2Dc;

  t_cuts->SetBranchAddress("BsvpvDisErr_2D", &bsvpvDisErr_2Dc);

  double weight2 = 1;

  for(int evt = 0; evt < t_cuts->GetEntries(); evt++){
    t_cuts->GetEntry(evt);

    weight2 = read_weights("BsvpvDisErr_2D", bsvpvDisErr_2Dc);
    weight2 *= pthat_cuts*weight_cuts;

    hist_passed_weights->Fill(bpt2, weight2);
    hist_passed_noweights->Fill(bpt2);
   }
  }


  else if(particle == 1){

  float bdt2_pt_5_10;
  float bdt2_pt_10_15;
  float bdt2_pt_15_20;
  float bdt2_pt_20_50;
   
  t_cuts->SetBranchAddress("BDT_pt_5_10", &bdt2_pt_5_10);
  t_cuts->SetBranchAddress("BDT_pt_10_15", &bdt2_pt_10_15);
  t_cuts->SetBranchAddress("BDT_pt_15_20", &bdt2_pt_15_20);
  t_cuts->SetBranchAddress("BDT_pt_20_50", &bdt2_pt_20_50);
 
 
  double weight2 = 1;
  double bdt2_total = 0;

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
        weight2 *= pthat_cuts*weight_cuts;
      }
      
      hist_passed_weights->Fill(bpt2, weight2);
      hist_passed_noweights->Fill(bpt2);
    }

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

  //GEN
  
  float bpt3;
  
  t_gen->SetBranchAddress("Gpt", &bpt3);
  
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

  //efficiency0 = events passing selection cuts+gen (SC)/ events passing gen (GC)(no weights)
  TEfficiency* efficiency0 = new TEfficiency(*hist_passed_noweights, *hist_tot_noweights);
  if(particle == 0){
    TCanvas c0;
    efficiency0->SetTitle("Nominal Efficiency #epsilon^{0};p_{T} (GeV);#epsilon^{0}");
    efficiency0->Draw("AP");
    c0.SaveAs("~/work2/BinQGP/results/Bu/efficiency/plots/efficiency0.pdf");
  }else if(particle == 1){
    TCanvas c0;
    efficiency0->SetTitle("Nominal Efficiency #epsilon^{0};p_{T} (GeV);#epsilon0");
    efficiency0->Draw("AP");
    c0.SaveAs("~/work2/BinQGP/results/Bs/efficiency/plots/efficiency0.pdf");
  }
  
  //efficiency1 with weights = events passing selection cuts / events passing gen (with weights)
  TEfficiency* efficiency1 = new TEfficiency(*hist_passed_weights, *hist_tot_weights);
  if(particle == 0){
    TCanvas c1;
    efficiency1->SetTitle("#epsilon^{1};p_{T} (GeV);#epsilon^{1}");
    efficiency1->Draw("AP");
    c1.SaveAs("~/work2/BinQGP/results/Bu/efficiency/plots/efficiency1.pdf");
  }else if(particle == 1){
    TCanvas c1;
    efficiency1->SetTitle("#epsilon^{1};p_{T} (GeV);#epsilon^{1}");
    efficiency1->Draw("AP");
    c1.SaveAs("~/work2/BinQGP/results/Bs/efficiency/plots/efficiency1.pdf");
  }
  
  //factor = events that pass (loose cuts) + gen (GC)/ gen (G)
  TEfficiency* acceptance = new TEfficiency(*hist_tot_noweights, *hist_tot_gen);
  if(particle == 0){
    TCanvas c2;
    acceptance->SetTitle("Factor that multiplied by eff0 gives eff*acc")
    acceptance->Draw("AP");
    c2.SaveAs("~/work2/BinQGP/results/Bu/efficiency/plots/acceptance.pdf");
  }else if(particle == 1){
    TCanvas c2;
    acceptance->SetTitle("Factor that multiplied by eff0 gives eff*acc")
    acceptance->Draw("AP");
    c2.SaveAs("~/work2/BinQGP/results/Bs/efficiency/plots/acceptance.pdf");
  }

  //Then efficiency * acceptance = efficiency0*acceptance = (SC/GC) * (GC*G) = SC/G = hist_passed_weights/hist_tot_gen

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
