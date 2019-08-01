#include <sstream>
#include <vector>
#include <cmath>

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
#include <TMath.h>

#include "/home/t3cms/ev19u032/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/interface/loop.h"
#include "/home/t3cms/ev19u032/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/interface/auxfunctions.h"
bool test=false;
bool optimized=true;
//pp:csyst=1, PbPb:csyst=2
void data_selection(TString fin1,TString data_selection_output_file, int csyst, int mc, int channel);


int main(int argc, char *argv[]){
  
  std::string collisionsystem="";
  int mc=0;
  int channel=4;
  
  for(int i=1 ; i<argc ; ++i){
    std::string argument = argv[i];
    std::stringstream convert;
    
    if(argument == "--colsyst"){
      convert << argv[++i];
      convert >> collisionsystem;
    }
    
    if(argument == "--mc"){
      convert << argv[++i];
      convert >> mc;
    }
    
    if(argument == "--channel"){
      convert << argv[++i];
      convert >> channel;
    }
    
  }
  
  
  TString input_file_pp = "/lstore/cms/julia/pp_files/ntB_EvtBase_20160420_BfinderData_pp_20160419_bPt0jpsiPt0tkPt0p5.root";
  TString output_file_pp= "/exper-sw/cmst3/cmssw/users/julia/CMSSW_7_5_8_patch5/src/UserCode/Bs_analysis/selected_data/selected_data_" + channel_to_ntuple_name(channel) + "_pp.root";
  TString input_file_PbPb = "/lstore/cms/julia/crab_Bfinder_20181220_HIDoubleMuon_HIRun2018A_PromptReco_v1v2_1031_NoJSON_skimhltBsize_ntKp.root";
  TString output_file_PbPb= "/home/t3cms/ev19u032/test/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/prefiltered_trees/selected_data_" + channel_to_ntuple_name(channel) + "_PbPb_2018.root";
  TString output_file_pp_mc= "/exper-sw/cmst3/cmssw/users/julia/CMSSW_7_5_8_patch5/src/UserCode/Bs_analysis/mc_for_central_code/selected_mc_ntphi_pp.root";
  TString input_file_pp_mc = "/lstore/cms/julia/pp_files/loop_Bs0_pthat5_bstojpsiphi_pp.root";
  TString output_file_PbPb_mc= "/home/t3cms/ev19u032/test/CMSSW_10_3_1_patch3/src/UserCode/BsinQGP/bin/prefiltered_trees/selected_mc_" + channel_to_ntuple_name(channel) + "_PbPb_2018_pthatweight.root";

  TString input_file_PbPb_mc = "/lstore/cms/julia/crab_Bfinder_20190221_Pythia8_BuToJpsiK_Bpt0p0_1032_NoJSON_pthatweight_hardcut_v2.root";
  
  if(collisionsystem=="pp" && mc==0){
    std::cout<<"Selecting data from pp collisions"<<std::endl;
    data_selection(input_file_pp,output_file_pp, 1, 0, channel);
  }
  if(collisionsystem=="PbPb" && mc==0){
    std::cout<<"Selecting data from PbPb collisions"<<std::endl;
    data_selection(input_file_PbPb,output_file_PbPb, 2, 0, channel);
  }
  if(collisionsystem=="pp" && mc==1){
    std::cout<<"Selecting mc from pp collisions"<<std::endl;
    data_selection(input_file_pp_mc,output_file_pp_mc, 1, 1, channel);
  }
  if(collisionsystem=="PbPb" && mc==1){
    std::cout<<"Selecting mc from PbPb collisions"<<std::endl;
    data_selection(input_file_PbPb_mc,output_file_PbPb_mc, 2, 1, channel);
  }
  
}

void data_selection(TString fin1, TString data_selection_output_file, int csyst, int mc, int channel){
  
  std::cout << "Reading Input TTree " << std::endl;
  
  
  TFile *f = new TFile(fin1);
  
  std::cout<< "Going to Bfinder directory " << std::endl;
  
  TTree* t1 = (TTree*)f->Get("Bfinder/ntKp");
  t1->AddFriend("skimanalysis/HltTree");
  t1->AddFriend("hiEvtAnalyzer/HiTree");
  
  
  readBranch(t1);
  
  std::cout << "Input TTree read " << std::endl;
  std::cout << "Creating new TTree " << std::endl;
  
  
  TFile *fout = new TFile(data_selection_output_file,"recreate");
  
  TTree *t4 = new TTree(channel_to_ntuple_name(channel),channel_to_ntuple_name(channel));
  
  Double_t mass;
  Double_t pt;
  Double_t y;
  Double_t trk1eta, trk2eta;
  Double_t chi2cl;
  Float_t svpvdistance, svpvdis_err;
  Float_t svpvdistance_2D, svpvdis_2D_err;
  Double_t alpha;
  Double_t trk1D0, trk2D0, trk1D0Err, trk2D0Err;
  Double_t trk1Dz, trk2Dz;
  Double_t trk1pt, trk2pt;
  Double_t mu1pt, mu2pt;
  Double_t mu1eta, mu2eta;
  Double_t d0, d0err;
  Double_t lxy;
  Float_t  Pthatweight;
  Int_t HiBin;
  std::cout << "New TTree Created " << std::endl;
  
  t4->Branch("Bmass", &mass);
  t4->Branch("Bpt", &pt);
  t4->Branch("By", &y);
  t4->Branch("Btrk1eta", &trk1eta);
  t4->Branch("Btrk2eta", &trk2eta);
  t4->Branch("Btrk1pt", &trk1pt);
  t4->Branch("Btrk2pt", &trk2pt);
  t4->Branch("Bmu1eta", &mu1eta);
  t4->Branch("Bmu2eta", &mu2eta);
  t4->Branch("Bmu1pt", &mu1pt);
  t4->Branch("Bmu2pt", &mu2pt);
  t4->Branch("Bchi2cl", &chi2cl);
  t4->Branch("BsvpvDistance_2D", &svpvdistance_2D);
  t4->Branch("BsvpvDistance_2D_Err", &svpvdis_2D_err);
  t4->Branch("BsvpvDistance", &svpvdistance);
  t4->Branch("BsvpvDistance_Err", &svpvdis_err);
  t4->Branch("Balpha", &alpha);
  t4->Branch("Btrk1D0", &trk1D0);
  t4->Branch("Btrk2D0", &trk2D0);
  t4->Branch("Btrk1D0Err", &trk1D0Err);
  t4->Branch("Btrk2D0Err", &trk2D0Err);
  t4->Branch("Btrk1Dz", &trk1Dz);
  t4->Branch("Btrk2Dz", &trk2Dz);
  t4->Branch("Bd0", &d0);
  t4->Branch("Bd0err", &d0err);
  t4->Branch("Blxy", &lxy);
  t4->Branch("Pthatweight", &Pthatweight);
  t4->Branch("HiBin", &HiBin);
  
  std::cout<<"Loop over events"<<std::endl;
  
  int count=0;
  for (int evt=0;evt<(t1->GetEntries());evt++) {
    
    t1->GetEntry(evt);
    
    for(int i=0; i<Bsize; i++){
      if(csyst==2){
        if(channel==1){
          //CUT FROM KIM
          if(Btrk1Pt[i]<=1.0) continue;
          if(Bpt[i]<5.0) continue;
          if((BsvpvDistance[i]/BsvpvDisErr[i])<=2.0) continue;
          if(Bchi2cl[i]<=0.05) continue;
          if(TMath::Abs(Btrk1Eta[i])>=2.4) continue;
          if(TMath::Abs(By[i])>=2.4) continue;
          if(TMath::Abs(PVz)>=15) continue;
          if(Bmass[i]<5) continue;
          if(Bmass[i]>6) continue;
          if(!Bmu1isGlobalMuon[i]) continue;
          if(!Bmu2isGlobalMuon[i]) continue;
          if(Bmu1dxyPV[i]>=0.3) continue;
          if(Bmu2dxyPV[i]>=0.3) continue;
          if(Bmu1dzPV[i]>=20) continue;
          if(Bmu2dzPV[i]>=20) continue;
          if(Bmu1InPixelLayer[i]<0) continue;
          if((Bmu1InPixelLayer[i]+Bmu1InStripLayer[i])<=5) continue;
          if(Bmu2InPixelLayer[i]<=0) continue;
          if((Bmu2InPixelLayer[i]+Bmu2InStripLayer[i])<=5) continue;
          if(((TMath::Abs(Bmu1eta[i])<1.2 && Bmu1pt[i]>3.5) || (TMath::Abs(Bmu1eta[i])>1.2 && TMath::Abs(Bmu1eta[i])<2.1 && Bmu1pt[i]>(5.77-1.8*TMath::Abs(Bmu1eta[i]))) || (TMath::Abs(Bmu1eta[i])>2.1 && TMath::Abs(Bmu1eta[i])<2.4 && Bmu1pt[i]>1.8))==0) continue;
          // if((abs(Bmu1eta[i])<1.2 || Bmu1pt[i]<=3.5) && (abs(Bmu1eta[i])>=1.2 || abs(Bmu1eta[i])<2.1 || Bmu1pt[i]<=(5.77-1.8*abs(Bmu1eta[i]))) && (abs(Bmu1eta[i])>=2.1 || abs(Bmu1eta[i])<2.4 || Bmu1pt[i]<=1.8)) continue;
          if(TMath::Abs(Bmumumass[i]-3.096900)>=0.15 ) continue;
          if(!Bmu1TMOneStationTight[i]) continue;
          if(!Bmu2TMOneStationTight[i]) continue;
          if(!Btrk1highPurity[i]) continue;
          if((Btrk1PixelHit[i]+Btrk1StripHit[i])<11) continue;
          if((Btrk1Chi2ndf[i]/(Btrk1nStripLayer[i]+Btrk1nPixelLayer[i]))>=0.18) continue;
          if(TMath::Abs(Btrk1PtErr[i]/Btrk1Pt[i])>=0.1) continue;
          
          if(mc==1){
            if(Bgen[i]!=23333) continue;
          }
        }
        
        if(optimized){
          if(Bpt[i]>5 && Bpt[i]<7){
            if((BsvpvDistance[i]/BsvpvDisErr[i])<10.497) continue;
            if(cos(Bdtheta[i])<0.922) continue;
            if(Btrk1Pt[i]<0.983) continue;
            if(Bchi2cl[i]<0.064) continue;
            if(TMath::Abs(Btrk1Eta[i])>2.380) continue;
          }
          
          if(Bpt[i]>7 && Bpt[i]<10){
            if((BsvpvDistance[i]/BsvpvDisErr[i])<15.208) continue;
            if(cos(Bdtheta[i])<0.966) continue;
            if(Btrk1Pt[i]<0.976) continue;
            if(Bchi2cl[i]<0.057) continue;
          }
          
          if(Bpt[i]>10 && Bpt[i]<15){
            if((BsvpvDistance[i]/BsvpvDisErr[i])<7.091) continue;
            if(cos(Bdtheta[i])<0.983) continue;
            if(Btrk1Pt[i]<1.540) continue;
            if(Bchi2cl[i]<0.052) continue;
          }
          
          if(Bpt[i]>15 && Bpt[i]<20){
            if((BsvpvDistance[i]/BsvpvDisErr[i])<5.884) continue;
            if(cos(Bdtheta[i])<0.616) continue;
            if(Btrk1Pt[i]<1.336) continue;
            if(Bchi2cl[i]<0.055) continue;
          }
          
          if(Bpt[i]>20 && Bpt[i]<30){
            if((BsvpvDistance[i]/BsvpvDisErr[i])<3.473) continue;
            if(cos(Bdtheta[i])<-0.535) continue;
            if(Btrk1Pt[i]<1.727) continue;
            if(Bchi2cl[i]<0.048) continue;
          }
          if(Bpt[i]>30 && Bpt[i]<50){
            if(cos(Bdtheta[i])<0.993) continue;
            if(Btrk1Pt[i]<1.577) continue;
            if(Bchi2cl[i]<0.050) continue;
            if(TMath::Abs(Btrk1Eta[i])>2.387) continue;
          }
          
          
        }
        
        
        if(channel==4){
          if(abs(Btktkmass[i]-PHI_MASS)>=0.2) continue;
          if(abs(Btrk2Eta[i])>=1.5) continue;//1.5 para Bs
          if(abs(Btrk1Eta[i])>=1.5) continue;//1.5 para Bs
          if(Btrk2Pt[i]<=1.5) continue;
          if(Btrk1Pt[i]<=1.5) continue;
          if(Bchi2cl[i]<=0.1) continue;//0.1
          if((Bd0[i]/Bd0Err[i])<=5.) continue;
          if(cos(Bdtheta[i])<=0.8) continue;
          if(Blxy[i]<=0.1) continue;
          if(Bpt[i]<=5) continue;
        }
        
        
      }
      // if(i==0) pvz=PVz;
      
      mass=Bmass[i];
      pt=Bpt[i];
      y=By[i];
      trk1eta=Btrk1Eta[i];
      trk2eta=Btrk2Eta[i];
      chi2cl=Bchi2cl[i];
      svpvdis_err=BsvpvDisErr[i];
      svpvdistance=BsvpvDistance[i];
      svpvdistance_2D=BsvpvDistance_2D[i];
      svpvdis_2D_err=BsvpvDisErr_2D[i];
      alpha=Balpha[i];
      trk1D0=Btrk1D0[i];
      trk2D0=Btrk2D0[i];
      trk1D0Err=Btrk1D0Err[i];
      trk2D0Err=Btrk2D0Err[i];
      trk1Dz=Btrk1Dz[i];
      trk2Dz=Btrk2Dz[i];
      trk1pt=Btrk1Pt[i];
      trk2pt=Btrk2Pt[i];
      mu1pt=Bmu1pt[i];
      mu2pt=Bmu2pt[i];
      mu1eta=Bmu1eta[i];
      mu2eta=Bmu2eta[i];
      d0= Bd0[i];
      d0err= Bd0Err[i];
      lxy= Blxy[i];
      HiBin=hiBin;
      if(mc==1) Pthatweight=pthatweight;
      if(mc==0) Pthatweight=0;
      t4->Fill();
      count++;
    }
  }
  std::cout<<"Selection Over"<<std::endl;
  
  t4->Write();
  f->Close();
  fout->Close();
  
}
