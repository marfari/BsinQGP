#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>

int main(){

  TFile* f_mc = new TFile("/lstore/cms/julia/corrected_samples/NewMCBPlus.root");
  //a
  TTree* t1 = (TTree*)f_mc->Get("Bfinder/ntKp");
  t1->AddFriend("skimanalysis/HltTree");
  t1->AddFriend("hiEvtAnalyzer/HiTree");

  TH1F* hist_data = new TH1F("By_total", "By", 100, -2.6, 2.6);
  TH1F* hist_cut = new TH1F("By_cut", "By", 100, -2.6, 2.6);
  
c  t1->Project("By_total", "By", "pthatweight");
  hist_data->SetTitle("Total By");
  
  t1->Project("By_cut", "By", "pthatweight*((pprimaryVertexFilter && phfCoincFilter2Th4 && pclusterCompatibilityFilter && Btrk1Pt>0.9 && Bpt>5.0 && (BsvpvDistance/BsvpvDisErr)>2.0 && Bchi2cl>0.05 && TMath::Abs(Btrk1Eta)<2.4 && TMath::Abs(By)<2.4 && TMath::Abs(PVz)<15 && Bmass>5 && Bmass<6 && TMath::Abs(Bmumumass-3.096900)<0.15 && Bmu1SoftMuID && Bmu2SoftMuID && ((TMath::Abs(Bmu1eta)<1.2 && Bmu1pt>3.5) || (TMath::Abs(Bmu1eta)>1.2 && TMath::Abs(Bmu1eta)<2.1 && Bmu1pt>5.47-1.89*TMath::Abs(Bmu1eta)) || (TMath::Abs(Bmu1eta)>2.1 && TMath::Abs(Bmu1eta)<2.4 && Bmu1pt>1.5)) && ((TMath::Abs(Bmu2eta)<1.2 && Bmu2pt>3.5) || (TMath::Abs(Bmu2eta)>1.2 && TMath::Abs(Bmu2eta)<2.1 && Bmu2pt>5.47-1.89*TMath::Abs(Bmu2eta)) || (TMath::Abs(Bmu2eta)>2.1 && TMath::Abs(Bmu2eta)<2.4 && Bmu2pt>1.5)) && Bmu1isTriggered && Bmu2isTriggered && (Btrk1PixelHit+Btrk1StripHit)>=11 && (Btrk1Chi2ndf/(Btrk1nStripLayer+Btrk1nPixelLayer))<0.18 && TMath::Abs(Btrk1PtErr/Btrk1Pt)<0.1)&&((Bpt>5 && Bpt<7 && (BsvpvDistance/BsvpvDisErr)>16.457 && cos(Bdtheta)>0.987 && TMath::Abs(Btrk1Dxy1/Btrk1DxyError1)>0.005 && Btrk1Pt>1.092 && Bchi2cl>0.065) || (Bpt>7 && Bpt<10 && (BsvpvDistance/BsvpvDisErr)>12.714 && cos(Bdtheta)>0.947 && TMath::Abs(Btrk1Dxy1/Btrk1DxyError1)>2.928 && Btrk1Pt>0.838 && Bchi2cl>0.053) || (Bpt>10 && Bpt<15 && (BsvpvDistance/BsvpvDisErr)>9.086 && cos(Bdtheta)>0.994 && TMath::Abs(Btrk1Dxy1/Btrk1DxyError1)>1.540 && Btrk1Pt>1.262 && Bchi2cl>0.055) || (Bpt>15 && Bpt<20 && (BsvpvDistance/BsvpvDisErr)>7.587 && cos(Bdtheta)>0.757 && TMath::Abs(Btrk1Dxy1/Btrk1DxyError1)>0.000 && Btrk1Pt>1.813 && Bchi2cl>0.056) || (Bpt>20 && Bpt<30 && (BsvpvDistance/BsvpvDisErr)>4.004 && cos(Bdtheta)>0.996 && TMath::Abs(Btrk1Dxy1/Btrk1DxyError1)>0.000 && Btrk1Pt>1.822 && Bchi2cl>0.050) || (Bpt>30 && Bpt<50 && (BsvpvDistance/BsvpvDisErr)>2.000 && cos(Bdtheta)>0.998 && TMath::Abs(Btrk1Dxy1/Btrk1DxyError1)>0.000 && Btrk1Pt>2.046 && Bchi2cl>0.050) || (Bpt>50 && Bpt<100 && (BsvpvDistance/BsvpvDisErr)>4.084 && cos(Bdtheta)>-0.112 && TMath::Abs(Btrk1Dxy1/Btrk1DxyError1)>0.000 && Btrk1Pt>1.645 && Bchi2cl>0.050)))");
  hist_cut->SetTitle("By after cuts");

}
