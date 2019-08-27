
void counter(){
TFile *f_mc = new TFile("/eos/user/j/jusilva/2018_data/newMCsamples/PthatMergedreweighted.root");

TTree* t1 = (TTree*)f_mc->Get("Bfinder/ntphi");
t1->AddFriend("hiEvtAnalyzer/HiTree");

TFile *f_data = new TFile("/eos/user/j/jusilva/2018_data/crab_Bfinder_20181220_HIDoubleMuon_HIRun2018A_PromptReco_v1v2_1031_NoJSON_skimhltBsize_ntphi.root");

TTree* t2 = (TTree*)f_data->Get("Bfinder/ntphi");
t2->AddFriend("hiEvtAnalyzer/HiTree");


//checking number of events
//gStyle->SetOptStat("neirmou");
//int nevts = t1->Draw("EvtNo", "(Btrk1Pt > 1 && Btrk2Pt > 1 && Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2 && Bpt > 5 && abs(Btrk1Eta-0.0) < 2.4  && abs(Btrk2Eta-0.0) < 2.4 && Bgen==23333 && Bpt>7 && Bpt<15)");
//int nevts_bin1 = t1->Draw("EvtNo", "(pthatweight/100)*(Btrk1Pt > 1 && Btrk2Pt > 1 && Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2 && Bpt > 5 && abs(Btrk1Eta-0.0) < 2.4  && abs(Btrk2Eta-0.0) < 2.4 && Bgen==23333 && Bpt>7 && Bpt<15)");
//int nevts_bin2 = t1->Draw("EvtNo", "(pthatweight/100)*(Btrk1Pt > 1 && Btrk2Pt > 1 && Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2 && Bpt > 5 && abs(Btrk1Eta-0.0) < 2.4  && abs(Btrk2Eta-0.0) < 2.4 && Bgen==23333 && Bpt>15 && Bpt<30)");
//int nevts_bin3 = t1->Draw("EvtNo", "(pthatweight/100)*(Btrk1Pt > 1 && Btrk2Pt > 1 && Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2 && Bpt > 5 && abs(Btrk1Eta-0.0) < 2.4  && abs(Btrk2Eta-0.0) < 2.4 && Bgen==23333 && Bpt>30 && Bpt<50)");

//drawing trk eta distributions in the same plot 
TH1F* h1 = new TH1F("BtrkEta1", "BtrkEta", 100, -2.4, 2.4);
TH1F* h2 = new TH1F("BtrkEta2", "BtrkEta", 100, -2.4, 2.4);
TH1F* h3 = new TH1F("BtrkEta3", "BtrkEta", 100, -2.4, 2.4);
TH1F* h7 = new TH1F("BtrkEta7", "BtrkEta", 100, -2.4, 2.4);


t1->Project("BtrkEta7","Btrk1Eta", "(pthatweight/100)*(Btrk1Pt > 1 && Btrk2Pt > 1 && Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2 && Bpt > 5 && abs(Btrk1Eta-0.0) < 2.4  && abs(Btrk2Eta-0.0) < 2.4 && Bgen==23333 && Bpt>5 && Bpt<7)");
h7->SetTitle("BtrkEta");
h7->SetMarkerColor(8);
h7->SetLineColor(8);
h7->Scale(1/(h7->Integral()));
h7->SetStats(0);
t1->Project("BtrkEta1","Btrk1Eta", "(pthatweight/100)*(Btrk1Pt > 1 && Btrk2Pt > 1 && Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2 && Bpt > 5 && abs(Btrk1Eta-0.0) < 2.4  && abs(Btrk2Eta-0.0) < 2.4 && Bgen==23333 && Bpt>7 && Bpt<15)");
h1->SetTitle("BtrkEta");
h1->SetMarkerColor(kBlue);
h1->SetLineColor(kBlue);
h1->Scale(1/(h1->Integral()));
h1->SetStats(0);
t1->Project("BtrkEta2","Btrk1Eta", "(pthatweight/100)*(Btrk1Pt > 1 && Btrk2Pt > 1 && Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2 && Bpt > 5 && abs(Btrk1Eta-0.0) < 2.4  && abs(Btrk2Eta-0.0) < 2.4 && Bgen==23333 && Bpt<30 && Bpt>15)");
h2->SetNameTitle("Btrk1Eta", "Btrk1Eta");
h2->SetMarkerColor(kBlack);
h2->SetLineColor(kBlack);
h2->Scale(1/(h2->Integral()));
h2->SetStats(0);
t1->Project("BtrkEta3","Btrk1Eta", "(pthatweight/100)*(Btrk1Pt > 1 && Btrk2Pt > 1 && Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2 && Bpt > 5 && abs(Btrk1Eta-0.0) < 2.4  && abs(Btrk2Eta-0.0) < 2.4 && Bgen==23333 && Bpt<50 && Bpt>30)");
h3->SetTitle("BtrkEta");
h3->SetMarkerColor(kRed);
h3->SetLineColor(kRed);
h3->Scale(1/(h3->Integral()));
h3->SetStats(0);
h2->GetYaxis()->SetRangeUser(0.,0.05);
h2->SetXTitle("#eta");
TCanvas c1;
h2->Draw("hist");
h3->Draw("hist same");
h1->Draw("hist same");
h7->Draw("hist same");
TLegend *leg1 = new TLegend (0.7, 0.7, 0.9, 0.9);
leg1->AddEntry("BtrkEta7", "5 < p_{T} < 7 GeV ", "l");
leg1->AddEntry("BtrkEta1", "7 < p_{T} < 15 GeV ", "l");
leg1->AddEntry("BtrkEta2", "15 < p_{T} < 30 GeV ", "l");
leg1->AddEntry("BtrkEta3", "30 < p_{T} < 50 GeV ", "l");
leg1->Draw("same");
c1.SaveAs("comparison_trk1eta.png");

TH1F* h4 = new TH1F("BtrkEta4", "BtrkEta", 100, -2.4, 2.4);
TH1F* h5 = new TH1F("BtrkEta5", "BtrkEta", 100, -2.4, 2.4);
TH1F* h6 = new TH1F("BtrkEta6", "BtrkEta", 100, -2.4, 2.4);
TH1F* h8 = new TH1F("BtrkEta8", "BtrkEta", 100, -2.4, 2.4);

t1->Project("BtrkEta8","Btrk2Eta", "(pthatweight/100)*(Btrk1Pt > 1 && Btrk2Pt > 1 && Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2 && Bpt > 5 && abs(Btrk1Eta-0.0) < 2.4  && abs(Btrk2Eta-0.0) < 2.4 && Bgen==23333 && Bpt>5 && Bpt<7)");
h8->SetTitle("BtrkEta");
h8->SetMarkerColor(8);
h8->SetLineColor(8);
h8->Scale(1/(h8->Integral()));
h8->SetStats(0);
t1->Project("BtrkEta4","Btrk2Eta", "(pthatweight/100)*(Btrk1Pt > 1 && Btrk2Pt > 1 && Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2 && Bpt > 5 && abs(Btrk1Eta-0.0) < 2.4  && abs(Btrk2Eta-0.0) < 2.4 && Bgen==23333 && Bpt>7 && Bpt<15)");
h4->SetTitle("BtrkEta");
h4->SetMarkerColor(kBlue);
h4->SetLineColor(kBlue);
h4->Scale(1/(h4->Integral()));
h4->SetStats(0);
t1->Project("BtrkEta5","Btrk2Eta", "(pthatweight/100)*(Btrk1Pt > 1 && Btrk2Pt > 1 && Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2 && Bpt > 5 && abs(Btrk1Eta-0.0) < 2.4  && abs(Btrk2Eta-0.0) < 2.4 && Bgen==23333 && Bpt<30 && Bpt>15)");
h5->SetNameTitle("Btrk1Eta", "Btrk2Eta");
h5->SetMarkerColor(kBlack);
h5->SetLineColor(kBlack);
h5->Scale(1/(h5->Integral()));
h5->SetStats(0);
t1->Project("BtrkEta6","Btrk2Eta", "(pthatweight/100)*(Btrk1Pt > 1 && Btrk2Pt > 1 && Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2 && Bpt > 5 && abs(Btrk1Eta-0.0) < 2.4  && abs(Btrk2Eta-0.0) < 2.4 && Bgen==23333 && Bpt<50 && Bpt>30)");
h6->SetTitle("BtrkEta");
h6->SetMarkerColor(kRed);
h6->SetLineColor(kRed);
h6->Scale(1/(h6->Integral()));
h6->SetStats(0);
h5->GetYaxis()->SetRangeUser(0.,0.05);
h5->SetXTitle("#eta");
TCanvas c2;
h5->Draw("hist");
h6->Draw("hist same");
h4->Draw("hist same");
h8->Draw("hist same");
leg1->Draw("same");
c2.SaveAs("comparison_trk2eta.png");


t1->Draw("Btrk1Eta>>htemp1","(Btrk1Pt > 1 && Btrk2Pt > 1 && Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2 && Bpt > 5 && abs(Btrk1Eta-0.0) < 2.4  && abs(Btrk2Eta-0.0) < 2.4 && Bgen==23333)");
TH1F* htemp1 = (TH1F*)gDirectory->Get("htemp1");
htemp1->SetTitle("BtrkEta");
htemp1->SetMarkerColor(kBlue);
htemp1->SetLineColor(kBlue);

TCanvas c3;
htemp1->Draw("");
c3.SaveAs("total_trkEta.png");


/*t2->Draw("Btrk1Eta>>htemp1", "Btrk1Pt > 1 && Btrk2Pt > 1 && Bchi2cl > 0.05 && BsvpvDistance/BsvpvDisErr > 2 && Bpt > 5 && abs(Btrk1Eta-0.0) < 2.4  && abs(Btrk2Eta-0.0) < 2.4 && Bpt>7 && Bpt<15");

TH1F* htemp1_data = (TH1F*)gDirectory->Get("htemp1");
htemp1_data->SetTitle("BtrkEta");
htemp1_data->SetMarkerColor(kBlack);
htemp1_data->SetLineColor(kBlack);*/


//std::cout<< "evts: "<< nevts_bin1 << std::endl;
//std::cout<< "evts: "<< nevts_bin2 << std::endl;
//std::cout<< "evts: "<< nevts_bin3 << std::endl;
//std::cout<< "evts: "<< nevts << std::endl;
 
}
