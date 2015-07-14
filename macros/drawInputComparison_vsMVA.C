{

  gStyle->SetOptStat(00000);

  TFile* file1 = new TFile("eos/cms/store/user/rgerosa/PUPPETAnalysis/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt25ns_MCRUN2_74_V9_last.root","OPEN");
  file1->cd();
  TDirectory* dir1 = (TDirectory*) file1->Get("PUPPET");
  TTree* t1 = (TTree*) dir1->Get("t");

  TFile* file2 = new TFile("/afs/desy.de/user/r/rfriese/public/mvaTrainingSample-DYSpring15-50ns-v2/mvaTrainingSample.root","OPEN");
  file2->cd();
  TTree* t2 = (TTree*) gDirectory->Get("Flat");
  
  TH1F* NVertex_1 = new TH1F("NVertex_1","",50,0,50);
  TH1F* NVertex_2 = new TH1F("NVertex_2","",50,0,50);

  t1->Draw("NVertex >> NVertex_1","Boson_daughter==13","goff");
  t2->Draw("nPV >> NVertex_2","z_m > 70 && z_m < 120","goff");
  
  NVertex_1->Scale(1./NVertex_1->Integral());
  NVertex_2->Scale(1./NVertex_2->Integral());

  NVertex_1->SetLineColor(kBlue);
  NVertex_2->SetLineColor(kRed);
  NVertex_1->SetLineWidth(2);
  NVertex_2->SetLineWidth(2);
  NVertex_2->SetLineStyle(7);

  TCanvas* c1 = new TCanvas("c1","",500,500);
  c1->cd();
  NVertex_1->Draw("hist");
  NVertex_2->Draw("hist same");
  c1->SaveAs("NVertex.png","png");


  TH1F* ZPT_1 = new TH1F("ZPT_1","",50,0,250);
  TH1F* ZPT_2 = new TH1F("ZPT_2","",50,0,250);

  t1->Draw("Boson_Pt >> ZPT_1","Boson_daughter==13","goff");
  t2->Draw("z_pT >> ZPT_2","z_m > 70 && z_m < 120","goff");
  
  ZPT_1->Scale(1./ZPT_1->Integral());
  ZPT_2->Scale(1./ZPT_2->Integral());

  ZPT_1->SetLineColor(kBlue);
  ZPT_2->SetLineColor(kRed);
  ZPT_1->SetLineWidth(2);
  ZPT_2->SetLineWidth(2);
  ZPT_2->SetLineStyle(7);

  c1->cd();
  ZPT_1->Draw("hist");
  ZPT_2->Draw("hist same");
  c1->SaveAs("ZPT.png","png");


  TH1F* ZM_1 = new TH1F("ZM_1","",40,70,120);
  TH1F* ZM_2 = new TH1F("ZM_2","",40,70,120);

  t1->Draw("Boson_M >> ZM_1","Boson_daughter==13","goff");
  t2->Draw("z_m >> ZM_2","z_m > 70 && z_m < 120","goff");
  
  ZM_1->Scale(1./ZM_1->Integral());
  ZM_2->Scale(1./ZM_2->Integral());

  ZM_1->SetLineColor(kBlue);
  ZM_2->SetLineColor(kRed);
  ZM_1->SetLineWidth(2);
  ZM_2->SetLineWidth(2);
  ZM_2->SetLineStyle(7);

  c1->cd();
  ZM_1->Draw("hist");
  ZM_2->Draw("hist same");
  c1->SaveAs("ZM.png","png");

		        
  TH1F* NCleanedJets_1 = new TH1F("NCleanedJets_1","",10,0,10);
  TH1F* NCleanedJets_2 = new TH1F("NCleanedJets_2","",10,0,10);

  t1->Draw("NCleanedJets >> NCleanedJets_1","Boson_daughter==13","goff");
  t2->Draw("nJets >> NCleanedJets_2","z_m > 70 && z_m < 120","goff");
  
  NCleanedJets_1->Scale(1./NCleanedJets_1->Integral());
  NCleanedJets_2->Scale(1./NCleanedJets_2->Integral());

  NCleanedJets_1->SetLineColor(kBlue);
  NCleanedJets_2->SetLineColor(kRed);
  NCleanedJets_1->SetLineWidth(2);
  NCleanedJets_2->SetLineWidth(2);
  NCleanedJets_2->SetLineStyle(7);

  c1->cd();
  NCleanedJets_1->Draw();
  NCleanedJets_2->Draw("same");
  c1->SaveAs("NCleanedJets.png","png");
		        
  TH1F* LeadingJet_Pt_1 = new TH1F("LeadingJet_Pt_1","",100,0,100);
  TH1F* LeadingJet_Pt_2 = new TH1F("LeadingJet_Pt_2","",100,0,100);

  t1->Draw("LeadingJet_Pt >> LeadingJet_Pt_1","Boson_daughter==13","goff");
  t2->Draw("jet1_pT >> LeadingJet_Pt_2","z_m > 70 && z_m < 120","goff");
  
  LeadingJet_Pt_1->Scale(1./LeadingJet_Pt_1->Integral());
  LeadingJet_Pt_2->Scale(1./LeadingJet_Pt_2->Integral());

  LeadingJet_Pt_1->SetLineColor(kBlue);
  LeadingJet_Pt_2->SetLineColor(kRed);
  LeadingJet_Pt_1->SetLineWidth(2);
  LeadingJet_Pt_2->SetLineWidth(2);
  LeadingJet_Pt_2->SetLineStyle(7);

  c1->cd();
  LeadingJet_Pt_1->Draw();
  LeadingJet_Pt_2->Draw("same");
  c1->SaveAs("LeadingJet_Pt.png","png");

  TH1F* TrailingJet_Pt_1 = new TH1F("TrailingJet_Pt_1","",100,0,200);
  TH1F* TrailingJet_Pt_2 = new TH1F("TrailingJet_Pt_2","",100,0,200);

  t1->Draw("TrailingJet_Pt >> TrailingJet_Pt_1","Boson_daughter==13","goff");
  t2->Draw("jet2_pT >> TrailingJet_Pt_2","z_m > 70 && z_m < 120","goff");
  
  TrailingJet_Pt_1->Scale(1./TrailingJet_Pt_1->Integral());
  TrailingJet_Pt_2->Scale(1./TrailingJet_Pt_2->Integral());

  TrailingJet_Pt_1->SetLineColor(kBlue);
  TrailingJet_Pt_2->SetLineColor(kRed);
  TrailingJet_Pt_1->SetLineWidth(2);
  TrailingJet_Pt_2->SetLineWidth(2);
  TrailingJet_Pt_2->SetLineStyle(7);

  c1->cd();
  TrailingJet_Pt_1->Draw();
  TrailingJet_Pt_2->Draw("same");
  c1->SaveAs("TrailingJet_Pt.png","png");		      
    
  TH1F* LeadingJet_Eta_1 = new TH1F("LeadingJet_Eta_1","",30,-5,5);
  TH1F* LeadingJet_Eta_2 = new TH1F("LeadingJet_Eta_2","",30,-5,5);

  t1->Draw("LeadingJet_Eta >> LeadingJet_Eta_1","Boson_daughter==13","goff");
  t2->Draw("jet1_eta >> LeadingJet_Eta_2","z_m > 70 && z_m < 120","goff");
  
  LeadingJet_Eta_1->Scale(1./LeadingJet_Eta_1->Integral());
  LeadingJet_Eta_2->Scale(1./LeadingJet_Eta_2->Integral());

  LeadingJet_Eta_1->SetLineColor(kBlue);
  LeadingJet_Eta_2->SetLineColor(kRed);
  LeadingJet_Eta_1->SetLineWidth(2);
  LeadingJet_Eta_2->SetLineWidth(2);
  LeadingJet_Eta_2->SetLineStyle(7);

  c1->cd();
  LeadingJet_Eta_1->Draw();
  LeadingJet_Eta_2->Draw("same");
  c1->SaveAs("LeadingJet_Eta.png","png");

  TH1F* TrailingJet_Eta_1 = new TH1F("TrailingJet_Eta_1","",30,-5,5);
  TH1F* TrailingJet_Eta_2 = new TH1F("TrailingJet_Eta_2","",30,-5,5);

  t1->Draw("TrailingJet_Eta >> TrailingJet_Eta_1","Boson_daughter==13","goff");
  t2->Draw("jet2_eta >> TrailingJet_Eta_2","z_m > 70 && z_m < 120","goff");
  
  TrailingJet_Eta_1->Scale(1./TrailingJet_Eta_1->Integral());
  TrailingJet_Eta_2->Scale(1./TrailingJet_Eta_2->Integral());

  TrailingJet_Eta_1->SetLineColor(kBlue);
  TrailingJet_Eta_2->SetLineColor(kRed);
  TrailingJet_Eta_1->SetLineWidth(2);
  TrailingJet_Eta_2->SetLineWidth(2);
  TrailingJet_Eta_2->SetLineStyle(7);

  c1->cd();
  TrailingJet_Eta_1->Draw();
  TrailingJet_Eta_2->Draw("same");
  c1->SaveAs("TrailingJet_Eta.png","png");		      

    
  TH1F* recoilPFPuppiMet_Pt_1  = new TH1F("recoilPFPuppiMet_Pt_1","",100,0,200);
  TH1F* recoilPFPuppiMet_Pt_2  = new TH1F("recoilPFPuppiMet_Pt_2","",100,0,200);

  t1->Draw("recoilPFPuppiMet_Pt >> recoilPFPuppiMet_Pt_1","Boson_daughter==13","goff");
  t2->Draw("particleFlow_U >> recoilPFPuppiMet_Pt_2","z_m > 70 && z_m < 120","goff");
  
  recoilPFPuppiMet_Pt_1->Scale(1./recoilPFPuppiMet_Pt_1->Integral());
  recoilPFPuppiMet_Pt_2->Scale(1./recoilPFPuppiMet_Pt_2->Integral());

  recoilPFPuppiMet_Pt_1->SetLineColor(kBlue);
  recoilPFPuppiMet_Pt_2->SetLineColor(kRed);
  recoilPFPuppiMet_Pt_1->SetLineWidth(2);
  recoilPFPuppiMet_Pt_2->SetLineWidth(2);
  recoilPFPuppiMet_Pt_2->SetLineStyle(7);

  c1->cd();
  recoilPFPuppiMet_Pt_1->Draw();
  recoilPFPuppiMet_Pt_2->Draw("same");
  c1->SaveAs("recoilPFPuppiMet_Pt.png","png");		      

  TH1F* recoilPFPuppiMet_sumEt_1  = new TH1F("recoilPFPuppiMet_sumEt_1","",200,0,1000);
  TH1F* recoilPFPuppiMet_sumEt_2  = new TH1F("recoilPFPuppiMet_sumEt_2","",200,0,1000);

  t1->Draw("recoilPFPuppiMet_sumEt >> recoilPFPuppiMet_sumEt_1","Boson_daughter==13","goff");
  t2->Draw("particleFlow_SumET >> recoilPFPuppiMet_sumEt_2","z_m > 70 && z_m < 120","goff");
  
  recoilPFPuppiMet_sumEt_1->Scale(1./recoilPFPuppiMet_sumEt_1->Integral());
  recoilPFPuppiMet_sumEt_2->Scale(1./recoilPFPuppiMet_sumEt_2->Integral());

  recoilPFPuppiMet_sumEt_1->SetLineColor(kBlue);
  recoilPFPuppiMet_sumEt_2->SetLineColor(kRed);
  recoilPFPuppiMet_sumEt_1->SetLineWidth(2);
  recoilPFPuppiMet_sumEt_2->SetLineWidth(2);
  recoilPFPuppiMet_sumEt_2->SetLineStyle(7);

  c1->cd();
  recoilPFPuppiMet_sumEt_1->Draw();
  recoilPFPuppiMet_sumEt_2->Draw("same");
  c1->SaveAs("recoilPFPuppiMet_sumEt.png","png");		      
  
  
  TH1F* recoilPFPuppiMet_ChargedPU_Pt_1  = new TH1F("recoilPFPuppiMet_ChargedPU_Pt_1","",100,0,400);
  TH1F* recoilPFPuppiMet_ChargedPU_Pt_2  = new TH1F("recoilPFPuppiMet_ChargedPU_Pt_2","",100,0,400);

  t1->Draw("recoilPFPuppiMet_ChargedPU_Pt >> recoilPFPuppiMet_ChargedPU_Pt_1","Boson_daughter==13","goff");
  t2->Draw("pileUp_MET >> recoilPFPuppiMet_ChargedPU_Pt_2","z_m > 70 && z_m < 120","goff");
  
  recoilPFPuppiMet_ChargedPU_Pt_1->Scale(1./recoilPFPuppiMet_ChargedPU_Pt_1->Integral());
  recoilPFPuppiMet_ChargedPU_Pt_2->Scale(1./recoilPFPuppiMet_ChargedPU_Pt_2->Integral());

  recoilPFPuppiMet_ChargedPU_Pt_1->SetLineColor(kBlue);
  recoilPFPuppiMet_ChargedPU_Pt_2->SetLineColor(kRed);
  recoilPFPuppiMet_ChargedPU_Pt_1->SetLineWidth(2);
  recoilPFPuppiMet_ChargedPU_Pt_2->SetLineWidth(2);
  recoilPFPuppiMet_ChargedPU_Pt_2->SetLineStyle(7);

  c1->cd();
  recoilPFPuppiMet_ChargedPU_Pt_1->Draw();
  recoilPFPuppiMet_ChargedPU_Pt_2->Draw("same");
  c1->SaveAs("recoilPFPuppiMet_ChargedPU_Pt.png","png");		      


  TH1F* recoilPFPuppiMet_ChargedPU_sumEt_1  = new TH1F("recoilPFPuppiMet_ChargedPU_sumEt_1","",200,0,1000);
  TH1F* recoilPFPuppiMet_ChargedPU_sumEt_2  = new TH1F("recoilPFPuppiMet_ChargedPU_sumEt_2","",200,0,1000);

  t1->Draw("recoilPFPuppiMet_ChargedPU_sumEt >> recoilPFPuppiMet_ChargedPU_sumEt_1","Boson_daughter==13","goff");
  t2->Draw("pileUp_SumET >> recoilPFPuppiMet_ChargedPU_sumEt_2","z_m > 70 && z_m < 120","goff");
  
  recoilPFPuppiMet_ChargedPU_sumEt_1->Scale(1./recoilPFPuppiMet_ChargedPU_sumEt_1->Integral());
  recoilPFPuppiMet_ChargedPU_sumEt_2->Scale(1./recoilPFPuppiMet_ChargedPU_sumEt_2->Integral());

  recoilPFPuppiMet_ChargedPU_sumEt_1->SetLineColor(kBlue);
  recoilPFPuppiMet_ChargedPU_sumEt_2->SetLineColor(kRed);
  recoilPFPuppiMet_ChargedPU_sumEt_1->SetLineWidth(2);
  recoilPFPuppiMet_ChargedPU_sumEt_2->SetLineWidth(2);
  recoilPFPuppiMet_ChargedPU_sumEt_2->SetLineStyle(7);

  c1->cd();
  recoilPFPuppiMet_ChargedPU_sumEt_1->Draw();
  recoilPFPuppiMet_ChargedPU_sumEt_2->Draw("same");
  c1->SaveAs("recoilPFPuppiMet_ChargedPU_sumEt.png","png");		      

  TH1F* recoilPFPuppiMet_ChargedPV_Pt_1  = new TH1F("recoilPFPuppiMet_ChargedPV_Pt_1","",100,0,200);
  TH1F* recoilPFPuppiMet_ChargedPV_Pt_2  = new TH1F("recoilPFPuppiMet_ChargedPV_Pt_2","",100,0,200);

  t1->Draw("recoilPFPuppiMet_ChargedPV_Pt >> recoilPFPuppiMet_ChargedPV_Pt_1","Boson_daughter==13","goff");
  t2->Draw("track_U >> recoilPFPuppiMet_ChargedPV_Pt_2","z_m > 70 && z_m < 120","goff");
  
  recoilPFPuppiMet_ChargedPV_Pt_1->Scale(1./recoilPFPuppiMet_ChargedPV_Pt_1->Integral());
  recoilPFPuppiMet_ChargedPV_Pt_2->Scale(1./recoilPFPuppiMet_ChargedPV_Pt_2->Integral());

  recoilPFPuppiMet_ChargedPV_Pt_1->SetLineColor(kBlue);
  recoilPFPuppiMet_ChargedPV_Pt_2->SetLineColor(kRed);
  recoilPFPuppiMet_ChargedPV_Pt_1->SetLineWidth(2);
  recoilPFPuppiMet_ChargedPV_Pt_2->SetLineWidth(2);
  recoilPFPuppiMet_ChargedPV_Pt_2->SetLineStyle(7);

  c1->cd();
  recoilPFPuppiMet_ChargedPV_Pt_1->Draw();
  recoilPFPuppiMet_ChargedPV_Pt_2->Draw("same");
  c1->SaveAs("recoilPFPuppiMet_ChargedPV_Pt.png","png");		      

  TH1F* recoilPFPuppiMet_ChargedPV_sumEt_1  = new TH1F("recoilPFPuppiMet_ChargedPV_sumEt_1","",200,0,200);
  TH1F* recoilPFPuppiMet_ChargedPV_sumEt_2  = new TH1F("recoilPFPuppiMet_ChargedPV_sumEt_2","",200,0,200);

  t1->Draw("recoilPFPuppiMet_ChargedPV_sumEt >> recoilPFPuppiMet_ChargedPV_sumEt_1","Boson_daughter==13","goff");
  t2->Draw("track_SumET >> recoilPFPuppiMet_ChargedPV_sumEt_2","z_m > 70 && z_m < 120","goff");
  
  recoilPFPuppiMet_ChargedPV_sumEt_1->Scale(1./recoilPFPuppiMet_ChargedPV_sumEt_1->Integral());
  recoilPFPuppiMet_ChargedPV_sumEt_2->Scale(1./recoilPFPuppiMet_ChargedPV_sumEt_2->Integral());

  recoilPFPuppiMet_ChargedPV_sumEt_1->SetLineColor(kBlue);
  recoilPFPuppiMet_ChargedPV_sumEt_2->SetLineColor(kRed);
  recoilPFPuppiMet_ChargedPV_sumEt_1->SetLineWidth(2);
  recoilPFPuppiMet_ChargedPV_sumEt_2->SetLineWidth(2);
  recoilPFPuppiMet_ChargedPV_sumEt_2->SetLineStyle(7);

  c1->cd();
  recoilPFPuppiMet_ChargedPV_sumEt_1->Draw();
  recoilPFPuppiMet_ChargedPV_sumEt_2->Draw("same");
  c1->SaveAs("recoilPFPuppiMet_ChargedPV_sumEt.png","png");		      

    
  TH1F* recoilPFPuppiMet_noPU_Pt_1  = new TH1F("recoilPFPuppiMet_noPU_Pt_1","",100,0,500);
  TH1F* recoilPFPuppiMet_noPU_Pt_2  = new TH1F("recoilPFPuppiMet_noPU_Pt_2","",100,0,500);

  t1->Draw("(recoilPFPuppiMet_NeutralPV_Pt + recoilPFPuppiMet_ChargedPV_Pt)  >> recoilPFPuppiMet_noPU_Pt_1","Boson_daughter==13","goff");
  t2->Draw("noPileUp_U >> recoilPFPuppiMet_noPU_Pt_2","z_m > 70 && z_m < 120","goff");
  
  recoilPFPuppiMet_noPU_Pt_1->Scale(1./recoilPFPuppiMet_noPU_Pt_1->Integral());
  recoilPFPuppiMet_noPU_Pt_2->Scale(1./recoilPFPuppiMet_noPU_Pt_2->Integral());

  recoilPFPuppiMet_noPU_Pt_1->SetLineColor(kBlue);
  recoilPFPuppiMet_noPU_Pt_2->SetLineColor(kRed);
  recoilPFPuppiMet_noPU_Pt_1->SetLineWidth(2);
  recoilPFPuppiMet_noPU_Pt_2->SetLineWidth(2);
  recoilPFPuppiMet_noPU_Pt_2->SetLineStyle(7);

  c1->cd();
  recoilPFPuppiMet_noPU_Pt_1->Draw();
  recoilPFPuppiMet_noPU_Pt_2->Draw("same");
  c1->SaveAs("recoilPFPuppiMet_noPU_Pt.png","png");		      


  TH1F* recoilPFPuppiMet_noPU_sumEt_1  = new TH1F("recoilPFPuppiMet_noPU_sumEt_1","",200,0,200);
  TH1F* recoilPFPuppiMet_noPU_sumEt_2  = new TH1F("recoilPFPuppiMet_noPU_sumEt_2","",200,0,200);

  t1->Draw("(recoilPFPuppiMet_NeutralPV_sumEt + recoilPFPuppiMet_ChargedPV_sumEt) >> recoilPFPuppiMet_noPU_sumEt_1","Boson_daughter==13","goff");
  t2->Draw("pileUp_SumET >> recoilPFPuppiMet_noPU_sumEt_2","z_m > 70 && z_m < 120","goff");
  
  recoilPFPuppiMet_noPU_sumEt_1->Scale(1./recoilPFPuppiMet_noPU_sumEt_1->Integral());
  recoilPFPuppiMet_noPU_sumEt_2->Scale(1./recoilPFPuppiMet_noPU_sumEt_2->Integral());

  recoilPFPuppiMet_noPU_sumEt_1->SetLineColor(kBlue);
  recoilPFPuppiMet_noPU_sumEt_2->SetLineColor(kRed);
  recoilPFPuppiMet_noPU_sumEt_1->SetLineWidth(2);
  recoilPFPuppiMet_noPU_sumEt_2->SetLineWidth(2);
  recoilPFPuppiMet_noPU_sumEt_2->SetLineStyle(7);

  c1->cd();
  recoilPFPuppiMet_noPU_sumEt_1->Draw();
  recoilPFPuppiMet_noPU_sumEt_2->Draw("same");
  c1->SaveAs("recoilPFPuppiMet_noPU_sumEt.png","png");		      

}
