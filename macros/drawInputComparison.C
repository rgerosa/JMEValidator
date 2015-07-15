{

  gStyle->SetOptStat(00000);

  TFile* file1 = new TFile("eos/cms/store/user/rgerosa/PUPPETAnalysis/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt25ns_MCRUN2_74_V9_last.root","OPEN");
  file1->cd();
  TDirectory* dir1 = (TDirectory*) file1->Get("PUPPET");
  TTree* t1 = (TTree*) dir1->Get("t");

  TFile* file2 = new TFile("eos/cms/store/user/rgerosa/PUPPETAnalysis/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Asympt50ns_MCRUN2_74_V9A_last_fixed.root","OPEN");
  file2->cd();
  TDirectory* dir2 = (TDirectory*) file2->Get("PUPPET");
  TTree* t2 = (TTree*) dir2->Get("t");
  
  TH1F* NVertex_1 = new TH1F("NVertex_1","",50,0,50);
  TH1F* NVertex_2 = new TH1F("NVertex_2","",50,0,50);

  t1->Draw("NVertex >> NVertex_1","Boson_daughter==13","goff");
  t2->Draw("NVertex >> NVertex_2","Boson_daughter==13","goff");
  
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

		        
  TH1F* NCleanedJets_1 = new TH1F("NCleanedJets_1","",10,0,10);
  TH1F* NCleanedJets_2 = new TH1F("NCleanedJets_2","",10,0,10);

  t1->Draw("NCleanedJets >> NCleanedJets_1","Boson_daughter==13","goff");
  t2->Draw("NCleanedJets >> NCleanedJets_2","Boson_daughter==13","goff");
  
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
  t2->Draw("LeadingJet_Pt >> LeadingJet_Pt_2","Boson_daughter==13","goff");
  
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
  t2->Draw("TrailingJet_Pt >> TrailingJet_Pt_2","Boson_daughter==13","goff");
  
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
  
  TH1F* LeadingJet_Phi_1 = new TH1F("LeadingJet_Phi_1","",30,-3.14,3.14);
  TH1F* LeadingJet_Phi_2 = new TH1F("LeadingJet_Phi_2","",30,-3.14,3.14);

  t1->Draw("LeadingJet_Phi >> LeadingJet_Phi_1","Boson_daughter==13","goff");
  t2->Draw("LeadingJet_Phi >> LeadingJet_Phi_2","Boson_daughter==13","goff");
  
  LeadingJet_Phi_1->Scale(1./LeadingJet_Phi_1->Integral());
  LeadingJet_Phi_2->Scale(1./LeadingJet_Phi_2->Integral());

  LeadingJet_Phi_1->SetLineColor(kBlue);
  LeadingJet_Phi_2->SetLineColor(kRed);
  LeadingJet_Phi_1->SetLineWidth(2);
  LeadingJet_Phi_2->SetLineWidth(2);
  LeadingJet_Phi_2->SetLineStyle(7);

  c1->cd();
  LeadingJet_Phi_1->Draw();
  LeadingJet_Phi_2->Draw("same");
  c1->SaveAs("LeadingJet_Phi.png","png");		      


  TH1F* TrailingJet_Phi_1 = new TH1F("TrailingJet_Phi_1","",30,-3.14,3.14);
  TH1F* TrailingJet_Phi_2 = new TH1F("TrailingJet_Phi_2","",30,-3.14,3.14);

  t1->Draw("TrailingJet_Phi >> TrailingJet_Phi_1","Boson_daughter==13","goff");
  t2->Draw("TrailingJet_Phi >> TrailingJet_Phi_2","Boson_daughter==13","goff");
  
  TrailingJet_Phi_1->Scale(1./TrailingJet_Phi_1->Integral());
  TrailingJet_Phi_2->Scale(1./TrailingJet_Phi_2->Integral());

  TrailingJet_Phi_1->SetLineColor(kBlue);
  TrailingJet_Phi_2->SetLineColor(kRed);
  TrailingJet_Phi_1->SetLineWidth(2);
  TrailingJet_Phi_2->SetLineWidth(2);
  TrailingJet_Phi_2->SetLineStyle(7);

  c1->cd();
  TrailingJet_Phi_1->Draw();
  TrailingJet_Phi_2->Draw("same");
  c1->SaveAs("TrailingJet_Phi.png","png"); 
  
  TH1F* LeadingJet_Eta_1 = new TH1F("LeadingJet_Eta_1","",30,-5,5);
  TH1F* LeadingJet_Eta_2 = new TH1F("LeadingJet_Eta_2","",30,-5,5);

  t1->Draw("LeadingJet_Eta >> LeadingJet_Eta_1","Boson_daughter==13","goff");
  t2->Draw("LeadingJet_Eta >> LeadingJet_Eta_2","Boson_daughter==13","goff");
  
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
  t2->Draw("TrailingJet_Eta >> TrailingJet_Eta_2","Boson_daughter==13","goff");
  
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
  t2->Draw("recoilPFPuppiMet_Pt >> recoilPFPuppiMet_Pt_2","Boson_daughter==13","goff");
  
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


  TH1F* recoilPFPuppiMet_Phi_1  = new TH1F("recoilPFPuppiMet_Phi_1","",25,-3.14,3.14);
  TH1F* recoilPFPuppiMet_Phi_2  = new TH1F("recoilPFPuppiMet_Phi_2","",25,-3.14,3.14);

  t1->Draw("recoilPFPuppiMet_Phi >> recoilPFPuppiMet_Phi_1","Boson_daughter==13","goff");
  t2->Draw("recoilPFPuppiMet_Phi >> recoilPFPuppiMet_Phi_2","Boson_daughter==13","goff");
  
  recoilPFPuppiMet_Phi_1->Scale(1./recoilPFPuppiMet_Phi_1->Integral());
  recoilPFPuppiMet_Phi_2->Scale(1./recoilPFPuppiMet_Phi_2->Integral());

  recoilPFPuppiMet_Phi_1->SetLineColor(kBlue);
  recoilPFPuppiMet_Phi_2->SetLineColor(kRed);
  recoilPFPuppiMet_Phi_1->SetLineWidth(2);
  recoilPFPuppiMet_Phi_2->SetLineWidth(2);
  recoilPFPuppiMet_Phi_2->SetLineStyle(7);

  c1->cd();
  recoilPFPuppiMet_Phi_1->Draw();
  recoilPFPuppiMet_Phi_2->Draw("same");
  c1->SaveAs("recoilPFPuppiMet_Phi.png","png");		      


  TH1F* recoilPFPuppiMet_sumEt_1  = new TH1F("recoilPFPuppiMet_sumEt_1","",200,0,200);
  TH1F* recoilPFPuppiMet_sumEt_2  = new TH1F("recoilPFPuppiMet_sumEt_2","",200,0,200);

  t1->Draw("recoilPFPuppiMet_sumEt >> recoilPFPuppiMet_sumEt_1","Boson_daughter==13","goff");
  t2->Draw("recoilPFPuppiMet_sumEt >> recoilPFPuppiMet_sumEt_2","Boson_daughter==13","goff");
  
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
  t2->Draw("recoilPFPuppiMet_ChargedPU_Pt >> recoilPFPuppiMet_ChargedPU_Pt_2","Boson_daughter==13","goff");
  
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


  TH1F* recoilPFPuppiMet_ChargedPU_Phi_1  = new TH1F("recoilPFPuppiMet_ChargedPU_Phi_1","",25,-3.14,3.14);
  TH1F* recoilPFPuppiMet_ChargedPU_Phi_2  = new TH1F("recoilPFPuppiMet_ChargedPU_Phi_2","",25,-3.14,3.14);

  t1->Draw("recoilPFPuppiMet_ChargedPU_Phi >> recoilPFPuppiMet_ChargedPU_Phi_1","Boson_daughter==13","goff");
  t2->Draw("recoilPFPuppiMet_ChargedPU_Phi >> recoilPFPuppiMet_ChargedPU_Phi_2","Boson_daughter==13","goff");
  
  recoilPFPuppiMet_ChargedPU_Phi_1->Scale(1./recoilPFPuppiMet_ChargedPU_Phi_1->Integral());
  recoilPFPuppiMet_ChargedPU_Phi_2->Scale(1./recoilPFPuppiMet_ChargedPU_Phi_2->Integral());

  recoilPFPuppiMet_ChargedPU_Phi_1->SetLineColor(kBlue);
  recoilPFPuppiMet_ChargedPU_Phi_2->SetLineColor(kRed);
  recoilPFPuppiMet_ChargedPU_Phi_1->SetLineWidth(2);
  recoilPFPuppiMet_ChargedPU_Phi_2->SetLineWidth(2);
  recoilPFPuppiMet_ChargedPU_Phi_2->SetLineStyle(7);

  c1->cd();
  recoilPFPuppiMet_ChargedPU_Phi_1->Draw();
  recoilPFPuppiMet_ChargedPU_Phi_2->Draw("same");
  c1->SaveAs("recoilPFPuppiMet_ChargedPU_Phi.png","png");		      

  TH1F* recoilPFPuppiMet_ChargedPU_sumEt_1  = new TH1F("recoilPFPuppiMet_ChargedPU_sumEt_1","",200,0,1000);
  TH1F* recoilPFPuppiMet_ChargedPU_sumEt_2  = new TH1F("recoilPFPuppiMet_ChargedPU_sumEt_2","",200,0,1000);

  t1->Draw("recoilPFPuppiMet_ChargedPU_sumEt >> recoilPFPuppiMet_ChargedPU_sumEt_1","Boson_daughter==13","goff");
  t2->Draw("recoilPFPuppiMet_ChargedPU_sumEt >> recoilPFPuppiMet_ChargedPU_sumEt_2","Boson_daughter==13","goff");
  
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
  t2->Draw("recoilPFPuppiMet_ChargedPV_Pt >> recoilPFPuppiMet_ChargedPV_Pt_2","Boson_daughter==13","goff");
  
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

  TH1F* recoilPFPuppiMet_ChargedPV_Phi_1  = new TH1F("recoilPFPuppiMet_ChargedPV_Phi_1","",25,-3.14,3.14);
  TH1F* recoilPFPuppiMet_ChargedPV_Phi_2  = new TH1F("recoilPFPuppiMet_ChargedPV_Phi_2","",25,-3.14,3.14);

  t1->Draw("recoilPFPuppiMet_ChargedPV_Phi >> recoilPFPuppiMet_ChargedPV_Phi_1","Boson_daughter==13","goff");
  t2->Draw("recoilPFPuppiMet_ChargedPV_Phi >> recoilPFPuppiMet_ChargedPV_Phi_2","Boson_daughter==13","goff");
  
  recoilPFPuppiMet_ChargedPV_Phi_1->Scale(1./recoilPFPuppiMet_ChargedPV_Phi_1->Integral());
  recoilPFPuppiMet_ChargedPV_Phi_2->Scale(1./recoilPFPuppiMet_ChargedPV_Phi_2->Integral());

  recoilPFPuppiMet_ChargedPV_Phi_1->SetLineColor(kBlue);
  recoilPFPuppiMet_ChargedPV_Phi_2->SetLineColor(kRed);
  recoilPFPuppiMet_ChargedPV_Phi_1->SetLineWidth(2);
  recoilPFPuppiMet_ChargedPV_Phi_2->SetLineWidth(2);
  recoilPFPuppiMet_ChargedPV_Phi_2->SetLineStyle(7);

  c1->cd();
  recoilPFPuppiMet_ChargedPV_Phi_1->Draw();
  recoilPFPuppiMet_ChargedPV_Phi_2->Draw("same");
  c1->SaveAs("recoilPFPuppiMet_ChargedPV_Phi.png","png");		      

  TH1F* recoilPFPuppiMet_ChargedPV_sumEt_1  = new TH1F("recoilPFPuppiMet_ChargedPV_sumEt_1","",200,0,200);
  TH1F* recoilPFPuppiMet_ChargedPV_sumEt_2  = new TH1F("recoilPFPuppiMet_ChargedPV_sumEt_2","",200,0,200);

  t1->Draw("recoilPFPuppiMet_ChargedPV_sumEt >> recoilPFPuppiMet_ChargedPV_sumEt_1","Boson_daughter==13","goff");
  t2->Draw("recoilPFPuppiMet_ChargedPV_sumEt >> recoilPFPuppiMet_ChargedPV_sumEt_2","Boson_daughter==13","goff");
  
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

    
  TH1F* recoilPFPuppiMet_NeutralPV_Pt_1  = new TH1F("recoilPFPuppiMet_NeutralPV_Pt_1","",100,0,200);
  TH1F* recoilPFPuppiMet_NeutralPV_Pt_2  = new TH1F("recoilPFPuppiMet_NeutralPV_Pt_2","",100,0,200);

  t1->Draw("recoilPFPuppiMet_NeutralPV_Pt >> recoilPFPuppiMet_NeutralPV_Pt_1","Boson_daughter==13","goff");
  t2->Draw("recoilPFPuppiMet_NeutralPV_Pt >> recoilPFPuppiMet_NeutralPV_Pt_2","Boson_daughter==13","goff");
  
  recoilPFPuppiMet_NeutralPV_Pt_1->Scale(1./recoilPFPuppiMet_NeutralPV_Pt_1->Integral());
  recoilPFPuppiMet_NeutralPV_Pt_2->Scale(1./recoilPFPuppiMet_NeutralPV_Pt_2->Integral());

  recoilPFPuppiMet_NeutralPV_Pt_1->SetLineColor(kBlue);
  recoilPFPuppiMet_NeutralPV_Pt_2->SetLineColor(kRed);
  recoilPFPuppiMet_NeutralPV_Pt_1->SetLineWidth(2);
  recoilPFPuppiMet_NeutralPV_Pt_2->SetLineWidth(2);
  recoilPFPuppiMet_NeutralPV_Pt_2->SetLineStyle(7);

  c1->cd();
  recoilPFPuppiMet_NeutralPV_Pt_1->Draw();
  recoilPFPuppiMet_NeutralPV_Pt_2->Draw("same");
  c1->SaveAs("recoilPFPuppiMet_NeutralPV_Pt.png","png");		      


  TH1F* recoilPFPuppiMet_NeutralPV_Phi_1  = new TH1F("recoilPFPuppiMet_NeutralPV_Phi_1","",25,-3.14,3.14);
  TH1F* recoilPFPuppiMet_NeutralPV_Phi_2  = new TH1F("recoilPFPuppiMet_NeutralPV_Phi_2","",25,-3.14,3.14);


  t1->Draw("recoilPFPuppiMet_NeutralPV_Phi >> recoilPFPuppiMet_NeutralPV_Phi_1","Boson_daughter==13","goff");
  t2->Draw("recoilPFPuppiMet_NeutralPV_Phi >> recoilPFPuppiMet_NeutralPV_Phi_2","Boson_daughter==13","goff");
  
  recoilPFPuppiMet_NeutralPV_Phi_1->Scale(1./recoilPFPuppiMet_NeutralPV_Phi_1->Integral());
  recoilPFPuppiMet_NeutralPV_Phi_2->Scale(1./recoilPFPuppiMet_NeutralPV_Phi_2->Integral());

  recoilPFPuppiMet_NeutralPV_Phi_1->SetLineColor(kBlue);
  recoilPFPuppiMet_NeutralPV_Phi_2->SetLineColor(kRed);
  recoilPFPuppiMet_NeutralPV_Phi_1->SetLineWidth(2);
  recoilPFPuppiMet_NeutralPV_Phi_2->SetLineWidth(2);
  recoilPFPuppiMet_NeutralPV_Phi_2->SetLineStyle(7);

  c1->cd();
  recoilPFPuppiMet_NeutralPV_Phi_1->Draw();
  recoilPFPuppiMet_NeutralPV_Phi_2->Draw("same");
  c1->SaveAs("recoilPFPuppiMet_NeutralPV_Phi.png","png");		      


  TH1F* recoilPFPuppiMet_NeutralPV_sumEt_1  = new TH1F("recoilPFPuppiMet_NeutralPV_sumEt_1","",200,0,200);
  TH1F* recoilPFPuppiMet_NeutralPV_sumEt_2  = new TH1F("recoilPFPuppiMet_NeutralPV_sumEt_2","",200,0,200);

  t1->Draw("recoilPFPuppiMet_NeutralPV_sumEt >> recoilPFPuppiMet_NeutralPV_sumEt_1","Boson_daughter==13","goff");
  t2->Draw("recoilPFPuppiMet_NeutralPV_sumEt >> recoilPFPuppiMet_NeutralPV_sumEt_2","Boson_daughter==13","goff");
  
  recoilPFPuppiMet_NeutralPV_sumEt_1->Scale(1./recoilPFPuppiMet_NeutralPV_sumEt_1->Integral());
  recoilPFPuppiMet_NeutralPV_sumEt_2->Scale(1./recoilPFPuppiMet_NeutralPV_sumEt_2->Integral());

  recoilPFPuppiMet_NeutralPV_sumEt_1->SetLineColor(kBlue);
  recoilPFPuppiMet_NeutralPV_sumEt_2->SetLineColor(kRed);
  recoilPFPuppiMet_NeutralPV_sumEt_1->SetLineWidth(2);
  recoilPFPuppiMet_NeutralPV_sumEt_2->SetLineWidth(2);
  recoilPFPuppiMet_NeutralPV_sumEt_2->SetLineStyle(7);

  c1->cd();
  recoilPFPuppiMet_NeutralPV_sumEt_1->Draw();
  recoilPFPuppiMet_NeutralPV_sumEt_2->Draw("same");
  c1->SaveAs("recoilPFPuppiMet_NeutralPV_sumEt.png","png");		      


  TH1F* recoilPFPuppiMet_NeutralPU_Pt_1  = new TH1F("recoilPFPuppiMet_NeutralPU_Pt_1","",100,0,200);
  TH1F* recoilPFPuppiMet_NeutralPU_Pt_2  = new TH1F("recoilPFPuppiMet_NeutralPU_Pt_2","",100,0,200);

  t1->Draw("recoilPFPuppiMet_NeutralPU_Pt >> recoilPFPuppiMet_NeutralPU_Pt_1","Boson_daughter==13","goff");
  t2->Draw("recoilPFPuppiMet_NeutralPU_Pt >> recoilPFPuppiMet_NeutralPU_Pt_2","Boson_daughter==13","goff");
  
  recoilPFPuppiMet_NeutralPU_Pt_1->Scale(1./recoilPFPuppiMet_NeutralPU_Pt_1->Integral());
  recoilPFPuppiMet_NeutralPU_Pt_2->Scale(1./recoilPFPuppiMet_NeutralPU_Pt_2->Integral());

  recoilPFPuppiMet_NeutralPU_Pt_1->SetLineColor(kBlue);
  recoilPFPuppiMet_NeutralPU_Pt_2->SetLineColor(kRed);
  recoilPFPuppiMet_NeutralPU_Pt_1->SetLineWidth(2);
  recoilPFPuppiMet_NeutralPU_Pt_2->SetLineWidth(2);
  recoilPFPuppiMet_NeutralPU_Pt_2->SetLineStyle(7);

  c1->cd();
  recoilPFPuppiMet_NeutralPU_Pt_1->Draw();
  recoilPFPuppiMet_NeutralPU_Pt_2->Draw("same");
  c1->SaveAs("recoilPFPuppiMet_NeutralPU_Pt.png","png");		      


  TH1F* recoilPFPuppiMet_NeutralPU_Phi_1  = new TH1F("recoilPFPuppiMet_NeutralPU_Phi_1","",25,-3.14,3.14);
  TH1F* recoilPFPuppiMet_NeutralPU_Phi_2  = new TH1F("recoilPFPuppiMet_NeutralPU_Phi_2","",25,-3.14,3.14);

  t1->Draw("recoilPFPuppiMet_NeutralPU_Phi >> recoilPFPuppiMet_NeutralPU_Phi_1","Boson_daughter==13","goff");
  t2->Draw("recoilPFPuppiMet_NeutralPU_Phi >> recoilPFPuppiMet_NeutralPU_Phi_2","Boson_daughter==13","goff");
  
  recoilPFPuppiMet_NeutralPU_Phi_1->Scale(1./recoilPFPuppiMet_NeutralPU_Phi_1->Integral());
  recoilPFPuppiMet_NeutralPU_Phi_2->Scale(1./recoilPFPuppiMet_NeutralPU_Phi_2->Integral());

  recoilPFPuppiMet_NeutralPU_Phi_1->SetLineColor(kBlue);
  recoilPFPuppiMet_NeutralPU_Phi_2->SetLineColor(kRed);
  recoilPFPuppiMet_NeutralPU_Phi_1->SetLineWidth(2);
  recoilPFPuppiMet_NeutralPU_Phi_2->SetLineWidth(2);
  recoilPFPuppiMet_NeutralPU_Phi_2->SetLineStyle(7);

  c1->cd();
  recoilPFPuppiMet_NeutralPU_Phi_1->Draw();
  recoilPFPuppiMet_NeutralPU_Phi_2->Draw("same");
  c1->SaveAs("recoilPFPuppiMet_NeutralPU_Phi.png","png");		      


  TH1F* recoilPFPuppiMet_NeutralPU_sumEt_1  = new TH1F("recoilPFPuppiMet_NeutralPU_sumEt_1","",200,0,2000);
  TH1F* recoilPFPuppiMet_NeutralPU_sumEt_2  = new TH1F("recoilPFPuppiMet_NeutralPU_sumEt_2","",200,0,2000);

  t1->Draw("recoilPFPuppiMet_NeutralPU_sumEt >> recoilPFPuppiMet_NeutralPU_sumEt_1","Boson_daughter==13","goff");
  t2->Draw("recoilPFPuppiMet_NeutralPU_sumEt >> recoilPFPuppiMet_NeutralPU_sumEt_2","Boson_daughter==13","goff");
  
  recoilPFPuppiMet_NeutralPU_sumEt_1->Scale(1./recoilPFPuppiMet_NeutralPU_sumEt_1->Integral());
  recoilPFPuppiMet_NeutralPU_sumEt_2->Scale(1./recoilPFPuppiMet_NeutralPU_sumEt_2->Integral());

  recoilPFPuppiMet_NeutralPU_sumEt_1->SetLineColor(kBlue);
  recoilPFPuppiMet_NeutralPU_sumEt_2->SetLineColor(kRed);
  recoilPFPuppiMet_NeutralPU_sumEt_1->SetLineWidth(2);
  recoilPFPuppiMet_NeutralPU_sumEt_2->SetLineWidth(2);
  recoilPFPuppiMet_NeutralPU_sumEt_2->SetLineStyle(7);

  c1->cd();
  recoilPFPuppiMet_NeutralPU_sumEt_1->Draw();
  recoilPFPuppiMet_NeutralPU_sumEt_2->Draw("same");
  c1->SaveAs("recoilPFPuppiMet_NeutralPU_sumEt.png","png");		      


}
