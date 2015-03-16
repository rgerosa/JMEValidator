#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "THStack.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TLegend.h"
#include "TTree.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TColor.h"
#include "TEllipse.h"

#include "JMEAnalysis/JMEValidator/interface/validator_Ntuple.h"
#include "JMEAnalysis/JMEValidator/interface/puppi_Ntuple.h"
#include "JMEAnalysis/JMEValidator/macros/circle_intersection.C"
#include "/home/aperloff/Scripts/tdrstyle_mod14.C"

#include <iostream>
#include <iomanip>
#include <vector>
#include <assert.h>
#include <algorithm>

using std::cout;
using std::endl;
using std::vector;
using std::max;

class JetCircle {
public:
   JetCircle(double tx, double ty, double tr, double tpt);

   void findAngles(JetCircle wholeCircle);
   double getX() const {return x;}
   double getY() const {return y;}
   double getRadius() const {return radius;}
   double getPt() const {return pt;}
   double getStartAngle() const {return startAngle;}
   double getEndAngle() const {return endAngle;}
private:
   double x,y,radius,pt;
   double intersection_point_1_x, intersection_point_1_y;
   double intersection_point_2_x, intersection_point_2_y;
   double startAngle, endAngle;
};
JetCircle::JetCircle(double tx, double ty, double tr, double tpt) {
   x=tx;
   y=ty;
   radius=tr;
   pt = tpt;
   startAngle = 0;
   endAngle = 360;
}
void JetCircle::findAngles(JetCircle wholeCircle) {
  circle_circle_intersection(x, y, radius,
                             wholeCircle.getX(), wholeCircle.getY(), wholeCircle.getRadius(),
                             &intersection_point_1_x, &intersection_point_1_y,
                             &intersection_point_2_x, &intersection_point_2_y);
  circle_cirlce_intersection_angle(x,y,radius,intersection_point_1_x,intersection_point_1_y,&startAngle);
  circle_cirlce_intersection_angle(x,y,radius,intersection_point_2_x,intersection_point_2_y,&endAngle);
}
std::ostream& operator<<(std::ostream& out, const JetCircle& jet){
   return out << "x=" << jet.getX() << ", y=" << jet.getY()
              << ", r=" << jet.getRadius() << ", pt=" << jet.getPt() << "\n"
              << "startAngle=" << jet.getStartAngle() << ", endAngle=" << jet.getEndAngle();
}

void loadbar2(unsigned int x, unsigned int n, unsigned int w = 50) {
   if ( (x != n) && (x % max((unsigned int)1,(n/100)) != 0) ) return;

   float ratio  =  x/(float)n;
   int   c      =  ratio * w;

   cout << setw(3) << (int)(ratio*100) << "% [";
   for (int i=0; i<c; i++) cout << "=";
   for (unsigned int i=c; i<w; i++) cout << " ";
   cout << "] (" << x << "/" << n << ")\r" << flush;
}

void set_plot_style() {
   const Int_t NRGBs = 5;
   const Int_t NCont = 999;

   Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
   Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
   Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
   Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
   TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
   gStyle->SetNumberContours(NCont);
}

void particleDrawer(TString filename, int entry = 0, int PVAssoc = 2) {

   cout << "particleDrawer::Setting the TDR style ... ";
   setTDRStyle();
   cout << "DONE" << endl;

   cout << "particleDrawer::Drawing the default (TDR) frame ... " << endl;
   TH1D* frame = new TH1D();
   frame->GetXaxis()->SetLimits(-5,5);
   frame->GetXaxis()->SetTitle("#eta");
   frame->GetYaxis()->SetRangeUser(-TMath::Pi(),TMath::Pi());
   frame->GetYaxis()->SetTitle("#phi");
   TCanvas* c = tdrCanvas("particleBasedEvent",frame,4,0,true);
   c->GetPad(0)->SetLogz();
   cout << "\r\r\r\r\r\r" << flush;
   cout << setw(52) << " " << "DONE" << endl << endl;

   cout << "particleDrawer::Opening the input file (" << filename << " ) ... ";
   TFile* inFile = TFile::Open(filename,"READ");
   assert(inFile!=NULL);
   cout << "DONE" << endl;

   cout << "particleDrawer::Getting the input trees ... ";
   TTree* puppiTree = (TTree*)inFile->Get("puppiReader/puppiTree");
   assert(puppiTree!=NULL);
   TTree* jetTree = (TTree*)inFile->Get("nt_AK4PFchs/t");
   assert(jetTree!=NULL);
   cout << "DONE" << endl;

   cout << "particleDrawer::Making the ntuples ... ";
   puppiNtuple*     pNtuple = new puppiNtuple(puppiTree);
   validatorNtuple* jNtuple = new validatorNtuple(jetTree);
   cout << "DONE" << endl;

   cout << "particleDrawer::Getting entry " << entry << " for puppiTree ... ";
   puppiTree->GetEntry(entry);
   cout << "DONE" << endl;

   cout << "particleDrawer::Filling the histograms ... ";
   TH2F* hPU = new TH2F("hPU","hPU",50,-5,5,60,-TMath::Pi(),TMath::Pi());
   TH2F* hHard = new TH2F("hHard","hHard",50,-5,5,60,-TMath::Pi(),TMath::Pi());

   for(unsigned int iparticle=0; iparticle<(*pNtuple->px).size(); iparticle++) {
      TLorentzVector tempVect((*pNtuple->px)[iparticle],(*pNtuple->py)[iparticle],
                              (*pNtuple->pz)[iparticle],(*pNtuple->e)[iparticle]);

      if((*pNtuple->fromPV)[iparticle]<PVAssoc)
         hPU->Fill(tempVect.Eta(),tempVect.Phi(),tempVect.Pt());
      else
         hHard->Fill(tempVect.Eta(),tempVect.Phi(),tempVect.Pt());
   }
   cout << "DONE" << endl;

   cout << "particleDrawer::Drawing the histograms ... ";
   //tdrDraw(hPU,"BOX",kFullSquare,kNone,kSolid,kGray,kNone,kNone);
   //tdrDraw(hHard,"colz");
   THStack* stack = new THStack("stack","stack");
   hPU->SetLineStyle(kSolid);
   hPU->SetLineColor(kGray);
   hPU->SetFillStyle(1001);
   hPU->SetFillColor(kNone);
   hPU->SetMarkerStyle(kFullSquare);
   hPU->SetMarkerColor(kNone);
   stack->Add(hHard,"colz");
   stack->Add(hPU,"BOX");
   tdrDraw(stack,"nostack");

   set_plot_style();
   c->Update();
   c->RedrawAxis();
   cout << "DONE" << endl;

   cout << "particleDrawer::Getting entry " << entry << " for jetTree ... ";
   jetTree->GetEntry(entry);
   cout << "DONE" << endl;

   cout << "particleDrawer::Drawing the jets ... " << endl;
   vector<JetCircle> jets;
   for(unsigned int ijet=0; ijet<jNtuple->nref; ijet++) {
      if((*jNtuple->jtpt)[ijet]<20) continue;
      double RJet = TMath::Sqrt((*jNtuple->jtarea)[ijet]/TMath::Pi());
      jets.push_back(JetCircle((*jNtuple->jteta)[ijet],(*jNtuple->jtphi)[ijet],RJet,(*jNtuple->jtpt)[ijet]));
   }
   for(unsigned int ijet=0; ijet<jets.size(); ijet++) {
      for(unsigned int jjet=ijet+1; jjet<jets.size(); jjet++) {
         if(check_overlap(jets[ijet].getX(),jets[ijet].getY(),jets[ijet].getRadius(),
                          jets[jjet].getX(),jets[jjet].getY(),jets[jjet].getRadius())) {
            cout << "Jet " << ijet << " overlaps with jet " << jjet << endl;
            if(jets[ijet].getPt()>jets[jjet].getPt()) {
               //find angle for jjet;
               jets[jjet].findAngles(jets[ijet]);               
            }
            else if(jets[ijet].getPt()<jets[jjet].getPt()) {
               //find angle for ijet
               jets[ijet].findAngles(jets[jjet]);
            }
            else {
               //must find angle for both jets
               //then must draw a straight line between the two intersection points
               jets[jjet].findAngles(jets[ijet]);
               jets[ijet].findAngles(jets[jjet]);
               //NEED TO COMPLETE THIS FUNCTION. CURRENTLY DOESNOT DRAW LINE BETWEEN THE JETS.
            }
            cout << "Jet " << ijet << ": \n" << jets[ijet] << endl;
            cout << "Jet " << jjet << ": \n" << jets[jjet] << endl;
         }
      }
   }
   for(unsigned int ijet=0; ijet<jets.size(); ijet++) {
      loadbar2(ijet+1, jets.size());
      TEllipse* cJet = new TEllipse(jets[ijet].getX(),jets[ijet].getY(),
                                    jets[ijet].getRadius(),jets[ijet].getRadius(),
                                    jets[ijet].getStartAngle(),jets[ijet].getEndAngle());
      cJet->SetFillStyle(0);
      cJet->SetFillColor(kNone);
      cJet->SetLineStyle(kSolid);
      cJet->SetLineColor(kRed);
      cJet->SetLineWidth(3);
      cJet->Draw("only same");
   }
   //cout << "\r\r\r\r" << flush;
   //cout << setw(37) << " " << "DONE" << endl << endl;

   cout << "particleDrawer::Saving the canvas ... ";
   c->SaveAs("particleMap.png");
   c->SaveAs("particleMap.pdf");
   c->SaveAs("particleMap.C");
   cout << "DONE" << endl;
}
