#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TColor.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TPaletteAxis.h"
#include <iostream>

#include "CMS_lumi.h"

//root -l 'AllWidthObs_PlotPFDijetDMV_qq_1D.C("xsecUL_Asymptotic_qq_PFDijet2016.root","AllWidthObs_updated1DPlotDMV_qq")'
void AllWidthObs_PlotPFDijetDMV_qq_1D(string inputFileName,string outputDIR) {
  system(("mkdir -p "+outputDIR).c_str());
  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  // Set the color palette
  bool useNicksPalette = true;
  int ncontours        = 999;

  gStyle->SetNumberContours(ncontours);

  TGraph* grobsW30 = new TGraph();
  TGraph* grobsW25 = new TGraph();
  TGraph* grobsW20 = new TGraph();
  TGraph* grobsW15 = new TGraph();
  TGraph* grobsW10 = new TGraph();
  TGraph* grobsW05 = new TGraph();
  TGraph* grobsW01 = new TGraph();
  TGraph* grobsW0p001 = new TGraph();
  TGraph* grexp = new TGraph();

  TGraph* grexp_1sigma_up   = new TGraphErrors();
  TGraph* grexp_2sigma_up   = new TGraphErrors();
  TGraph* grexp_1sigma_dw   = new TGraphErrors();
  TGraph* grexp_2sigma_dw   = new TGraphErrors();

  int expcounter          = 0;
  int exp_up_counter_1s   = 0;
  int exp_down_counter_1s = 0;
  int exp_up_counter_2s   = 0;
  int exp_down_counter_2s = 0;
  int obscounter          = 0;

  int medMin = 1600;
  int medMax = 8000;

vector<int> medMassList;
string width[]={"0p3","0p25","0p2","0p15","0p1","0p05","0p01"};
//double width_[]={0.01,0.05,0.1,0.15,0.2,0.25,0.3};

  for(int w=0; w<7; w++) {
  string inputFileName1= "cards_qq_DMV_"+width[w]+"/"+inputFileName; 
  // This is where all the plots are made
  TFile *file = TFile::Open(inputFileName1.c_str(),"READ");
  TTree *tree = (TTree*)file->Get("xsecTree");
        double mass;
        double xsecULObs_PFDijet2016;
        double xsecULExpPlus2_PFDijet2016;
        double xsecULExpPlus_PFDijet2016;
        double xsecULExp_PFDijet2016;
        double xsecULExpMinus_PFDijet2016;
        double xsecULExpMinus2_PFDijet2016;
        tree->SetBranchAddress("mass", &mass);
        tree->SetBranchAddress("xsecULObs_PFDijet2016", &xsecULObs_PFDijet2016);
        tree->SetBranchAddress("xsecULExpPlus2_PFDijet2016", &xsecULExpPlus2_PFDijet2016);
        tree->SetBranchAddress("xsecULExpPlus_PFDijet2016", &xsecULExpPlus_PFDijet2016);
        tree->SetBranchAddress("xsecULExp_PFDijet2016", &xsecULExp_PFDijet2016);
        tree->SetBranchAddress("xsecULExpMinus_PFDijet2016", &xsecULExpMinus_PFDijet2016);
        tree->SetBranchAddress("xsecULExpMinus2_PFDijet2016", &xsecULExpMinus2_PFDijet2016);

  for (int i = 0; i < tree->GetEntries(); i++){
     tree->GetEntry(i);    
     if(mass<1600 or mass>8000) continue; 
      grexp->SetPoint(expcounter, double(mass), xsecULExp_PFDijet2016);
      expcounter++;
      cout<<"mass =" <<mass<<" xsecULExp ="<<xsecULExp_PFDijet2016<<" xsecULObs ="<<xsecULObs_PFDijet2016<<endl;
      // find max and min for frame
      medMassList.push_back(mass);

      
      if (width[w]=="0p3" and mass<=5000) { grobsW30->SetPoint(obscounter, double(mass),xsecULObs_PFDijet2016 ); obscounter++;}
      if (width[w]=="0p25" and mass<=5500) { grobsW25->SetPoint(obscounter, double(mass),xsecULObs_PFDijet2016 ); obscounter++;}
      if (width[w]=="0p2" and mass<=6000) { grobsW20->SetPoint(obscounter, double(mass),xsecULObs_PFDijet2016 ); obscounter++;}
      if (width[w]=="0p15" and mass<=6500) { grobsW15->SetPoint(obscounter, double(mass),xsecULObs_PFDijet2016 ); obscounter++;}
      if (width[w]=="0p1" and mass<=7000) { grobsW10->SetPoint(obscounter, double(mass),xsecULObs_PFDijet2016 ); obscounter++;}
      if (width[w]=="0p05" and mass<=7500) { grobsW05->SetPoint(obscounter, double(mass),xsecULObs_PFDijet2016 ); obscounter++;}
      if (width[w]=="0p01" and mass<=8000) { grobsW01->SetPoint(obscounter, double(mass),xsecULObs_PFDijet2016 ); obscounter++;}

    // 1 sigma dw
    grexp_1sigma_dw->SetPoint(exp_down_counter_1s, double(mass),xsecULExpMinus_PFDijet2016 );      
      exp_down_counter_1s++;
    // 1 sigma up
      grexp_1sigma_up->SetPoint(exp_up_counter_1s, double(mass),xsecULExpPlus_PFDijet2016 );      
      exp_up_counter_1s++;

    // 2 sigma dw
      grexp_2sigma_dw->SetPoint(exp_down_counter_2s, double(mass),xsecULExpMinus2_PFDijet2016 );      
      exp_down_counter_2s++;
    // 2 sigma up
      grexp_2sigma_up->SetPoint(exp_up_counter_2s, double(mass),xsecULExpPlus2_PFDijet2016 );      
      exp_up_counter_2s++;
  }
  tree->ResetBranchAddresses();
 }

 // Make 1 and 2 sigma brazilian bands
 TGraphAsymmErrors* graph_1sigma_band = new TGraphAsymmErrors();
 TGraphAsymmErrors* graph_2sigma_band = new TGraphAsymmErrors();

 if(exp_up_counter_1s == exp_down_counter_1s and exp_down_counter_1s == expcounter){
   for(int iPoint = 0; iPoint < exp_up_counter_1s; iPoint++){
     double x_central, y_central;
     grexp->GetPoint(iPoint,x_central,y_central);
     graph_1sigma_band->SetPoint(iPoint,x_central,y_central);
     double y_up, y_dw;
     grexp_1sigma_dw->GetPoint(iPoint,x_central,y_dw);
     grexp_1sigma_up->GetPoint(iPoint,x_central,y_up);
     float rangeDw = 0;
     float rangeUp = 0;
     if(iPoint == 0){
       rangeUp = (medMassList.at(iPoint+1)-medMassList.at(iPoint))/2;
     }
     else if(iPoint == exp_up_counter_1s-1){
       rangeDw = (medMassList.at(iPoint)-medMassList.at(iPoint-1))/2;
     }
     else{
       rangeUp = (medMassList.at(iPoint+1)-medMassList.at(iPoint))/2;
       rangeDw = (medMassList.at(iPoint)-medMassList.at(iPoint-1))/2;
     }

     double x_obs, y_obs;
     grobsW30->GetPoint(iPoint,x_obs,y_obs);
     graph_1sigma_band->SetPointError(iPoint,rangeDw,rangeUp,fabs(y_dw-y_central),fabs(y_up-y_central));      
     }
   }
 else {
   cerr<<"Number of expected limits value: mediat, 1-sigma up and 1-sigma down don't match --> skip "<<endl;
   return;
 }

 if(exp_up_counter_2s == exp_down_counter_2s and exp_down_counter_2s == expcounter){

   for(int iPoint = 0; iPoint < exp_up_counter_2s; iPoint++){
     double x_central, y_central;
     grexp->GetPoint(iPoint,x_central,y_central);
     graph_2sigma_band->SetPoint(iPoint,x_central,y_central);
     double y_up, y_dw;
     grexp_2sigma_dw->GetPoint(iPoint,x_central,y_dw);
     grexp_2sigma_up->GetPoint(iPoint,x_central,y_up);
     float rangeDw = 0;
     float rangeUp = 0;
     if(iPoint == 0){
       rangeUp = (medMassList.at(iPoint+1)-medMassList.at(iPoint))/2;
     }
     else if(iPoint == exp_up_counter_1s-1){
       rangeDw = (medMassList.at(iPoint)-medMassList.at(iPoint-1))/2;
     }
     else{
       rangeUp = (medMassList.at(iPoint+1)-medMassList.at(iPoint))/2;
       rangeDw = (medMassList.at(iPoint)-medMassList.at(iPoint-1))/2;
     }
     graph_2sigma_band->SetPointError(iPoint,rangeDw,rangeUp,fabs(y_dw-y_central),fabs(y_up-y_central));      
   }
 }
 else {
   cerr<<"Number of expected limits value: mediat, 2-sigma up and 2-sigma down don't match --> skip "<<endl;
   return;
 }

 cout<<"graph_2sigma_band "<<graph_2sigma_band->GetN()<<endl; 
 //////////// All the plotting and cosmetics
 TCanvas* canvas = new TCanvas("canvas", "canvas",600,600);
 TH1* frame = canvas->DrawFrame(medMin,TMath::MinElement(graph_2sigma_band->GetN(),graph_2sigma_band->GetY())*0.0001,
       			 medMax,TMath::MaxElement(graph_2sigma_band->GetN(),graph_2sigma_band->GetY())*4.5, "");
 //frame->GetYaxis()->CenterTitle();
 frame->GetXaxis()->SetTitle("Resonance Mass [GeV]");
 frame->GetYaxis()->SetTitle("#sigma #times #it{B} #times #it{A} [pb]");
 frame->GetXaxis()->SetTitleOffset(1.15);
 frame->GetYaxis()->SetTitleOffset(1.10);  
 frame->Draw();
 //CMS_lumi(canvas,"35.9");
 CMS_lumi(canvas,"36",false,true,false,0.,0);  

 int lineColor[] = {kRed,kGreen,kGreen+3,kMagenta+3,kCyan-4,kPink-9,kTeal+4,kMagenta-9,kSpring-9};
 grobsW30->SetLineColor(kRed);
 grobsW30->SetLineWidth(2);
 grobsW30->Draw("Lsame*");
 grobsW30->SetMarkerColor(kRed);
 grobsW30->SetMarkerStyle(20);

 grobsW25->SetLineColor(kGreen+3);
 grobsW25->SetLineWidth(2);
 grobsW25->Draw("Lsame*");
 grobsW25->SetMarkerColor(kGreen+3);
 grobsW25->SetMarkerStyle(22);

 grobsW20->SetLineColor(kGreen);
 grobsW20->SetLineWidth(2);
 grobsW20->Draw("Lsame*");
 grobsW20->SetMarkerColor(kGreen);
 grobsW20->SetMarkerStyle(21);

 grobsW15->SetLineColor(kMagenta+3);
 grobsW15->SetLineWidth(2);
 grobsW15->Draw("Lsame*");
 grobsW15->SetMarkerColor(kMagenta+3);
 grobsW15->SetMarkerStyle(24);

 grobsW10->SetLineColor(kBlue);
 grobsW10->SetLineWidth(2);
 grobsW10->Draw("Lsame*");
 grobsW10->SetMarkerColor(kBlue);
 grobsW10->SetMarkerStyle(23);

 grobsW05->SetLineColor(kPink-9);
 grobsW05->SetLineWidth(2);
 grobsW05->Draw("Lsame*");
 grobsW05->SetMarkerColor(kPink-9);
 grobsW05->SetMarkerStyle(25);

 grobsW01->SetLineColor(kBlack);
 grobsW01->SetLineWidth(2);
 grobsW01->Draw("Lsame*");
 grobsW01->SetMarkerColor(kBlack);
 grobsW01->SetMarkerStyle(26);

 grobsW0p001->SetLineColor(kMagenta-9);
 grobsW0p001->SetLineWidth(2);
 grobsW0p001->Draw("Lsame*");
 grobsW0p001->SetMarkerColor(kMagenta-9);
 grobsW0p001->SetMarkerStyle(27);

 TF1* line = new TF1 ("line","1",medMin,medMax);
 line->SetLineColor(kBlue);
 line->SetLineWidth(2);
 //  line->Draw("L same");

 TLegend *leg = new TLegend(0.60,0.55,0.99,0.86,NULL,"brNDC");
 leg->AddEntry(grobsW30,"#Gamma/M = 30%","PL");
 leg->AddEntry(grobsW25,"#Gamma/M = 25%","PL");
 leg->AddEntry(grobsW20,"#Gamma/M = 20%","PL");
 leg->AddEntry(grobsW15,"#Gamma/M = 15%","PL");
 leg->AddEntry(grobsW10,"#Gamma/M = 10%","PL");
 leg->AddEntry(grobsW05,"#Gamma/M = 5%","PL");
 leg->AddEntry(grobsW01,"#Gamma/M = 1%","PL");
 leg->SetFillColor(0);
 leg->SetFillStyle(0);
 leg->SetBorderSize(0);
 leg->Draw("SAME");
 
 TLatex * tex = new TLatex();
 tex->SetNDC();
 tex->SetTextFont(42);
 tex->SetLineWidth(2);
 tex->SetTextSize(0.042);
 tex->Draw();
 tex->DrawLatex(0.175,0.80,("quark-quark"));
 tex->DrawLatex(0.175,0.75,("Spin 1"));
 tex->DrawLatex(0.61,0.87,("#bf{Observed 95% CL}"));
 gPad->RedrawAxis();
 gPad->Modified(); 
 gPad->Update();
 
// canvas->SaveAs((outputDIR+"/DMV_qq_limit_AllWidth_Observed_UpdatedFinal_35p9fb_1D_"+to_string(int(1))+".pdf").c_str(),"pdf");
// canvas->SaveAs((outputDIR+"/DMV_qq_limit_AllWidth_Observed_UpdatedFinal_35p9fb_1D_"+to_string(int(1))+".png").c_str(),"pdf");

 canvas->SetLogy();
 frame->GetYaxis()->SetRangeUser(TMath::MinElement(graph_2sigma_band->GetN(),graph_2sigma_band->GetY())*0.1,
       			  TMath::MaxElement(graph_2sigma_band->GetN(),graph_2sigma_band->GetY())*200);
 canvas->SaveAs((outputDIR+"/DMV_qq_limit_AllWidth_Observed_UpdatedFinal_36fb_1D_"+to_string(int(1))+"_log.pdf").c_str(),"pdf");
 canvas->SaveAs((outputDIR+"/DMV_qq_limit_AllWidth_Observed_UpdatedFinal_36fb_1D_"+to_string(int(1))+"_log.png").c_str(),"pdf");
 canvas->SaveAs((outputDIR+"/DMV_qq_limit_AllWidth_Observed_UpdatedFinal_36fb_1D_"+to_string(int(1))+"_log.C").c_str(),"pdf");
}

