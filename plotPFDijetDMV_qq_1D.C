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

//run 
//root -l 'plotPFDijetDMV_qq_1D.C("xsecUL_Asymptotic_qq_PFDijet2016.root","PlotDM_V",30)'

void plotPFDijetDMV_qq_1D(string inputFileName, string outputDIR,int w) {
 
  system(("mkdir -p "+outputDIR).c_str());
  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  // Set the color palette
  bool useNicksPalette = true;
  int ncontours        = 999;
  string inputFileName1="";
   string width[]={"0p3","0p25","0p2","0p15","0p1","0p05","0p01"};
   if(w==30) inputFileName1= "cards_qq_DMV_"+width[0]+"/"+inputFileName;
   if(w==25) inputFileName1= "cards_qq_DMV_"+width[1]+"/"+inputFileName;
   if(w==20) inputFileName1= "cards_qq_DMV_"+width[2]+"/"+inputFileName;
   if(w==15) inputFileName1= "cards_qq_DMV_"+width[3]+"/"+inputFileName;
   if(w==10) inputFileName1= "cards_qq_DMV_"+width[4]+"/"+inputFileName;
   if(w==5)  inputFileName1= "cards_qq_DMV_"+width[5]+"/"+inputFileName;
   if(w==1)  inputFileName1= "cards_qq_DMV_"+width[6]+"/"+inputFileName;


  TGraph* grobsNarrow = new TGraph();
  double narrowLimit[]={2.07e-01,2.99e-01,2.62e-01,1.61e-01,1.08e-01,7.56e-02,4.90e-02,2.86e-02,2.80e-02,3.05e-02,3.47e-02,3.88e-02,3.87e-02,3.53e-02,2.68e-02,1.36e-02,1.64e-02,
			1.78e-02,1.69e-02,1.48e-02,1.19e-02,1.01e-02,9.02e-03,7.72e-03,6.29e-03,5.17e-03,4.52e-03,4.61e-03,5.35e-03,5.65e-03,5.55e-03,5.26e-03,4.79e-03,3.88e-03,
			2.85e-03,2.14e-03,1.73e-03,1.45e-03,1.22e-03,1.01e-03,8.54e-04,7.88e-04,8.00e-04,8.09e-04,7.91e-04,7.45e-04,6.84e-04,6.26e-04,5.75e-04,5.21e-04,4.72e-04,
			4.30e-04,4.06e-04,4.00e-04,3.98e-04,3.94e-04,3.86e-04,3.74e-04,3.57e-04,3.33e-04,3.10e-04,2.84e-04,2.50e-04,2.20e-04,1.99e-04,1.97e-04};
  
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
  TGraph* grobs = new TGraph();
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
  for (int i = 0; i < tree->GetEntries(); i++){
    tree->GetEntry(i);    
     if(mass<1600 and mass>8000) continue; 
      grexp->SetPoint(expcounter, double(mass), xsecULExp_PFDijet2016);
      expcounter++;
      cout<<"mass =" <<mass<<" xsecULExp ="<<setprecision(9)<<xsecULExp_PFDijet2016<<" 		xsecULObs ="<<setprecision(9)<<xsecULObs_PFDijet2016<< endl;

      // find max and min for frame
/*      if(mass < medMin)
	medMin = mass;

      if(mass > medMax)
	medMax = mass;
*/
      medMassList.push_back(mass);

      grobs->SetPoint(obscounter, double(mass),xsecULObs_PFDijet2016 );
      grobsNarrow->SetPoint(i-6, double(mass),narrowLimit[i-6]);
      obscounter++;

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
      //}

      double x_obs, y_obs;
      grobs->GetPoint(iPoint,x_obs,y_obs);
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
  
  //////////// All the plotting and cosmetics
  TCanvas* canvas = new TCanvas("canvas", "canvas",600,600);
  TH1* frame = canvas->DrawFrame(medMin,TMath::MinElement(graph_2sigma_band->GetN(),graph_2sigma_band->GetY())*0.5,
				 medMax,TMath::MaxElement(graph_2sigma_band->GetN(),graph_2sigma_band->GetY())*1.5, "");
  //frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->SetTitle("Resonance Mass [GeV]");
  frame->GetYaxis()->SetTitle("#sigma #times #it{B} #times #it{A} (pb)");
  frame->GetXaxis()->SetTitleOffset(1.15);
  frame->GetYaxis()->SetTitleOffset(1.10);  
  frame->Draw();
  CMS_lumi(canvas,"36",false,true,false,0.37,0);

  graph_2sigma_band->SetFillColor(kOrange);
  graph_1sigma_band->SetFillColor(kGreen+1);
  graph_2sigma_band->SetLineColor(kOrange);
  graph_1sigma_band->SetLineColor(kGreen+1);
  
  graph_2sigma_band->Draw("3same");
  graph_1sigma_band->Draw("3same");

  grexp->SetLineColor(kBlack);
  grexp->SetLineStyle(2);
  grexp->SetLineWidth(2);
  grexp->Draw("Lsame");

  grobs->SetLineColor(kRed);
  grobs->SetLineWidth(2);
  grobs->Draw("Lsame");

  grobsNarrow->SetLineColor(kBlue);
  grobsNarrow->SetLineWidth(2);
  //grobsNarrow->Draw("Lsame");

  TF1* line = new TF1 ("line","1",medMin,medMax);
  line->SetLineColor(kBlue);
  line->SetLineWidth(2);
//  line->Draw("L same");

  TLegend *leg = new TLegend(0.53,0.57,0.92,0.84,NULL,"brNDC");
  leg->AddEntry(grobs,"Observed 95% CL","L");
  leg->AddEntry(grexp,"Median expected 95% CL","L");
  leg->AddEntry(graph_1sigma_band,"68% expected","F");
  leg->AddEntry(graph_2sigma_band,"95% expected","F");
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->Draw("SAME");
  
  TLatex * tex = new TLatex();
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetLineWidth(2);
  tex->SetTextSize(0.030);
  tex->Draw();
  tex->DrawLatex(0.175,0.80,("#bf{quark-quark}"));
  tex->DrawLatex(0.175,0.75,("#bf{Spin 1}"));

  if(w==30 )tex->DrawLatex(0.535,0.54,("#bf{ #Gamma/M = "+to_string(int(w))+"%}").c_str()); 
  if(w==25 )tex->DrawLatex(0.535,0.54,("#bf{ #Gamma/M = "+to_string(int(w))+"%}").c_str());
  if(w==20 )tex->DrawLatex(0.535,0.54,("#bf{ #Gamma/M = "+to_string(int(w))+"%}").c_str());
  if(w==15 )tex->DrawLatex(0.535,0.54,("#bf{ #Gamma/M = "+to_string(int(w))+"%}").c_str());
  if(w==10 )tex->DrawLatex(0.535,0.54,("#bf{ #Gamma/M = "+to_string(int(w))+"%}").c_str());
  if(w==5 )tex->DrawLatex(0.535,0.54,("#bf{ #Gamma/M = "+to_string(int(w))+"%}").c_str());
  if(w==1 )tex->DrawLatex(0.535,0.54,("#bf{ #Gamma/M = "+to_string(int(w))+"%}").c_str());
  tex->SetTextColor(kBlue);
  //tex->DrawLatex(0.175,0.60,("#bf{ Narrow Resonance Limits}"));
  gPad->RedrawAxis();
  gPad->Modified(); 
  gPad->Update();
  
  canvas->SaveAs((outputDIR+"/PFDMV_qq_limit_observed_UpdatedFinal_1D_W_"+to_string(int(w))+".pdf").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/PFDMV_qq_limit_observed_UpdatedFinal_1D_W_"+to_string(int(w))+".png").c_str(),"pdf");

  canvas->SetLogy();
  frame->GetYaxis()->SetRangeUser(TMath::MinElement(graph_2sigma_band->GetN(),graph_2sigma_band->GetY())*0.1,
				  TMath::MaxElement(graph_2sigma_band->GetN(),graph_2sigma_band->GetY())*200);
  canvas->SaveAs((outputDIR+"/PFDMV_qq_limit_observed_UpdatedFinal_1D_W_"+to_string(int(w))+"_log.pdf").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/PFDMV_qq_limit_observed_UpdatedFinal_1D_W_"+to_string(int(w))+"_log.png").c_str(),"pdf");

}

