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
#include "iostream"
#include "math.h"
#include "CMS_lumi.h"

//run 
//root -l 'plotPFDijetDMV_qq_1D_limit_vs_width.C("xsecUL_Asymptotic_qq_PFDijet2016.root","newPlotDM_V",1500,0,0,0)'
//arg1 first line  x and y location (0,1,2,3,4,5,6,7,8,9,10) for cross section 
//arg2 second line x and y psition (0,1,2,3,4,5,6) for linmit
//arg3 0=obs limit, 1=expPlus2, 2=expPlus, 3=exp, 4=expMinus,5=expMinus2,

void plotPFDijetDMV_qq_1D_limit_vs_width(string inputFileName, string outputDIR,double zm, int arg1,int arg2,int arg3) {
  system(("mkdir -p "+outputDIR).c_str());
  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  // Set the color palette
  bool useNicksPalette = true;
  int ncontours        = 999;

  if (useNicksPalette) {    
    TColor::InitializeColors();
    Double_t stops[9] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000};
    Double_t red[9]   = { 243./255., 243./255., 240./255., 240./255., 241./255., 239./255., 186./255., 151./255., 129./255.};
    Double_t green[9] = {   0./255.,  46./255.,  99./255., 149./255., 194./255., 220./255., 183./255., 166./255., 147./255.};
    Double_t blue[9]  = {   6./255.,   8./255.,  36./255.,  91./255., 169./255., 235./255., 246./255., 240./255., 233./255.};
    TColor::CreateGradientColorTable(9, stops, red, green, blue, ncontours);
  }
  else gStyle->SetPalette(70);  
  gStyle->SetNumberContours(ncontours);
 
  TGraph* grXsec = new TGraph();
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

  int medMin = 0;
  double medMax = 0.6;
  vector<int> medMassList;
 
  // This is where all the plots are made
string width[]={"0p01","0p05","0p1","0p15","0p2","0p25","0p3"};
double width_[]={0.01,0.05,0.1,0.15,0.2,0.25,0.3};
float width_limit[7][7];
for (int i=0; i<7; i++){
for (int j=0; j<7; j++){
	width_limit[i][j]=0;
}
}

for(int w=0; w<7; w++) {
  string inputFileName1= "cards_qq_DMV_"+width[w]+"/"+inputFileName;
//  cout<<" file name = "<<inputFileName1<<endl;
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
      if (mass!=zm) continue;
      grexp->SetPoint(expcounter, double(width_[w]), xsecULExp_PFDijet2016);
      expcounter++;
      //cout<<"mass =" <<mass<<endl;
      cout<<" X = "<<width_[w]<<" Y xsecULObs ="<<xsecULObs_PFDijet2016<<endl;
      cout<<" X = "<<width_[w]<<" Y xsecULExpplus2 ="<<xsecULExpPlus2_PFDijet2016 << endl;
      cout<<" X = "<<width_[w]<<" Y xsecULExpPlus "<<xsecULExpPlus_PFDijet2016<<endl;
      cout<<" X = "<<width_[w]<<" Y xsecULExp ="<<xsecULExp_PFDijet2016<<endl;
      cout<<" X = "<<width_[w]<<" Y xsecULExpMinus ="<<xsecULExpMinus_PFDijet2016<<endl;
      cout<<" X = "<<width_[w]<<" Y xsecULExpMinus2 ="<<xsecULExpMinus2_PFDijet2016<<endl;
	width_limit[w][0]=xsecULObs_PFDijet2016;
        width_limit[w][1]=xsecULExpPlus2_PFDijet2016;
	width_limit[w][2]=xsecULExpPlus_PFDijet2016;
	width_limit[w][3]=xsecULExp_PFDijet2016;
	width_limit[w][4]=xsecULExpMinus_PFDijet2016;
	width_limit[w][5]=xsecULExpMinus2_PFDijet2016;
      // find max and min for frame
/*      if(mass < medMin)
	medMin = mass;

      if(mass > medMax)
	medMax = mass;
*/
      medMassList.push_back(mass);

      grobs->SetPoint(obscounter, double(width_[w]),xsecULObs_PFDijet2016 );
      obscounter++;

    // 1 sigma dw
      grexp_1sigma_dw->SetPoint(exp_down_counter_1s, double(width_[w]),xsecULExpMinus_PFDijet2016 );      
      exp_down_counter_1s++;
    // 1 sigma up
      grexp_1sigma_up->SetPoint(exp_up_counter_1s, double(width_[w]),xsecULExpPlus_PFDijet2016 );      
      exp_up_counter_1s++;

    // 2 sigma dw
      grexp_2sigma_dw->SetPoint(exp_down_counter_2s, double(width_[w]),xsecULExpMinus2_PFDijet2016 );      
      exp_down_counter_2s++;
    // 2 sigma up
      grexp_2sigma_up->SetPoint(exp_up_counter_2s, double(width_[w]),xsecULExpPlus2_PFDijet2016 );      
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
 ifstream input_1 ("DMV_DM1_WithoutMjjCut.txt");
 if (input_1.fail()) {
   exit(1);
 }
  int zMass = 0; double xsec=0.;  double widthF=0.; double coup=0.; double cutXsec=0; int evnt=0; int evntCut=0;
  int j=0;
  int k=-1;
  double arCoup[10];
  double arXsec[10];
  double arCutXsec[10];
  double arWidthF[10];
  for (int s=0; s<10; ++s){
  arCoup[s]=0;
  arXsec[s]=0;
  arCutXsec[s]=0;
  arWidthF[s]=0;
  }
//REading txt file
  while (!input_1.eof()) {
    input_1 >> zMass >> xsec >> widthF >> coup >> cutXsec >>evnt >>evntCut ;
    if (!input_1.fail()) {
        k++;
        if( zMass==zm){
        arCoup[j]=coup;
        arXsec[j]=xsec;
        arCutXsec[j]=cutXsec;
        arWidthF[j]=widthF;
         j++;
        }
        }
    }

	//sorted g_q coupling constant from 0.1 to 1     
        for (int kk=0; kk<10; ++kk){
                for (int nn=kk+1; nn<10; ++nn){
                        if(arCoup[nn]<arCoup[kk] && arCoup[nn]!=0){
                                double c=arCoup[kk];
                                double x=arXsec[kk];
				double cutx=arCutXsec[kk];
                                double ww=arWidthF[kk];
                                arCoup[kk]=arCoup[nn];
                                arCoup[nn]=c;
                                arXsec[kk]=arXsec[nn];
 				arCutXsec[kk]=arCutXsec[nn];
                                arXsec[nn]=x;
				arCutXsec[nn]=cutx;
                                arWidthF[kk]=arWidthF[nn];
                                arWidthF[nn]=ww;
                        }
                }
        }
     // seted xsec acording to q_q
     for (int tt=0; tt<10; ++tt){
                if (arCutXsec[tt]!=0) grXsec->SetPoint(tt,arWidthF[tt]/zm,arCutXsec[tt]);
                cout <<" arcopling ="<<arCoup[tt]<<" X ="<<arWidthF[tt]/zm <<" xsec Y ="<<arCutXsec[tt]<<endl;

     } 
 
  //////////// All the plotting and cosmetics
  TCanvas* canvas = new TCanvas("canvas", "canvas",600,600);
  TH1* frame = canvas->DrawFrame(medMin,TMath::MinElement(graph_2sigma_band->GetN(),graph_2sigma_band->GetY())*0.01,
				 medMax,TMath::MaxElement(graph_2sigma_band->GetN(),graph_2sigma_band->GetY())*4.5, "");
  //frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->SetTitle("#Gamma/M");
  frame->GetYaxis()->SetTitle("#sigma #times #it{B} #times #it{A} (pb)");
  frame->GetXaxis()->SetTitleOffset(1.15);
  frame->GetYaxis()->SetTitleOffset(1.10);  
  frame->Draw();
  CMS_lumi(canvas,"35.9",false,true,false,0.0,0.0);

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

  //grXsec->SetMarkerStyle(24);
  grXsec->SetLineColor(kBlack);
  grXsec->SetLineWidth(2);
  grXsec->Draw("Lsame*");
  grXsec->SetMarkerStyle(24);

  TF1* line = new TF1 ("line","1",medMin,medMax);
  line->SetLineColor(kBlue);
  line->SetLineWidth(2);
//  line->Draw("L same");

  TLegend *leg = new TLegend(0.175,0.5,0.57,0.77);  
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
  tex->SetTextColor(kBlack);
  tex->SetTextFont(42);
  tex->SetLineWidth(2);
  tex->SetTextSize(0.030);
  tex->Draw();
  tex->DrawLatex(0.175,0.80,("xSec for Z' Mass = "+to_string(int(zm))+" GeV").c_str()); 
  tex->DrawLatex(0.775,0.80,("Spin 1"));
  gPad->RedrawAxis();
  gPad->Modified(); 
  gPad->Update();
  
  canvas->SaveAs((outputDIR+"/PFDMV_qq_limit_observed_Updated35p9fb_Final_1D_M"+to_string(int(zm))+".pdf").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/PFDMV_qq_limit_observed_Updated35p9fb_Final_1D_M"+to_string(int(zm))+".png").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/PFDMV_qq_limit_observed_Updated35p9fb_Final_1D_M"+to_string(int(zm))+".C").c_str(),"pdf");

  canvas->SetLogy();
  frame->GetYaxis()->SetRangeUser(TMath::MinElement(graph_2sigma_band->GetN(),graph_2sigma_band->GetY())*0.01,
				  TMath::MaxElement(graph_2sigma_band->GetN(),graph_2sigma_band->GetY())*100000);
  canvas->SaveAs((outputDIR+"/PFDMV_qq_limit_observed_Updated35p9fb_Final_1D_M"+to_string(int(zm))+"_log.pdf").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/PFDMV_qq_limit_observed_Updated35p9fb_Final_1D_M"+to_string(int(zm))+"_log.png").c_str(),"pdf");

 for (int i=0; i<7; i++){
 	for (int j=0; j<6; j++){
         	cout <<width_[i] <<" = "<< width_limit[i][j]<<endl;
 	}
 }

    double m1=0, c1=0, m2=0, c2=0;
    double x1=0, y1=0, x2=0, y2=0;
    double dx=0, dy=0,dxi=0, dyi=0, c_d=0, in_d=0, x=0;
    double intersection_X=0, intersection_Y=0;
    double x_1=0, y_1=0;
    std::cout << " Program to find the intersecting point of two lines:\n";
    std::cout << "Enter Line1 - X1: "<<arWidthF[arg1]/zm<<endl;
    x1=arWidthF[arg1]/zm;
    x_1=x1;

    std::cout << "Enter Line1 - Y1: "<<arCutXsec[arg1]<<endl;
    y1 =arCutXsec[arg1];
    y_1 = y1;

    std::cout << "Enter Line1 - X2: "<<arWidthF[arg1+1]/zm<<endl;
    x2=arWidthF[arg1+1]/zm;

    std::cout << "Enter Line1 - Y2: "<<arCutXsec[arg1+1]<<endl;
    y2 = arCutXsec[arg1+1];
    dx = x2 - x1;
    dy = y2 - y1;
    m1 = dy / dx;
    c_d= sqrt(dx*dx + dy*dy);
    // y = mx + c

    // intercept c = y - mx
    c1 = y1 - m1 * x1; // which is same as y2 - slope * x2

    x1=0, y1=0, x2=0, y2=0;
    dx=0, dy=0,dxi=0, dyi=0,in_d=0;

    std::cout << "Enter Line2 - X1: "<<width_[arg2]<<endl;
    x1= width_[arg2];

    std::cout << "Enter Line2 - Y1: "<<width_limit[arg2][arg3]<<endl;
    y1=width_limit[arg2][arg3];

    std::cout << "Enter Line2 - X2: "<<width_[arg2+1]<<endl;
    x2= width_[arg2+1];

    std::cout << "Enter Line2 - Y2: "<<width_limit[arg2+1][arg3]<<endl;
    y2=width_limit[arg2+1][arg3];

    dx = x2 - x1;
    dy = y2 - y1;
    m2 = dy / dx;
    c2 = y2 - m2 * x2; // which is same as y2 - slope * x2

    if( (m1 - m2) == 0)
        std::cout << "No Intersection between the lines\n";
    else{
        intersection_X = (c2 - c1) / (m1 - m2);
        intersection_Y = m1 * intersection_X + c1;
        std::cout << "Intersecting Point: = ";
        std::cout << intersection_X;
        std::cout << ",";
        std::cout << intersection_Y<<std::endl;
        std::cout << "coupling d ="<<c_d<<std::endl;
        dxi = intersection_X - x_1;
        dyi = intersection_Y - y_1;
        in_d = sqrt(dxi*dxi + dyi*dyi);
        std::cout << "crossing d ="<<in_d<<std::endl;
        x = (0.1 * in_d)/c_d;
        std::cout << " final coupling value ="<<arCoup[arg1]+x<<std::endl;


    }
}

