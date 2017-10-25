from ROOT import *
import ROOT,sys,os,math
from math import *
import numpy as np
#Run
#python wide_plot_gq.py

#from mass_plot import *

gROOT.ForceStyle()
gStyle.SetLegendBorderSize(0)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(0)
gStyle.SetPadLeftMargin(0.10)
gStyle.SetPadBottomMargin(0.12)
gStyle.SetPadTopMargin(0.075)
gStyle.SetPadRightMargin(0.10)

def SetCurveStyle(gr,color,lineStyle,lineWidth,markerSize):
    gr.SetMarkerColor(color)
    gr.SetMarkerSize(markerSize)
    gr.SetLineColor(color)
    gr.SetLineWidth(lineWidth)
    gr.SetFillStyle(lineStyle)
    gr.SetFillColor(color)
    return 

def g_q_res(gq):
    #return gq*math.pow(0.5+math.pow(0.25+1/(18*(gq**2)),0.5),0.5)
    return math.pow(3*gq**2+math.pow(gq**2*(9*gq**2+2),0.5),0.5)/math.pow(6,0.5)

def medWidth(gq):
  return 9*gq**2/(4*3.141592653)+1/(12*3.141592653)

def medWidth_res(gq):
    return 6*gq**2/(4*3.141592653)

if __name__=="__main__":

    isGQ=False  ### Should always be false. If true, plot gq instead of gq_prime. ###

    isLogx=False   

    doBoundary=False  ### plot full boundaries for the exclusion limits, not useful now, should always be false ###

    file_wide=TFile("FinalPlotDM_V/grFinal_PFDMV-g_q-wideWithNarrow-limit_35p9fb.root")

    if isGQ:
        med_min=600
        med_max=3700
    else:
        med_min=1600
        med_max=4100
    med_step=10
    
    exp_wide=file_wide.Get("graph_expected")
    obs_wide=file_wide.Get("graph_observed")    
    wide_All=file_wide.Get("graph_All")
    wide_All_leg=file_wide.Get("My_leg")
    
    ### dijet wide resonance search ###

    new_exp_wide=TGraph(0)
    for mmed in np.linspace(med_min,med_max,num=int((med_max-med_min)/med_step+1)):
        new_exp_wide.SetPoint(new_exp_wide.GetN(),mmed,g_q_res(exp_wide.Eval(mmed)))
    SetCurveStyle(exp_wide,kBlack,2,402,0.1)
    SetCurveStyle(new_exp_wide,kBlack,2,402,0.1)

    new_obs_wide=TGraph(0)
    for mmed in np.linspace(med_min,med_max,num=int((med_max-med_min)/med_step+1)):
        new_obs_wide.SetPoint(new_obs_wide.GetN(),mmed,g_q_res(obs_wide.Eval(mmed)))
    SetCurveStyle(obs_wide,kRed,1,402,0.1)
    SetCurveStyle(new_obs_wide,kRed,1,402,0.1)

    ### Plotting ###
    
    canvas=TCanvas("myCanvas","myCanvas",0,0,640,610)
    xmin=1600
    xmax=4100
    if isLogx:
        canvas.SetLogx()
        xmin=40
        xmax=8000

    mg=TMultiGraph()

    if isGQ:
        mg.Add(new_exp_wide,"l")
        mg.Add(new_obs_wide,"l")
    else:
        mg.Add(exp_wide,"l")
        mg.Add(obs_wide,"l")

    ymin=0.
    if isGQ:
        ymax=1.42788
    else:
        ymax=0.755
    
    #mg.Draw("apl")
    #mg.GetYaxis().SetRangeUser(ymin,ymax)
    wide_All.Draw("3Lsame")
    wide_All_leg.Draw("same")
    wide_All.GetYaxis().SetRangeUser(ymin,ymax) 
    if isGQ:
        mg.GetYaxis().SetTitle("g_{q}")
    else:
        wide_All.GetYaxis().SetTitle("g_{q}")
    wide_All.GetYaxis().SetTitleOffset(0.8)
    wide_All.GetYaxis().SetTitleSize(0.05)
    wide_All.GetYaxis().SetLabelSize(0.04)
    
    wide_All.GetXaxis().SetLimits(xmin,xmax)
    wide_All.GetXaxis().SetTitle("M_{Med} [GeV]")
    wide_All.GetXaxis().SetTitleOffset(1.)
    wide_All.GetXaxis().SetTitleSize(0.05)
    wide_All.GetXaxis().SetLabelSize(0.04)

    if isLogx:
        #mg.GetXaxis().SetNdivisions(10)
        mg.GetXaxis().SetMoreLogLabels()
        mg.GetXaxis().SetMoreLogLabels()
        #mg.GetXaxis().GetLabels().ChangeLabel(7000,"7000")
        leglabel=TLatex(6900,-0.066,"8000")
        leglabel.SetTextSize(0.04)
        leglabel.SetTextFont(42)
        leglabel.Draw()
        mg.GetXaxis().SetNdivisions(10)
        mg.GetXaxis().SetNoExponent(True); 
    
    minwidth=medWidth(ymin)
    print minwidth,
    myFunc=TF1("myFunc","pow((x-1/(12*3.141592653))*(4*3.141592653)/6,0.5)",minwidth,0.3)
    y2=TGaxis(xmax, ymin, xmax, ymax,"myFunc",510, "+L")
    y2.SetTitle("#Gamma/M_{Med}")
    y2.SetLabelSize(0.035)
    y2.SetTitleSize(0.035)
    y2.SetTitleOffset(1.4)
    y2.Draw()
    canvas.RedrawAxis();
    canvas.Modified();
    canvas.Update();
    ### plot full boundaries for the exclusion limits, not useful now ###

    if doBoundary:
    
        upper_exp_wide=TGraph(0)
        upper_exp_wide.SetPoint(upper_exp_wide.GetN(),600,exp_wide.Eval(600))
        upper_exp_wide.SetPoint(upper_exp_wide.GetN(),600,0.5)
        upper_exp_wide.SetPoint(upper_exp_wide.GetN(),3700,0.5)
        #upper_exp_wide.SetPoint(upper_exp_wide.GetN(),3700,exp_wide.Eval(3700))
        #SetCurveStyle(upper_exp_wide,kRed-10,3004,-402,0.1)
        upper_exp_wide.Draw()

        upper_obs_wide=TGraph(0)
        upper_obs_wide.SetPoint(upper_obs_wide.GetN(),600,obs_wide.Eval(600))
        upper_obs_wide.SetPoint(upper_obs_wide.GetN(),600,0.5)
        upper_obs_wide.SetPoint(upper_obs_wide.GetN(),3700,0.5)
        #upper_obs_wide.SetPoint(upper_obs_wide.GetN(),3700,obs_wide.Eval(3700))
        #SetCurveStyle(upper_obs_wide,kRed+2,3004,-402,0.1)
        upper_obs_wide.Draw("same")

    if isLogx:
        position=0.12
    else:
        position=0.84
        
    if isGQ:
        l0p5=TLine(xmin,0.3,xmax,0.3)
        l0p5.SetLineColor(kGray+1)
        l0p5.SetLineStyle(kDashed)
        l0p5.Draw("same")
        l0p5T=TLatex((xmax-xmin)*position+xmin,0.5+0.05,"g_{q}=0.5, #Gamma/M_{Med}=15%")
        l0p5T.SetTextSize(0.035)
        l0p5T.SetTextColor(kGray+1)
        l0p5T.Draw("same")
    else:
        l0p5=TLine(xmin,0.4,xmax,0.4)
        l0p5.SetLineColor(kGray+1)
        l0p5.SetLineStyle(kDashed)
        #l0p5.Draw("same")
        l0p5T=TLatex((xmax-xmin)*position+xmin,0.4+0.03," #Gamma/M_{Med}=10%")
        l0p5T.SetTextSize(0.035)
        l0p5T.SetTextColor(kGray+1)
        #l0p5T.Draw("same")


    l1p0=TLine(xmin,0.77,xmax,0.77)
    l1p0.SetLineColor(kGray+1)
    l1p0.SetLineStyle(kDotted)
    #l1p0.Draw("same")
    if isGQ:
        l1p0T=TLatex((xmax-xmin)*position+xmin,1+0.05,"g_{q}=1.0, #Gamma/M_{Med}=50%")
    else:
        l1p0T=TLatex((xmax-xmin)*position+xmin,0.77+0.03," #Gamma/M_{Med}=30%")
    l1p0T.SetTextSize(0.035)
    l1p0T.SetTextColor(kGray+1)
    #l1p0T.Draw("same")

    #leg=TLegend(0.11,0.43,0.5,0.78,"CMS 95% CL Upper Limits")
    leg=TLegend(0.18,0.835,0.45,0.78," 95% CL Upper Limits")
    leg.SetFillStyle(0)
    leg.SetTextSize(0.042)
    
    leg.Draw()

    if isLogx:
        offset=10
    else:
        offset=200
    if isGQ:
        leg1=TLatex(xmin+offset,ymax-0.18,"#splitline{Vector/Axial-Vector Mediator}{m_{DM} = 1 GeV, g_{DM} = 1.0}")
    else:
        #leg1=TLatex(xmin+offset,ymax-0.18,"#splitline{Vector/Axial-Vector Mediator}{#Gamma_{DM} = 0}")
        leg1=TLatex(xmin+offset,ymax-0.12,"Vector/Axial-Vector Mediator")
        leg1p=TLatex(xmin+offset,ymax-0.21,"#Gamma_{DM} = 0")
        leg1p.SetTextSize(0.04)
        #leg1p.Draw("same")
    #leg1.SetTextFont(42)
    leg1.SetTextSize(0.04)
    #leg1.Draw("same")

    # CMS
    leg2=TLatex(xmin+300,ymax-0.07,"#bf{CMS}")
    leg2.SetTextFont(42)
    leg2.SetTextSize(0.045)

    if isLogx:
        lumioffset=5250
    else:
        lumioffset=810
    # lumi
    leg3=TLatex(xmax-lumioffset,ymax+0.01,"36 fb^{-1} (13 TeV)")
    leg3.SetTextFont(42)
    leg3.SetTextSize(0.04)
    leg2.Draw("same")
    leg3.Draw("same")
    
    
    if isGQ:
        if isLogx:
            canvas.SaveAs("dijet_combined_gq_logx.pdf")
        else:
            canvas.SaveAs("dijet_combined_gq.pdf")
    else:
        if isLogx:
            canvas.SaveAs("dijet_combined_gq_logx.pdf")
        else:
            canvas.SaveAs("wide-g_q-withNarrow-limit_35p9fb.pdf")

    
    
