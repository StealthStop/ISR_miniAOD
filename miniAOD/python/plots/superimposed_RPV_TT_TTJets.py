import ROOT
import math
import sys
import numpy as np
import array as arr
from ROOT import TFile, gROOT, gStyle

debug = False 

def print_db(input):
    if (debug):
        print input

def main():
    gROOT.SetBatch(True)
    gStyle.SetOptStat(0)
    #gStyle.SetOptStat(1)    

    # ---------------------------------
    # root path & years & histograms
    # ---------------------------------
    pathTT     = "/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/ISR_miniAOD/CMSSW_10_2_9/src/ISR_miniAOD/miniAOD/python/" 
    pathTTJets = "/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/ISR_miniAOD/CMSSW_10_2_9/src/ISR_miniAOD/miniAOD/python/"
    pathRPV    = "/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/ISR_miniAOD/CMSSW_10_2_9/src/ISR_miniAOD/miniAOD/python/"   
 
    years = [
        "2016" ,
        #"2017" ,
        #"2018pre" ,
        #"2018post" ,
    ]

    datasets = [
        "TT" ,
        "TTJets"  ,
        "RPV_350" ,
        #"RPV_550" ,
        #"RPV_850" ,
    ]

    names = {
        "TT"      : "miniAOD_TT" ,
        "TTJets"  : "miniAOD_TTJets" ,
        "RPV_350" : "miniAOD_RPV350" ,
        #"RPV_550" : "miniAOD_RPV550" ,
        #"RPV_850" : "miniAOD_RPV850" ,
    }

    histnames = [
        "miniAOD_ISR_Eta"    ,
        #"miniAOD_ISR_Eta_gg" ,
        #"miniAOD_ISR_Eta_qg" ,
        #"miniAOD_ISR_Eta_gq" ,
        #"miniAOD_ISR_Pt" ,
        #"miniAOD_ISR_Pt_gg" ,
        #"miniAOD_ISR_Pt_qg" ,
        #"miniAOD_ISR_Pt_gq" ,

        "miniAOD_ISR_Eta_g"     ,
        #"miniAOD_ISR_Eta_gg_g"  ,
        #"miniAOD_ISR_Eta_gg_q"  ,
        #"miniAOD_ISR_Eta_qq_g"  , 
        #"miniAOD_ISR_Eta_qq_q"  , 
        #"miniAOD_ISR_Eta_gq_gq" ,
        #"miniAOD_ISR_Pt_g"     ,
        #"miniAOD_ISR_Pt_gg_g"  ,
        #"miniAOD_ISR_Pt_gg_q"  ,
        #"miniAOD_ISR_Pt_qq_g"  ,
        #"miniAOD_ISR_Pt_qq_q"  ,
        #"miniAOD_ISR_Pt_gq_gq" ,


        # dau-mom
        #"miniAOD_ISR_Eta_gudsc" ,
        #"miniAOD_ISR_Eta_gb" ,
        #"miniAOD_ISR_Eta_bg" ,
        #"miniAOD_ISR_Eta_gg" ,
        #"miniAOD_ISR_Eta_udscg" , 
        #"miniAOD_ISR_Eta_PtDiv_0to100_udscg" ,
        #"miniAOD_ISR_Eta_PtDiv_100to200_udscg" ,
        #"miniAOD_ISR_Eta_PtDiv_g200_udscg" ,
        #"miniAOD_ISR_Eta_neg_udscg" ,
        #"miniAOD_ISR_Eta_anti_ug" ,
        #"miniAOD_ISR_Eta_anti_dg" ,
        #"miniAOD_ISR_Eta_anti_sg" ,
        #"miniAOD_ISR_Eta_anti_cg" ,
        #"miniAOD_ISR_Eta_pos_udscg" ,
        #"miniAOD_ISR_Eta_ug" ,
        #"miniAOD_ISR_Eta_dg" ,
        #"miniAOD_ISR_Eta_sg" ,
        #"miniAOD_ISR_Eta_cg" ,

        #"miniAOD_ISR_Pt_udscg" ,
        #"miniAOD_ISR_divPt_udscg" ,
    ]
    
    varList = [
        "Eta" ,
        #"Eta_gg" ,
        #"Eta_qg" ,
        #"Eta_gq" ,
        #"Pt" ,
        #"Pt_gg" ,
        #"Pt_qg" ,
        #"Pt_gq" ,

        "Eta_g" ,
        #"Eta_g_gg" ,
        #"Eta_g_qg" ,
        #"Eta_g_gq" ,
        #"Pt_g" ,
        #"Pt_g_gg" ,
        #"Pt_g_qg" ,
        #"Pt_g_gq" ,


        #"Eta_gudsc" ,
        #"Eta_gb",
        #"Eta_bg" ,
        #"Eta_gg" , 
        #"Eta_udscg" ,  
        #"Eta_PtDiv_0to100_udscg" ,
        #"Eta_PtDiv_100to200_udscg" ,
        #"Eta_PtDiv_g200_udscg" ,
        #"Eta_neg_udscg" ,
        #"Eta_anti_ug" ,
        #"Eta_anti_dg" ,
        #"Eta_anti_sg" ,
        #"Eta_anti_cg" ,   
        #"Eta_pos_udscg" ,
        #"Eta_ug" ,
        #"Eta_dg" ,
        #"Eta_sg" ,
        #"Eta_cg" ,

        #"Pt_udscg" ,
        #"divPt_udscg" ,
    ]

    # -----------------------------
    # loop over & get histograms 
    # -----------------------------
    for year in years:
        print_db("Processing year " + year)
            
        f = 0; f1 = 0; f2 = 0; f3 = 0; f4 = 0

        for dataset in datasets:
            print_db("Processing dataset " + dataset)    
    
            if dataset == "TT":
                filename = pathTT + year + "_" + names[dataset] + ".root"
                print_db("Opening file " + filename)
                f = ROOT.TFile.Open(filename, "READ")

            if dataset == "TTJets":
                filename1 = pathTTJets + year + "_" + names[dataset] + ".root"
                print_db("Opening file " + filename1)
                f1 = ROOT.TFile.Open(filename1, "READ")          
 
            elif dataset == "RPV_350": 
                filename2 = pathRPV + year + "_" + names[dataset] + ".root"
                print_db("Opening file " + filename2)
                f2 = ROOT.TFile.Open(filename2, "READ")           
 
            #elif dataset == "RPV_550":
            #    filename3 = pathRPV + year + "_" + names[dataset] + ".root"
            #    print_db("Opening file " + filename3)
            #    f3 = ROOT.TFile.Open(filename3, "READ")            

            #elif dataset == "RPV_850":
            #    filename4 = pathRPV + year + "_" + names[dataset] + ".root"
            #    print_db("Opening file " + filename4)
            #    f4 = ROOT.TFile.Open(filename4, "READ")    

        # get the histograms in separate loop
        for histname in histnames:
            print_db("Processing histogram " + histname)
            
            for var in varList:
                ROOT.TH1.AddDirectory(0)
            
                canvas     = ROOT.TCanvas("c_%s"%(var), "c_%s"%(var), 1400, 1400)
                #new_legend = ROOT.TLegend(0.60, 0.8, 0.9, 0.90)
                new_legend = ROOT.TLegend(0.50, 0.8, 0.9, 0.90)

                # open the histograms
                hISR_TT     = f.Get("demo/" + "miniAOD_ISR_%s"%(var))
                hISR_TTJets = f1.Get("demo/" + "miniAOD_ISR_%s"%(var))            
                hISR_RPV350 = f2.Get("demo/" + "miniAOD_ISR_%s"%(var))
                #hISR_RPV550 = f3.Get("demo/" + "miniAOD_ISR_%s"%(var))
                #hISR_RPV850 = f4.Get("demo/" + "miniAOD_ISR_%s"%(var))
 

                # define the histograms
                hISR_TT.SetLineColor(40)
                hISR_TT.SetLineWidth(3)
                hISR_TT.SetFillColor(40)
                hISR_TT.SetFillStyle(3007)
                hISR_TT.RebinX(2)
                hISR_TT.Scale(1.0/hISR_TT.Integral())

                hISR_TTJets.SetLineColor(38)
                hISR_TTJets.SetLineWidth(3)
                hISR_TTJets.SetFillColor(38)
                hISR_TTJets.SetFillStyle(3007)
                hISR_TTJets.RebinX(2)
                hISR_TTJets.Scale(1.0/hISR_TTJets.Integral())

                hISR_RPV350.SetLineColor(42) # 29
                hISR_RPV350.SetLineWidth(3)
                hISR_RPV350.SetFillColor(42)
                hISR_RPV350.SetFillStyle(3002)
                hISR_RPV350.RebinX(2) 
                hISR_RPV350.Scale(1.0/hISR_RPV350.Integral())
                
                #hISR_RPV550.SetLineColor(42)
                #hISR_RPV550.SetLineWidth(3)
                #hISR_RPV550.SetFillColor(42)
                #hISR_RPV550.SetFillStyle(3002)
                #hISR_RPV550.RebinX(2)
                #hISR_RPV550.Scale(1.0/hISR_RPV550.Integral())

                #hISR_RPV850.SetLineColor(42) # 41
                #hISR_RPV850.SetLineWidth(3)
                #hISR_RPV850.SetFillColor(42)
                #hISR_RPV850.SetFillStyle(3002)
                #hISR_RPV850.RebinX(2)
                #hISR_RPV850.Scale(1.0/hISR_RPV850.Integral())


                # change the hist if the axis range larger for which one
                hISR_RPV350.GetXaxis().SetTitle(var)
                #hISR_TTJets.GetXaxis().SetTitle(var)
                #hISR_TT.GetXaxis().SetRangeUser(0,600) # mass 100, pt 500
                #hISR_TTJets.SetTitle(year + " " + dataset + " " + var + " ")
                hISR_RPV350.SetTitle(year + "_miniAOD_GenISR" + " ")
                #hISR_TTJets.GetXaxis().SetBinLabel(1, "0-100 GeV")
                #hISR_TTJets.GetXaxis().SetBinLabel(2, "100-200 GeV") 
                #hISR_TTJets.GetXaxis().SetBinLabel(3, ">200 GeV")                

                # get the legend
                new_legend.AddEntry(hISR_TT,     "GenISR - TT - MEAN: %f" %(hISR_TT.GetMean()), "L")
                new_legend.AddEntry(hISR_TTJets, "GenISR - TTJets - MEAN: %f" %(hISR_TTJets.GetMean()), "L")
                new_legend.AddEntry(hISR_RPV350, "GenISR - RPV350 - MEAN: %f" %(hISR_RPV350.GetMean()), "L")
                #new_legend.AddEntry(hISR_RPV550, "GenISR - RPV550 - MEAN: %f" %(hISR_RPV550.GetMean()), "L")
                #new_legend.AddEntry(hISR_RPV850, "GenISR - RPV850 - MEAN: %f" %(hISR_RPV850.GetMean()), "L")    
                new_legend.SetTextSize(0.018)

                #hISR_TTJets.SetStats(ROOT.kTRUE)
                #hISR_RPV350.SetStats(ROOT.kTRUE)

                # change the order if the axis range is larger for which one
                canvas.cd()
                hISR_RPV350.Draw("E HIST")
                hISR_TTJets.Draw("SAME E HIST")
                #ROOT.gPad.Modified()
                #ROOT.gPad.Update()
                #st = hISR_TTJets.GetListOfFunctions().FindObject("stats")
                #st.__class__ = ROOT.TPaveStats
                #st.SetX1NDC(0.78)
                #st.SetX2NDC(0.98)
                #st.SetY1NDC(0.6)
                #st.SetY2NDC(0.7)
                #st.SetBorderSize(3)
                #st.SetLineColor(38)
                hISR_TT.Draw("SAMES E HIST")
                #ROOT.gPad.Modified()
                #ROOT.gPad.Update()
                #st = hISR_RPV350.GetListOfFunctions().FindObject("stats")
                #st.__class__ = ROOT.TPaveStats
                #st.SetX1NDC(0.78)
                #st.SetX2NDC(0.98)
                #st.SetY1NDC(0.49)
                #st.SetY2NDC(0.59)
                #st.SetBorderSize(3)
                #st.SetLineColor(42)
                #hISR_RPV550.Draw("SAME E HIST")
                #hISR_RPV850.Draw("SAME E HIST")
                new_legend.Draw("SAME")
                
                # save the canvas
                canvas.SaveAs("miniAOD_ISR_Eta_brief/" + year + "_miniAOD_GenISR_" + var + "_RPV350_TTJets_TT_normalized" + ".pdf")
                #canvas.SaveAs("miniAOD_ISR_eta/" + year + "_" + var + "_RPV_TTJets" + ".pdf")
                #canvas.SaveAs("miniAOD_ISR_eta/normalized/" + year + "_" + var + "_RPV_TTJets_normalized" + ".pdf")
                canvas.Close()

        f.Close()
        f1.Close()
        f2.Close()
        #f3.Close()
        #f4.Close()
 
if __name__ == "__main__":
    main()

