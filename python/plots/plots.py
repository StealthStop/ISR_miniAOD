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
    #gStyle.SetOptStat(0)
    gStyle.SetOptStat(1)    
 
    # ---------------------------------
    # root path & years & histograms
    # --------------------------------- 
    histnames = [
        # mother:proton - daughter:quark - status:71
        #"NumDau_pq_73" ,
        #"MomPt_pq_73"  ,
        #"DauPt_pq_73"  ,
        
        # mother:proton - daughter:gluon - status:71
        #"NumDau_pg_73" ,
        #"MomPt_pg_73"  ,
        #"DauPt_pg_73"  ,

        # mother:quark - daughter:quark - status:23,44,71
        #"NumDau_qq_23" ,
        #"MomPt_qq_23"  ,
        #"DauPt_qq_23"  ,
        #"NumDau_qq_44" ,
        #"MomPt_qq_44"  ,
        #"DauPt_qq_44"  ,
        #"NumDau_qq_73" ,
        #"MomPt_qq_73"  ,
        #"DauPt_qq_73"  ,        

        ## mother:quark - daughter:gluon - status:23,71
        #"NumDau_qg_23" ,
        #"MomPt_qg_23"  ,
        #"DauPt_qg_23"  ,
        #"NumDau_qg_73" ,
        #"MomPt_qg_73"  ,
        #"DauPt_qg_73"  ,

        ## mother:gluon - daughter:gluon - status:23,44,71
        #"NumDau_gg_23" ,
        "MomPt_gg_23"  ,
        #"DauPt_gg_23"  ,
        #"NumDau_gg_44" ,
        #"MomPt_gg_44"  ,
        #"DauPt_gg_44"  ,
        #"NumDau_gg_73" ,
        #"MomPt_gg_73"  ,
        #"DauPt_gg_73"  ,

        ## mother:gluon - daughter:quark - status:23,71
        #"NumDau_gq_23" ,
        "MomPt_gq_23"  ,
        #"DauPt_gq_23"  ,   
        #"NumDau_gq_73" ,
        #"MomPt_gq_73"  ,
        #"DauPt_gq_73"  ,
    ]
    
    # -----------------------------
    # loop over & get histograms 
    # -----------------------------
    f = ROOT.TFile.Open("../miniAOD.root", "READ")

    for histname in histnames:
           
        canvas = ROOT.TCanvas("c1", "c1", 1400, 1000)
        canvas.SetLogy()
        ROOT.TH1.AddDirectory(0)
       
        hist = f.Get("demo/" + histname)
 
        #xmin = 999.0
        #xmax = -0.999
        #if (xmin < xmax):
        #    hist.GetXaxis().SetRangeUser(xmin, xmax)
   
        #ymin = 999.0
        #ymax = -0.999
        #if (ymin < ymax):
        #    hist.GetYaxis().SetRangeUser(ymin, ymax)            
   
        hist.SetLineColor(30) 
        hist.SetLineWidth(2)
        #hist.RebinX(20)
        hist.GetXaxis().SetRangeUser(0, 50)
        hist.Draw("hist E")
        hist.SetTitle(histname)
        hist.GetXaxis().SetTitle("Number of Daughter")
        hist.GetYaxis().SetTitle("Events")

        canvas.SaveAs("plot/" + histname + ".pdf")
        canvas.Close()

    f.Close()


 
if __name__ == "__main__":
    main()

