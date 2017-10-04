#!/usr/bin/env python
# encoding: utf-8

# File        : Reader.py
# Author      : Zhenbin Wu
# Contact     : zhenbin.wu@gmail.com
# Date        : 2017 Sep 26
#
# Description : 

from __future__ import print_function
from __future__ import division
from rootpy.plotting import Hist,Canvas
from rootpy.plotting import Hist,Canvas
from rootpy.plotting.style import get_style
from rootpy.interactive import wait
import ROOT

class AK8JetReader():
    def __init__(self, tree, ctype):
        self.tree_ = tree
        self.folder_ = ctype
        self.hismap = {}
        self.BookHistogram()
        self.Jets = []

    def BookHistogram(self):
        self.hismap["FatJetPt"]=     Hist(100, 0, 1000, title="FatJetPt;Jet Pt;Event")
        self.hismap["FatJetMass"]=   Hist(100, 0, 200,  title="FatJetMass;Jet Mass;Event")
        self.hismap["Tau32"]=        Hist(10,  0, 1,    title="Tau32;Tau32;No.  of Jets")
        self.hismap["Tau32Pt200"]=   Hist(10,  0, 1,    title="Tau32; Tau32 (p_T > 200GeV);No. of Jets")
        self.hismap["Tau32Pt400"]=   Hist(10,  0, 1,    title="Tau32; Tau32 (p_T > 400GeV);No. of Jets")
        self.hismap["Tau21Pt200"]=   Hist(10,  0, 1,    title="Tau21; Tau21 (p_T > 200GeV);No. of Jets")
        self.hismap["t_Tau32"]=      Hist(10,  0, 1,    title="Tau32;Tau32_t;No.  of Jets")
        self.hismap["t_Tau32Pt400"]= Hist(10,  0, 1,    title="Tau32; Tau32_t (p_T > 400GeV);No. of Jets")
        self.hismap["w_Tau21"]=      Hist(10,  0, 1,    title="Tau21;Tau21_w;No.  of Jets")
        self.hismap["w_Tau21Pt200"]= Hist(10,  0, 1,    title="Tau21; Tau21_w (p_T > 200GeV);No. of Jets")

    def Loop(self):
        for event in self.tree_:
            self.Run(event)

    def Run(self):
        e = self.tree_
        self.Jets = []
        if len(e.FatJet_pt) > 0:
            self.hismap["FatJetPt"].Fill(e.FatJet_pt.at(0))
            self.hismap["FatJetMass"].Fill(e.FatJet_mass.at(0))

        for i in range(len(e.FatJet_pt)):
            j = ROOT.TLorentzVector(0, 0, 0, 0)
            j.SetPtEtaPhiM(self.tree_.FatJet_pt.at(i), self.tree_.FatJet_eta.at(i),
                           self.tree_.FatJet_phi.at(i), self.tree_.FatJet_mass.at(i))
            self.Jets.append(j)
            if self.tree_.Tau2.at(i) != 0:
                self.hismap["Tau32"].Fill(e.Tau3.at(i)/ e.Tau2.at(i))
                if e.FatJet_pt.at(i) > 200:
                    self.hismap["Tau32Pt200"].Fill(e.Tau3.at(i)/ e.Tau2.at(i))
                if e.FatJet_pt.at(i) > 400:
                    self.hismap["Tau32Pt400"].Fill(e.Tau3.at(i)/ e.Tau2.at(i))
            if self.tree_.Tau1.at(i) != 0:
                if e.FatJet_pt.at(i) > 200:
                    self.hismap["Tau21Pt200"].Fill(e.Tau2.at(i)/ e.Tau1.at(i))

    def CheckWithGenJet(self, genObjjet, genlepjet):
        topjet = set()
        for i in range(len(self.tree_.FatJet_pt)):
            j = self.Jets[i]
            for top in genObjjet:
                if top.DeltaR(j) < 0.4:
                    topjet.add(i)
            for lep in genlepjet:
                if lep.DeltaR(j) < 0.4 and i in topjet:
                    topjet.remove(i)
        return topjet

    def PlotTopMatchedJet(self, genTopjet, genlepjet):
        for i in self.CheckWithGenJet(genTopjet, genlepjet):
            j = self.Jets[i]
            self.hismap["t_Tau32"].Fill(self.tree_.Tau3.at(i)/ self.tree_.Tau2.at(i)
                                        if self.tree_.Tau2.at(i)  else 0)
            if j.Pt() > 400:
                self.hismap["t_Tau32Pt400"].Fill(self.tree_.Tau3.at(i)/ self.tree_.Tau2.at(i)
                                        if self.tree_.Tau2.at(i)  else 0)


    def PlotWMatchedJet(self, genWjet, genlepjet):
        for i in self.CheckWithGenJet(genWjet, genlepjet):
            j = self.Jets[i]
            self.hismap["w_Tau21"].Fill(self.tree_.Tau2.at(i)/ self.tree_.Tau1.at(i)
                                        if self.tree_.Tau1.at(i)  else 0)
            if j.Pt() > 200:
                self.hismap["w_Tau21Pt200"].Fill(self.tree_.Tau2.at(i)/ self.tree_.Tau1.at(i)
                                        if self.tree_.Tau1.at(i)  else 0)


    def Test(self):
        print(self.tree_.mc_pt)

    def Draw(self):
        get_style('ATLAS')
        canvas = Canvas(width=700, height=500)
        canvas.SetLeftMargin(0.15)
        canvas.SetBottomMargin(0.15)
        canvas.SetTopMargin(0.10)
        canvas.SetRightMargin(0.05)
        for k,v in self.hismap.items():
            canvas.Clear()
            v.Draw()
            wait()

    def GetHist(self):
        return self.hismap
