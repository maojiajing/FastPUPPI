#!/usr/bin/env python
# encoding: utf-8

# File        : GenJetReader.py
# Author      : Zhenbin Wu
# Contact     : zhenbin.wu@gmail.com
# Date        : 2017 Sep 28
#
# Description : 

from __future__ import print_function
from rootpy.plotting import Hist,Canvas
import ROOT 

class GenJetReader():
    def __init__(self, tree, ctype='ntuple'):
        self.tree_ = tree
        self.folder_ = ctype
        self.hismap = {}
        self.BookHistogram()
        self.topjets = []
        self.Wjets = []
        self.Lepjets = []

    def BookHistogram(self):
        self.hismap["GenTopEta"]=  Hist(100, -5, 5, title= "GenTopEta;GenTopEta;Event" )
        self.hismap["GenTopEta_Pt200"]=  Hist(100, -5, 5,  title="GenTopEta200;GenTopEta Pt>200;Event" )
        self.hismap["GenTopPt"]=  Hist(100, 0, 1000)
        self.hismap["Tau32"]=  Hist(100, 0, 10)
        self.hismap["Tau32_Pt200"]=  Hist(100, 0, 10)

    def Run(self):
        self.topjets = self.GetGenJets([6])
        self.Wjets = self.GetGenJets([24])
        self.Lepjets = self.GetGenJets([11, 13, 15])
        self.PlotTop()


    def GetGenJets(self, pdgid):
        idx = []
        rejets = []
        for i in range(len(self.tree_.GenPar_ID)):
            if abs(self.tree_.GenPar_ID.at(i)) in pdgid:
                idx.append(i)
        for i in idx:
            j = ROOT.TLorentzVector(0, 0, 0, 0)
            j.SetPtEtaPhiM(self.tree_.GenPar_pt.at(i),
                           self.tree_.GenPar_eta.at(i),
                           self.tree_.GenPar_phi.at(i), 0)
            rejets.append(j)
        return rejets
    
    def PlotTop(self):
        if len(self.topjets) == 0:
            return None
        for j in self.topjets:
            self.hismap["GenTopEta"].Fill(j.Eta())
            self.hismap["GenTopPt"].Fill(j.Pt())
            if j.Pt() > 200:
                self.hismap["GenTopEta_Pt200"].Fill(j.Eta())

    def GetHist(self):
        return self.hismap
