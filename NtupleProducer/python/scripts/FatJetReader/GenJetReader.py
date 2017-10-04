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

    def BookHistogram(self):
        self.hismap["FatJetPt"]=  Hist(100, 0, 1000)
        self.hismap["Tau32"]=  Hist(100, 0, 10)
        self.hismap["Tau32_Pt200"]=  Hist(100, 0, 10)

    def Run(self):
        return None
        e = self.tree_
        print('-----')
        [print(i) for i in e.GenPar_ID]
        print('-----')


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
            # j.SetPtEtaPhiM(self.tree_.GenPar_pt.at(i), self.tree_.GenPar_eta.at(i), self.tree_.GenPar_phi.at(i), self.tree_.GenPar_mass.at(i))
            rejets.append(j)
        return rejets
