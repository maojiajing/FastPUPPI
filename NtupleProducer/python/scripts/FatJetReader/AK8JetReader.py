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
        self.Nevent = 0
        self.ptorder =[]

    def BookHistogram(self):
        ## Rate study
        self.hismap["NEvent"]=       Hist(2, 0, 2, title="Event;Event;Event")
        self.hismap["RatePT"]=       Hist(100, 0, 1000, title="FatJetPt Rate;Jet Pt;Rate [kHz]")
        self.hismap["RatePT_M100"]= Hist(100, 0, 1000, title="FatJetPt Rate;Jet Pt, M > 100; Rate [kHz]")
        self.hismap["RatePT_M40"]= Hist(100, 0, 1000, title="FatJetPt Rate;Jet Pt, M > 40; Rate [kHz]")
        self.hismap["RatePT_Eta24"]= Hist(100, 0, 1000, title="FatJetPt Rate;Jet Pt, |eta| < 2.4; Rate [kHz]")
        self.hismap["RatePT_Eta3"]= Hist(100, 0, 1000, title="FatJetPt Rate;Jet Pt, |eta| < 3; Rate [kHz]")
        self.hismap["RatePT_Eta3_M40"]= Hist(100, 0, 1000, title="FatJetPt Rate;Jet Pt, |eta| < 3, M > 40; Rate [kHz]")

        ##                           Jet study
        self.hismap["NJetPt30"]      =   Hist(100, 0,  100,  title = "NFatJetPt30;No. Of Jet Pt>30;Event")
        self.hismap["FirstJetPt"]    =   Hist(100, 0,  1000, title = "FirstJetPt;Leading Jet Pt;Event")
        self.hismap["FirstJetMass"]  =   Hist(100, 0,  200,  title = "FirstJetMass;Leading Jet Mass;Event")
        self.hismap["FirstJetEta"]   =   Hist(100, -5, 5,    title = "FirstJetEta;Leading Jet Eta;Event")
        self.hismap["SecondJetPt"]   =   Hist(100, 0,  1000, title = "SecondJetPt;2nd Leading Jet Pt;Event")
        self.hismap["SecondJetMass"] =   Hist(100, 0,  200,  title = "SecondJetMass;2nd Leading Jet Mass;Event")
        self.hismap["SecondJetEta"]  =   Hist(100, -5, 5,    title = "SecondJetEta;2nd Leading Jet Eta;Event")
        self.hismap["ThirdJetPt"]    =   Hist(100, 0,  1000, title = "ThirdJetPt;3rd Leading Jet Pt;Event")
        self.hismap["ThirdJetMass"]  =   Hist(100, 0,  200,  title = "ThirdJetMass;3rd Leading Jet Mass;Event")
        self.hismap["ThirdJetEta"]   =   Hist(100, -5, 5,    title = "ThirdJetEta;3rd Leading Jet Eta;Event")
        self.hismap["ForthJetPt"]    =   Hist(100, 0,  1000, title = "ForthJetPt;4th Leading Jet Pt;Event")
        self.hismap["ForthJetMass"]  =   Hist(100, 0,  200,  title = "ForthJetMass;4th Leading Jet Mass;Event")
        self.hismap["ForthJetEta"]   =   Hist(100, -5, 5,    title = "ForthJetEta;4th Leading Jet Eta;Event")

        self.hismap["FirstJdPt"]    =   Hist(100, -100,  100, title = "FirstJdPt;Leading Jet dPt;Event")
        self.hismap["FirstJdMass"]  =   Hist(100, -50,  50,  title = "FirstJdMass;Leading Jet dMass;Event")
        self.hismap["FirstJdEta"]   =   Hist(100, -1, 1,    title = "FirstJdEta;Leading Jd Eta;Event")
        self.hismap["FirstJdR"]   =   Hist(60, 0, 6,    title = "FirstJdR;Leading Jd dR;Event")

        self.hismap["SecondJdPt"]   =   Hist(100, -100,  100, title = "SecondJdPt;2nd Leading Jet dPt;Event")
        self.hismap["SecondJdMass"] =   Hist(100, -50,  50,  title = "SecondJdMass;2nd Leading Jet dMass;Event")
        self.hismap["SecondJdEta"]  =   Hist(100, -1, 1,    title = "SecondJdEta;2nd Leading Jet dEta;Event")
        self.hismap["SecondJdR"]  =   Hist(60, 0, 6,    title = "SecondJdEta;2nd Leading Jet dR;Event")

        self.hismap["ThirdJdPt"]    =   Hist(100, -100,  100, title = "ThirdJdPt;3rd Leading Jet dPt;Event")
        self.hismap["ThirdJdMass"]  =   Hist(100, -50,  50,  title = "ThirdJdMass;3rd Leading Jet dMass;Event")
        self.hismap["ThirdJdEta"]   =   Hist(100, -1, 1,    title = "ThirdJdEta;3rd Leading Jet dEta;Event")
        self.hismap["ThirdJdR"]   =   Hist(60, 0, 6,    title = "ThirdJdEta;3rd Leading Jet dR;Event")

        self.hismap["ForthJdPt"]    =   Hist(100, -100,  100, title = "ForthJdPt;4th Leading Jet dPt;Event")
        self.hismap["ForthJdMass"]  =   Hist(100,-50,  50,  title = "ForthJdMass;4th Leading Jet dMass;Event")
        self.hismap["ForthJdEta"]   =   Hist(100, -1, 1,    title = "ForthJdEta;4th Leading Jet dEta;Event")
        self.hismap["ForthJdR"]   =   Hist(60, 0, 6,    title = "ForthJdEta;4th Leading Jet dR;Event")


        ## NJettiness study
        # self.hismap["Tau32"]        = Hist(10,  0, 1,    title = "Tau32;Tau32;No.  of Jets")
        # self.hismap["Tau32Pt200"]   = Hist(10,  0, 1,    title = "Tau32; Tau32 (p_T > 200GeV);No. of Jets")
        # self.hismap["Tau32Pt400"]   = Hist(10,  0, 1,    title = "Tau32; Tau32 (p_T > 400GeV);No. of Jets")
        # self.hismap["Tau21Pt200"]   = Hist(10,  0, 1,    title = "Tau21; Tau21 (p_T > 200GeV);No. of Jets")
        # self.hismap["t_Tau32"]      = Hist(10,  0, 1,    title = "Tau32;Tau32_t;No.  of Jets")
        # self.hismap["t_Tau32Pt400"] = Hist(10,  0, 1,    title = "Tau32; Tau32_t (p_T > 400GeV);No. of Jets")
        # self.hismap["t_Mass"]      = Hist(50, 0, 300, title = "Mass;JetMass_t;No.  of Jets")
        # self.hismap["t_MassPt400"] = Hist(50, 0, 300, title = "Mass; JetMass_t (p_T > 400GeV);No. of Jets")
        # self.hismap["w_Tau21"]      = Hist(10,  0, 1,    title = "Tau21;Tau21_w;No.  of Jets")
        # self.hismap["w_Tau21Pt200"] = Hist(10,  0, 1,    title = "Tau21; Tau21_w (p_T > 200GeV);No. of Jets")
        # self.hismap["w_Mass"]      = Hist(100,  0, 200,    title = "Mass;JetMass_w;No.  of Jets")
        # self.hismap["w_MassPt200"] = Hist(100,  0, 200,    title = "Mass;JetMass_w (p_T > 200GeV);No. of Jets")

    def Loop(self):
        for event in self.tree_:
            self.Run(event)

    def Run(self):
        ## Rate study
        self.Nevent += 1
        self.FillRate()

        e = self.tree_
        if len(e.pt) ==0:
            return False
        self.Jets = []
        self.ptorder = sorted(range(len(self.tree_.pt)), key=lambda k: self.tree_.pt[k])
        lead_idx = self.ptorder[-1]
        self.hismap["FirstJetPt"].Fill(e.pt[lead_idx])
        self.hismap["FirstJetMass"].Fill(e.mass[lead_idx])
        self.hismap["FirstJetEta"].Fill(e.eta[lead_idx])

        self.hismap["SecondJetPt"].Fill(e.pt[self.ptorder[-2]])
        self.hismap["SecondJetMass"].Fill(e.mass[self.ptorder[-2]])
        self.hismap["SecondJetEta"].Fill(e.eta[self.ptorder[-2]])
        self.hismap["ThirdJetPt"].Fill(e.pt[self.ptorder[-3]])
        self.hismap["ThirdJetMass"].Fill(e.mass[self.ptorder[-3]])
        self.hismap["ThirdJetEta"].Fill(e.eta[self.ptorder[-3]])

        self.hismap["ForthJetPt"].Fill(e.pt[self.ptorder[-4]])
        self.hismap["ForthJetMass"].Fill(e.mass[self.ptorder[-4]])
        self.hismap["ForthJetEta"].Fill(e.eta[self.ptorder[-4]])


        self.hismap["NJetPt30"].Fill(sum(1 if x >= 30 else 0 for x in e.pt))
        for i in range(len(e.pt)):
            j = ROOT.TLorentzVector(0, 0, 0, 0)
            j.SetPtEtaPhiM(self.tree_.pt[i], self.tree_.eta[i],
                           self.tree_.phi[i], self.tree_.mass[i])
            self.Jets.append(j)
            # if self.tree_.Tau2[i] != 0:
                # self.hismap["Tau32"].Fill(e.Tau3[i]/ e.Tau2[i])
            # if self.tree_.Tau2[i] != 0:
                # self.hismap["Tau32"].Fill(e.Tau3[i]/ e.Tau2[i])
                # if e.pt[i] > 200:
                    # self.hismap["Tau32Pt200"].Fill(e.Tau3[i]/ e.Tau2[i])
                # if e.pt[i] > 400:
                    # self.hismap["Tau32Pt400"].Fill(e.Tau3[i]/ e.Tau2[i])
            # if self.tree_.Tau1[i] != 0:
                # if e.pt[i] > 200:
                    # self.hismap["Tau21Pt200"].Fill(e.Tau2[i]/ e.Tau1[i])

    def CheckWithGenJet(self, genObjjet, genlepjet):
        topjet = set()
        for i in range(len(self.tree_.pt)):
            j = self.Jets[i]
            for top in genObjjet:
                if top.DeltaR(j) < 0.4:
                    topjet.add(i)
            for lep in genlepjet:
                if lep.DeltaR(j) < 0.4 and i in topjet:
                    topjet.remove(i)
        return topjet

    def PlotTopMatchedJet(self, genTopjet, genlepjet):
        if len(genTopjet) == 0:
            return None
        for i in self.CheckWithGenJet(genTopjet, genlepjet):
            j = self.Jets[i]
            self.hismap["t_Mass"].Fill(j.M())
            self.hismap["t_Tau32"].Fill(self.tree_.Tau3[i]/ self.tree_.Tau2[i]
                                        if self.tree_.Tau2[i]  else 0)
            if j.Pt() > 400:
                self.hismap["t_MassPt400"].Fill(j.M())
                self.hismap["t_Tau32Pt400"].Fill(self.tree_.Tau3[i]/ self.tree_.Tau2[i]
                                        if self.tree_.Tau2[i]  else 0)

    def PlotWMatchedJet(self, genWjet, genlepjet):
        if len(genWjet) == 0:
            return None
        for i in self.CheckWithGenJet(genWjet, genlepjet):
            j = self.Jets[i]
            self.hismap["w_Mass"].Fill(j.M())
            self.hismap["w_Tau21"].Fill(self.tree_.Tau2[i]/ self.tree_.Tau1[i]
                                        if self.tree_.Tau1[i]  else 0)
            if j.Pt() > 200:
                self.hismap["w_MassPt200"].Fill(j.M())
                self.hismap["w_Tau21Pt200"].Fill(self.tree_.Tau2[i]/ self.tree_.Tau1[i]
                                        if self.tree_.Tau1[i]  else 0)

    def FillRate(self):
        self.hismap["NEvent"].Fill(1)
        if len(self.tree_.pt) == 0:
            return None
        lead_jet = max( self.tree_.pt)
        lead_idx = list(self.tree_.pt).index(lead_jet)
        # print(lead_jet, lead_idx, self.tree_.pt[lead_idx])
        for i in range(self.hismap["RatePT"].GetNbinsX()):
            if lead_jet >= self.hismap["RatePT"].GetBinLowEdge(i):
                self.hismap["RatePT"].Fill(self.hismap["RatePT"].GetBinLowEdge(i))
                if abs(self.tree_.eta[lead_idx]) <= 2.4:
                    self.hismap["RatePT_Eta24"].Fill(self.hismap["RatePT_Eta24"].GetBinLowEdge(i))
                if abs(self.tree_.eta[lead_idx]) <= 3:
                    self.hismap["RatePT_Eta3"].Fill(self.hismap["RatePT_Eta3"].GetBinLowEdge(i))


        M40 = [i for i in xrange(len(self.tree_.mass)) if self.tree_.mass[i] > 40] 
        M100 = [i for i in xrange(len(self.tree_.mass)) if self.tree_.mass[i] > 100] 
        M40Eta3 = [i for i in xrange(len(self.tree_.mass)) if
                   self.tree_.mass[i] > 40 and abs(self.tree_.eta[i]) <= 3.0] 

        M40_jet = max([self.tree_.pt[i] for i in M40] or [None])
        M100_jet = max([self.tree_.pt[i] for i in M100] or [None])
        M40Eta3_jet = max([self.tree_.pt[i] for i in M40Eta3] or [None])

        for i in range(self.hismap["RatePT"].GetNbinsX()):
            if M40_jet and M40_jet >= self.hismap["RatePT"].GetBinLowEdge(i):
                self.hismap["RatePT_M40"].Fill(self.hismap["RatePT_M40"].GetBinLowEdge(i))
            if M100_jet and M100_jet >= self.hismap["RatePT"].GetBinLowEdge(i):
                self.hismap["RatePT_M100"].Fill(self.hismap["RatePT_M100"].GetBinLowEdge(i))
            if M40Eta3_jet and M40Eta3_jet >= self.hismap["RatePT"].GetBinLowEdge(i):
                self.hismap["RatePT_Eta3_M40"].Fill(self.hismap["RatePT_Eta3_M40"].GetBinLowEdge(i))


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
        for k, h in self.hismap.items():
            if "Rate" in k:
                h.Scale(11246 * 2808 / self.Nevent/1000)
        return self.hismap

    def GetJet(self, rank):
        idx = self.ptorder[-1*rank]
        j = ROOT.TLorentzVector(0, 0, 0, 0)
        j.SetPtEtaPhiM(self.tree_.pt[idx], self.tree_.eta[idx],  self.tree_.phi[idx], self.tree_.mass[idx])
        return j

    def PlotRelation(self, rank, refJ):
        prefix= ""
        if rank == 1:
            prefix = "FirstJ"
        if rank == 2:
            prefix = "SecondJ"
        if rank == 3:
            prefix = "ThirdJ"
        if rank == 4:
            prefix = "ForthJ"
        idx = self.ptorder[-1*rank]
        j = ROOT.TLorentzVector(0, 0, 0, 0)
        j.SetPtEtaPhiM(self.tree_.pt[idx], self.tree_.eta[idx],  self.tree_.phi[idx], self.tree_.mass[idx])
        self.hismap["%sdPt" % prefix].Fill(j.Pt() - refJ.Pt())
        self.hismap["%sdMass" % prefix].Fill(j.M() - refJ.M())
        self.hismap["%sdEta" % prefix].Fill(j.Eta() - refJ.Eta())
        self.hismap["%sdR" % prefix].Fill(j.DeltaR(refJ))
