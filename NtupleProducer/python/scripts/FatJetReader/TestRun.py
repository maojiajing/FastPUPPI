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
from rootpy.io import root_open
from rootpy.tree import Tree,TreeChain
from rootpy.interactive import wait
from rootpy.plotting import Hist,Canvas
from rootpy.plotting.style import get_style
from AK8JetReader import AK8JetReader
from GenJetReader import GenJetReader
from collections import defaultdict
from Ploter import *
import rootpy.stl as stl
import rootpy
import argparse
import ROOT
import pprint
import glob
import pickle


def GetTree(f, treename):
    folder = getattr(f, treename, None)
    if folder is not None:
        return folder.tree

def Test():
    tt = Tree()
    tt.glob()

def LoadTree():
    global treemap
    global ak8map
    global genmap
    for k, fname in filename.items():
        t = {}
        ak8map[k] = {}
        for g, v in treename.items():
            t[g] = ROOT.TChain('%s/tree' % v)
            for j in glob.glob(fname):
                t[g].Add(j)
            if g != 'gen':
                ak8map[k][g] = AK8JetReader(t[g], v)
            if g == 'gen':
                genmap[k] = GenJetReader(t[g], v)
        treemap[k] = t


if __name__ == "__main__":
    folder = "/Users/benwu/Data/Dataset/L1PFInputs/FatJetNtuple/Oct02withGen_v1/"
    filename = {
        "QCD_PU0"     : "%s/QCD_PU0*.root"       % folder,
        "QCD_PU140"   : "%s/QCD_PU140*.root"     % folder,
        "TTbar_PU140" : "%s/TTbar_PU140*.root"   % folder,
        "TTbar_PU0"   : "%s/TTbar_PU0*.root"     % folder,
    }

    treename     = {
        'gen'        : 'GenJet',
        'ak8TKVtx'   : 'myak8TKVtx',
        'ak8TK'      : 'myak8TK',
        'ak8Puppi'   : 'myak8Puppi',
        'ak8PF'      : 'myak8PF',
        'ak8Calo'    : 'myak8Calo',
        'ak8RawCalo' : 'myak8RawCalo',
    }

    # Create a map of TTree
    treemap = defaultdict(dict)
    ak8map = defaultdict(dict)
    genmap = defaultdict(dict)
    LoadTree()

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('process', nargs='*', default='TTbar_PU0', help='Select the process for run')
    args = parser.parse_args()

    processes = args.process

    for process in processes:
        for i in xrange(treemap[process].values()[0].GetEntries()):
            # if i > 30:
                # break
            ## No Accessing to Tchain friend in PyROOT :-(
            [ tr.GetEntry(i) for tr in  treemap[process].values() ]
            genmap[process].Run()
            genTops = genmap[process].GetGenJets([6])
            genWs = genmap[process].GetGenJets([24])
            genLeps = genmap[process].GetGenJets([11, 13, 15])
            for tr in ak8map[process].keys():
                ak8map[process][tr].Run()
                ak8map[process][tr].PlotTopMatchedJet(genTops, genLeps)
                ak8map[process][tr].PlotWMatchedJet(genWs, genLeps)



                # ak8map[process][tr].Draw()

    # Getting back the historgram
    histmap = defaultdict(dict)
    for process in processes:
            for tr in  ak8map[process].keys():
                histmap[process][tr] = ak8map[process][tr].GetHist()

    pprint.pprint(histmap)
    pickle.dump(histmap, open("%s.p" % args.process[0], "wb"))
