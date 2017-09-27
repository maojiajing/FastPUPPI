import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing ('analysis')

# setup any defaults you want
options.outputFile = 'test.root'
options.inputFiles= ''
options.maxEvents = -1 # -1 means all events

# get and parse the command line arguments
options.parseArguments()

process = cms.Process("OWNPARTICLES")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

ifiles = []
with open(options.inputFiles[0]) as files:
    for j in files.readlines():
        ifiles.append(j.strip())

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        # 'file:inputs_17D.root'
        ifiles
    )
)
process.TFileService = cms.Service("TFileService", fileName = cms.string(options.outputFile))

process.source.duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
process.load('FastPUPPI.NtupleProducer.l1tPFHGCalProducerFrom3DTPsEM_cfi')
process.load('FastPUPPI.NtupleProducer.caloNtupleProducer_cfi')
process.load('FastPUPPI.NtupleProducer.ntupleProducer_cfi')
process.CaloInfoOut.outputName = ""; # turn off Ntuples
process.InfoOut.outputName = ""; # turn off Ntuples
process.pp = cms.Sequence(process.l1tPFHGCalProducerFrom3DTPsEM + process.CaloInfoOut + process.InfoOut)

# InputSrc = [ 'InfoOut:Puppi']
InputSrc = [ 'InfoOut:RawCalo', 'InfoOut:Calo', 'InfoOut:TK', 'InfoOut:TKVtx', 'InfoOut:PF', 'InfoOut:Puppi' ]

#============================================================================#
#-------------------------------     Jets     -------------------------------#
#============================================================================#
from RecoJets.Configuration.RecoPFJets_cff import ak8PFJetsCHSSoftDrop
from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness

def AddJetCollection(cone, inputs):
    names = []
    ## Default jets
    attrname = "ak%d%s"  % ( cone *10, inputs.split(":")[-1])
    names.append(attrname)
    setattr(process, attrname,   ak4PFJets.clone(rParam       = cms.double(cone),  src = inputs) )
    # ## SoftDrop
    # attrname = "ak%d%sSD"  % ( cone *10, inputs.split(":")[-1])
    # names.append(attrname)
    # setattr(process, attrname,   ak8PFJetsCHSSoftDrop.clone(rParam       = cms.double(cone),  src = inputs) )
    ## NJetti
    attrname = "ak%d%sNJ"  % ( cone *10, inputs.split(":")[-1])
    names.append(attrname)
    setattr(process, attrname,   Njettiness.clone(src = "ak%d%s"  % ( cone *10, inputs.split(":")[-1])))
    return names


attrnames = []
for j in InputSrc:
    attrnames += AddJetCollection(0.8, j )

strf = "+".join(["process.%s" % attrname for attrname in attrnames])
process.jets =cms.Sequence(eval(strf))
#============================================================================#
#-----------------------------     Producer     -----------------------------#
#============================================================================#
pp = []
strf = ""
for j in attrnames:
    if "NJ" in j:
        continue
    attrname = "my" + j
    setattr(process, attrname,   cms.EDAnalyzer('FatJet',
                                                FatJetTag = cms.InputTag(j),
                                                Tau1Tag = cms.InputTag (j+"NJ", "tau1"),
                                                Tau2Tag = cms.InputTag (j+"NJ", "tau2"),
                                                Tau3Tag = cms.InputTag (j+"NJ", "tau3"),
                                                Tau4Tag = cms.InputTag (j+"NJ", "tau4"),
                                                ))
    pp.append(getattr(process, attrname, None))
    pp.append(getattr(process, attrname, None))
    strf += "process.%s + " % attrname
process.ppp =cms.Sequence(eval(strf[:-3]))


from RecoJets.Configuration.RecoGenJets_cff import ak8GenJetsNoNu
from RecoJets.Configuration.GenJetParticles_cff import genParticlesForJetsNoNu
process.genParticlesForJetsNoNu2 = genParticlesForJetsNoNu.clone(src = "genParticles")
process.ak8GenJetsNoNu2 = ak8GenJetsNoNu.clone(src="genParticlesForJetsNoNu2")

process.ntuple = cms.EDAnalyzer("ResponseNTuplizer",
    genJets = cms.InputTag("ak8GenJetsNoNu2"),
    genParticles = cms.InputTag("genParticles"),
    isParticleGun = cms.bool(False),
    doRandom = cms.bool(False),
    objects = cms.PSet(),
    copyUInts = cms.VInputTag()
)

process.p = cms.Path(process.pp + process.jets + process.ppp+
                     process.genParticlesForJetsNoNu2 + process.ak8GenJetsNoNu2+process.ntuple)


# process.out = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string('myOutputFile.root'))
# process.e = cms.EndPath(process.pp+process.out)
