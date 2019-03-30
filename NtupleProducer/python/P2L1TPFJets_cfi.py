import FWCore.ParameterSet.Config as cms

SeedJet  =  cms.EDProducer('SeedJetProducer',
                           src = cms.InputTag(""),
                           seedPtCut    = cms.double(3.0),
                           cone    = cms.double(0.4),
                           seedID    = cms.uint32(999),
                          )

P2L1TJet  =  cms.EDProducer('P2L1TPFJetProducer',
                            src          = cms.InputTag('particleFlow'),
                            rParam       = cms.double(0.4),
                            groupR       = cms.double(0.4),
                            jetPtMin     = cms.double(5.0),
                            inputEtMin   = cms.double(0.0),
                            ktsign = cms.int32(-1),
                            metricKt = cms.bool(False),
                            metricdR = cms.bool(False),
                            mergingE = cms.bool(False),
                            mergingWTA = cms.bool(False),
                            resumWTA = cms.bool(False),
                            N2Tile = cms.bool(False),
                            N2Group = cms.bool(False),
                          )
