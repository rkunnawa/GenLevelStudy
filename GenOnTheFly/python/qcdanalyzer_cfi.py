import FWCore.ParameterSet.Config as cms

qcdAna = cms.EDAnalyzer('QCDAnalyzer',
                        genJetSrc = cms.InputTag("ak3GenJets"),
                        genParticleSrc = cms.InputTag("genParticles"),
                        doFlavor = cms.bool(False),
                        doSpecies = cms.bool(False),
                        onlyDS = cms.bool(False),
                        flavorSrc = cms.InputTag("flavourByRef"),
                        flavorId = cms.vint32(21), # gluons
                        speciesId = cms.vint32(11), # pi0
                        useRapidity = cms.bool(False),
                        jetEtaMin = cms.double(-1.0),
                        jetEtaMax = cms.double(1.0),
                        hEtaMin = cms.double(-0.7),
                        hEtaMax = cms.double(0.7),
                        jetRadius = cms.double(0.3),
                        pthatMin = cms.double(0.),
                        pthatMax = cms.double(20.),
                        jetPtBins = cms.vdouble( 3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 
                                                 47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 
                                                 196, 220, 245, 272, 300, 429, 692, 1000 ),
                        jetPhiBins = cms.vdouble(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.15),
                        hPtBins = cms.vdouble( 3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 
                                               47, 55, 64, 74, 84, 97, 114, 133, 153, 174, 
                                               196, 220, 245, 272, 300, 429, 692, 1000 ),
                        qScalePtBins = cms.vdouble( 3, 4, 5, 7, 9, 12, 15, 18, 22, 27, 33, 39, 
                                                    47, 55, 64, 74, 84, 97, 114, 133, 153, 174,
                                                    196, 220, 245, 272, 300, 429, 692, 1000 ),
                        etaBins = cms.vdouble( -2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5 ) 
)
