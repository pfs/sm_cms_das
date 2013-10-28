import FWCore.ParameterSet.Config as cms

process = cms.Process("DataAna")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = gtag

process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True),
                                        SkipEvent = cms.untracked.vstring('ProductNotFound')
                                        ) 


process.source = cms.Source("PoolSource",
                            inputList = cms.untracked.vstring('/store/data//Run2012A/DoubleMu/AOD//22Jan2013-v1/20000/F4C34C30-B581-E211-8269-003048FFD7A2.root') 
                            )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

if(isMC) : process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.TFileService = cms.Service("TFileService", fileName = cms.string('DataAnalysis.root'))



##-------------------- Import the JEC services -----------------------
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

##-------------------- Import the Jet RECO modules -----------------------
process.load('RecoJets.Configuration.RecoPFJets_cff')

process.kt6PFJets.doRhoFastjet = True
process.ak5PFJets.doAreaFastjet = True

#apply a good vertex selector and filter out scraping
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter( "PrimaryVertexObjectFilter",
                                                   filterParams = pvSelector.clone( minNdof = cms.double(4.0),
                                                                                    maxZ = cms.double(24.0),
                                                                                    maxd0 = cms.double(2.0)
                                                                                    ),
                                                   src=cms.InputTag('offlinePrimaryVertices')
                                                   )
process.goodVertexFilter = cms.EDFilter("GoodVertexFilter",
                                        vertexCollection = cms.InputTag('goodOfflinePrimaryVertices'),
                                        minimumNDOF = cms.uint32(4),
                                        maxAbsZ = cms.double(24),
                                        maxd0 = cms.double(2)
                                        )

process.noscraping = cms.EDFilter("FilterOutScraping",
                                  applyfilter = cms.untracked.bool(True),
                                  debugOn = cms.untracked.bool(False),
                                  numtrack = cms.untracked.uint32(10),
                                  thresh = cms.untracked.double(0.25)
                               )

# optional MET filters
# cf.https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFilters
process.load('RecoMET.METFilters.metFilters_cff')
process.hcalLaserEventFilter.taggingMode   = cms.bool(True)
process.EcalDeadCellTriggerPrimitiveFilter.taggingMode=cms.bool(True)
process.eeBadScFilter.taggingMode           = cms.bool(True)
process.ecalLaserCorrFilter.taggingMode     = cms.bool(True)
process.trackingFailureFilter.VertexSource  = cms.InputTag('goodOfflinePrimaryVertices')
process.trackingFailureFilter.taggingMode   = cms.bool(True)
process.manystripclus53X.taggedMode         = cms.untracked.bool(True)
process.manystripclus53X.forcedValue        = cms.untracked.bool(False)
process.toomanystripclus53X.taggedMode      = cms.untracked.bool(True)
process.toomanystripclus53X.forcedValue     = cms.untracked.bool(False)
process.logErrorTooManyClusters.taggedMode  = cms.untracked.bool(True)
process.logErrorTooManyClusters.forcedValue = cms.untracked.bool(False)  

process.metFilteringTaggers = cms.Sequence(process.HBHENoiseFilter*
                                           process.hcalLaserEventFilter *
                                           process.EcalDeadCellTriggerPrimitiveFilter *
                                           process.eeBadScFilter *
                                           process.ecalLaserCorrFilter *
                                           process.trackingFailureFilter *
                                           process.trkPOGFilters)

#PF2PAT
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.tools.pfTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *

postfix = "PFlow"
jetAlgo="AK5"
jecLevels=['L1FastJet', 'L2Relative', 'L3Absolute']
if(not isMC): jecLevels.append('L2L3Residual')

usePF2PAT(process,
          runPF2PAT=True,
          jetAlgo=jetAlgo,
          runOnMC=isMC,
          postfix=postfix,
          jetCorrections=('AK5PFchs', jecLevels),
          pvCollection=cms.InputTag('goodOfflinePrimaryVertices'),
          typeIMetCorrections=False)

#setup trigger matching
#from UserCode.llvv_fwk.triggerMatching_cfg import *
#addTriggerMatchingTo(process)

#custom electrons
useGsfElectrons(process,postfix=postfix,dR="03")
process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi')
process.eidMVASequence = cms.Sequence(  process.mvaTrigV0 + process.mvaNonTrigV0 )
process.patElectronsPFlow.electronIDSources.mvaTrigV0    = cms.InputTag("mvaTrigV0")
process.patElectronsPFlow.electronIDSources.mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0")

#custom muons
process.patMuonsPFlow.pfMuonSource = cms.InputTag("pfSelectedMuonsPFlow")
process.muonMatchPFlow.src = cms.InputTag("pfSelectedMuonsPFlow")

#custom jets for CHS
process.pfPileUpPFlow.checkClosestZVertex = cms.bool(False)
process.pfPileUpIsoPFlow.checkClosestZVertex = cms.bool(False)
getattr(process,"pfNoPileUp"+postfix).enable = True
getattr(process,"pfNoMuon"+postfix).enable = False     # to use muon-cleaned electron collection set to True (check iso)
getattr(process,"pfNoElectron"+postfix).enable = False # to use electron-cleaned tau collection set to True (check iso)
getattr(process,"pfNoTau"+postfix).enable = False      # to use tau-cleaned jet collection set to True (check what is a tau)
getattr(process,"pfNoJet"+postfix).enable = True       # this i guess it's for photons...      

#compute rho from central pf candidates only
from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
process.kt6PFJetsCentral = kt4PFJets.clone( rParam = cms.double(0.6),
                                            doAreaFastjet = cms.bool(True),
                                            doRhoFastjet = cms.bool(True),
                                            Rho_EtaMax = cms.double(2.5),
                                            Ghost_EtaMax = cms.double(2.5) )

# cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMetAnalysis
process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
process.load("JetMETCorrections.Type1MET.pfMETCorrectionType0_cfi")
process.pfType1CorrectedMet.applyType0Corrections = cms.bool(False)
process.pfType1CorrectedMet.srcType1Corrections = cms.VInputTag( cms.InputTag('pfMETcorrType0'),
                                                                 cms.InputTag('pfJetMETcorr', 'type1')
                                                                 )

process.dataAnalyzer = cms.EDAnalyzer( "DataAnalyzer",
                                       triggerSource = cms.InputTag("TriggerResults::HLT"),
                                       triggerPaths = cms.vstring('HLT_Mu17_Mu8_v',
                                                                  'HLT_Mu17_TkMu8_v',
                                                                  'HLT_IsoMu24_eta2p1_v'),
                                       triggerCats  = cms.vint32(1313,
                                                                 1313,
                                                                 13),
                                       genSource       = cms.InputTag("genParticles"),
                                       vtxSource       = cms.InputTag("goodOfflinePrimaryVertices"),
                                       beamSpotSource  = cms.InputTag("offlineBeamSpot"),
                                       muonSource      = cms.InputTag("selectedPatMuonsTriggerMatch"),
                                       electronSource  = cms.InputTag("selectedPatElectrons"),
                                       conversionSource= cms.InputTag("allConversions"),
                                       ebrechitsSource = cms.InputTag("reducedEcalRecHitsEB"),
                                       eerechitsSource = cms.InputTag("reducedEcalRecHitsEE"),
                                       jetSource       = cms.InputTag("selectedPatJetsPFlow"),
                                       metSource       = cms.VInputTag("pfMETPFlow","pfMet","pfType1CorrectedMet","pfType1p2CorrectedMet"),
                                       rhoSource       = cms.InputTag("kt6PFJets:rho")
                                       )
                
#counters for specific filters
process.startCounter = cms.EDProducer("EventCountProducer")
process.scrapCounter = process.startCounter.clone()
process.vtxCounter   = process.startCounter.clone()
process.metCounter   = process.startCounter.clone() 
process.p = cms.Path( process.startCounter
                      *process.noscraping
                      *process.scrapCounter
                      *process.goodOfflinePrimaryVertices
                      *process.goodVertexFilter
                      *process.vtxCounter
                      *process.metFilteringTaggers
		      *process.metCounter
                      *process.eidMVASequence
                      *getattr(process,"patPF2PATSequence"+postfix)
                      *process.kt6PFJetsCentral
                      *process.type0PFMEtCorrection*process.producePFMETCorrections
                      *process.selectedPatElectronsWithTrigger*process.selectedPatElectronsPFlowHeep
                      *process.selectedPatMuonsTriggerMatch 
                      *process.dataAnalyzer
                      )

	




