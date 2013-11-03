import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

InputFileName = "results/DYToMuMu.root"
OutputFilePrefix = "efficiency-"

EfficiencyBinningSpecification = cms.PSet(  UnbinnedVariables = cms.vstring('mass'),
                                            BinnedVariables = cms.PSet( pt=cms.vdouble(20,40,200),
                                                                        eta=cms.vdouble(-2.5,-1.5,0.0,1.5,2.5)
                                                                        ),
                                            BinToPDFmap = cms.vstring('pdfSplusB')
                                            )


####
# Muon -> id+iso efficiency
####
process.TagProbeFitTreeAnalyzer = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
                                                 InputFileNames = cms.vstring(InputFileName),
                                                 InputDirectoryName = cms.string("tandp"),
                                                 InputTreeName = cms.string("tandp"),
                                                 OutputFileName = cms.string(OutputFilePrefix+"MuonSelection.root"),
                                                 NumCPU = cms.uint32(1),
                                                 SaveWorkspace = cms.bool(True),
                                                 floatShapeParameters = cms.bool(True),
                                                 Variables = cms.PSet( mass = cms.vstring("Mass", "60.0", "120.0", "[GeV]"),
                                                                       pt   = cms.vstring("Transverse momentum", "0", "1000", "[GeV]"),
                                                                       eta  = cms.vstring("Pseudo-rapidity", "-2.5", "2.5", ""),                
                                                                       ),
                                                 Categories = cms.PSet( passLoose    = cms.vstring("passLoose",    "dummy[pass=1,fail=0]"),
                                                                        passLooseIso = cms.vstring("passLooseIso", "dummy[pass=1,fail=0]"),
                                                                        passTight    = cms.vstring("passTight",    "dummy[pass=1,fail=0]"),
                                                                        passTightIso = cms.vstring("passTightIso", "dummy[pass=1,fail=0]"),
                                                                        ),
                                                 PDFs = cms.PSet( pdfSplusB = cms.vstring( 'Gaussian::signal(mass, mean[91.0,90.0,92.0], sigma[5,0,20])',
                                                                                           'Chebychev::backgroundPass(mass, cPass[0,-1,1])',
                                                                                           'Chebychev::backgroundFail(mass, cFail[0,-1,1])',
                                                                                           'efficiency[0.8,0,1]',
                                                                                           'signalFractionInPassing[1.0]'     
                                                                                           ),
                                                                  ),
                                                 Efficiencies = cms.PSet(  Loose = cms.PSet( EfficiencyBinningSpecification,
                                                                                             EfficiencyCategoryAndState = cms.vstring("passLoose","pass","passLooseIso","pass")
                                                                                             ),
                                                                           Tight = cms.PSet( EfficiencyBinningSpecification,
                                                                                             EfficiencyCategoryAndState = cms.vstring("passTight","pass","passTightIso","pass")
                                                                                             )
                                                                           )                              
                                                 )

# run the analyizer 
process.p = cms.Path( process.TagProbeFitTreeAnalyzer )
