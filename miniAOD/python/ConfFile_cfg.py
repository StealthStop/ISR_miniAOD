import FWCore.ParameterSet.Config as cms

process = cms.Process("miniAOD")

# -----------------
# get the cfi files
# -----------------
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi") 

# ----------------
# number of events
# ---------------- 
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) ) # -1 for all events

# -----------------
# get the root file
# -----------------
process.source = cms.Source("PoolSource",
                                # --------------------------------------------------------------
                                # ------------------------ 2016 RPV_350 ------------------------
                                # --------------------------------------------------------------  
                                #filesNames = (cms.untracked.vstring('root://cmsxrootd.fnal.gov//store/mc/RunIISummer16MiniAODv3/RPV_2t6j_mStop-350_mN1-100_TuneCUEP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v1/10000/28D1E302-C3F4-E811-8546-901B0E542756.root')
                                #fileNames = cms.untracked.vstring('file:/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/ISR_miniAOD/CMSSW_10_2_9/src/ISR_miniAOD/miniAOD/python/28D1E302-C3F4-E811-8546-901B0E542756.root')

                                # --------------------------------------------------------------
                                # ------------------------ 2016 RPV_550 ------------------------
                                # -------------------------------------------------------------- 
                                #fileNames = cms.untracked.vstring('root://cmsxrootd.fnal.gov//store/mc/RunIISummer16MiniAODv3/RPV_2t6j_mStop-550_mN1-100_TuneCUEP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v1/90000/16D5849D-50F3-E811-89F4-0CC47AFC3C86.root')
                                #fileNames = cms.untracked.vstring('file:/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/ISR_miniAOD/CMSSW_10_2_9/src/ISR_miniAOD/miniAOD/python/16D5849D-50F3-E811-89F4-0CC47AFC3C86.root')

                                # --------------------------------------------------------------
                                # ------------------------ 2016 RPV_850 ------------------------
                                # -------------------------------------------------------------- 
                                #fileNames = cms.untracked.vstring('root://cmsxrootd.fnal.gov//store/mc/RunIISummer16MiniAODv3/RPV_2t6j_mStop-850_mN1-100_TuneCUEP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v1/20000/2646E8AE-4DF1-E811-A7D4-549F358EB7D7.root')
                                fileNames = cms.untracked.vstring('file:/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/ISR_miniAOD/CMSSW_10_2_9/src/ISR_miniAOD/miniAOD/python/2646E8AE-4DF1-E811-A7D4-549F358EB7D7.root')
                                
                                # -------------------------------------------------------------
                                # ------------------------ 2016 TTJets ------------------------
                                # ------------------------------------------------------------- 
                                #fileNames = cms.untracked.vstring('root://cmsxrootd.fnal.gov//store/mc/RunIISummer16MiniAODv3/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v1/40000/0A4EAAB1-9223-E911-B512-A4BF01283A8B.root')
                                #fileNames = cms.untracked.vstring('file:/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/ISR_miniAOD/CMSSW_10_2_9/src/ISR_miniAOD/miniAOD/python/0A4EAAB1-9223-E911-B512-A4BF01283A8B.root')
                            
                                # ---------------------------------------------------------    
                                # ------------------------ 2016 TT ------------------------
                                # --------------------------------------------------------- 
                                #fileNames = cms.untracked.vstring('root://cmsxrootd.fnal.gov//store/mc/RunIISummer16MiniAODv3/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_backup_94X_mcRun2_asymptotic_v3-v2/00000/000708C6-53F1-E811-A279-0CC47A4D75EE.root')
                                #fileNames = cms.untracked.vstring('file:/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/ISR_miniAOD/CMSSW_10_2_9/src/ISR_miniAOD/miniAOD/python/000708C6-53F1-E811-A279-0CC47A4D75EE.root')
                                
                           )

# ---------------------
# get the gen particles
# ---------------------
process.demo = cms.EDAnalyzer('ISR_miniAOD',
                                genParticles = cms.untracked.InputTag("prunedGenParticles"))


# -------------------
# print particle list
# -------------------
#process.printAllTree = cms.EDAnalyzer("ParticleListDrawer",
#                                    maxEventsToPrint         = cms.untracked.int32(10), # 1
#                                    printVertex              = cms.untracked.bool(False),
#                                    printOnlyHardInteraction = cms.untracked.bool(False),
#                                    src                      = cms.InputTag("prunedGenParticles")
#                                  )
#process.p = cms.Path(process.printAllTree)

# -----------
# print decay
# -----------
#process.printDecay = cms.EDAnalyzer("ParticleDecayDrawer",
#                                    src           = cms.InputTag("prunedGenParticles"),
#                                    printP4       = cms.untracked.bool(False),
#                                    printPtEtaPhi = cms.untracked.bool(False),
#                                    printVertex   = cms.untracked.bool(False)
#                                   )
#process.p = cms.Path(process.printDecay)

# ----------
# print tree
# ----------
#process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
#                                   src           = cms.InputTag("prunedGenParticles"),                       
#                                  #printP4       = cms.untracked.bool(False),
#                                  #printPtEtaPhi = cms.untracked.bool(False),
#                                  #printVertex   = cms.untracked.bool(False),
#                                  #printStatus   = cms.untracked.bool(False),
#                                  #printIndex    = cms.untracked.bool(False),
#                                  #status        = cms.untracked.vint32( 3 )
#                                  )
#process.p = cms.Path(process.printTree)

#process.p = cms.Path(process.printAllTree
#                     *process.printDecay
#                     * process.printTree
#                    )

# ---------------------------
# create the output root file
# ---------------------------
#process.out = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string("output.root"))
process.TFileService = cms.Service("TFileService", fileName = cms.string("2016_miniAOD_RPV850.root") )
#process.TFileService = cms.Service("TFileService", fileName = cms.string("TreeMaker.root") )

process.p  = cms.Path(process.demo)
#process.ep = cms.EndPath(process.out)

process.schedule = cms.Schedule( process.p )


