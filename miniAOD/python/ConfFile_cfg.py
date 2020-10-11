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
                                fileNames = cms.untracked.vstring('root://cmsxrootd.fnal.gov//store/mc/RunIISummer16MiniAODv3/RPV_2t6j_mStop-550_mN1-100_TuneCUEP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v1/90000/16D5849D-50F3-E811-89F4-0CC47AFC3C86.root')
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
process.TFileService = cms.Service("TFileService", fileName = cms.string("miniAOD.root") )


process.p  = cms.Path(process.demo)
#process.ep = cms.EndPath(process.out)

process.schedule = cms.Schedule( process.p )


