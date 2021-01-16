// ----------------------------------------------------------------------------- //
// Package:    Demo/DemoAnalyzer                                                 //
// Class:      DemoAnalyzer                                                      //
//                                                                               //
// \class DemoAnalyzer DemoAnalyzer.cc Demo/DemoAnalyzer/plugins/DemoAnalyzer.cc //
//                                                                               //    
// Original Author:  Semra Turkcapar                                             //
//         Created:  Mon, 14 Sep 2020 16:50:48 GMT                               //
// ----------------------------------------------------------------------------- //

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

// from TreeMaker/Utils/src/ISRJetProducer.cc 
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "TFile.h"
#include "TTree.h"

// --------------------
// class declaration
// --------------------
class ISR_miniAOD : public edm::one::EDAnalyzer<edm::one::SharedResources>  
{
    public:
        explicit ISR_miniAOD(const edm::ParameterSet&);
        ~ISR_miniAOD();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        // -----------
        // Member data
        // ----------- 
        edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;        
   
        // --------------------------------------- 
        // ISR definition - dau-mom: qg, gq, gg   
        // ---------------------------------------  
        TH1F * miniAOD_ISR_Eta; 

        // -------------------------------
        // Intermediate ISR definitions 
        // -------------------------------
        // dau-mom: gluon - udsc 
        TH1F * miniAOD_ISR_Eta_gudsc;

        // dau-mom: gluon - b
        TH1F * miniAOD_ISR_Eta_gb;

        // dau-mom: gluon - t
        TH1F * miniAOD_ISR_Eta_gt;
        
        // dau-mom: udsc - gluon
        TH1F * miniAOD_ISR_Eta_udscg;
        
        // dau-mom: b - gluon
        TH1F * miniAOD_ISR_Eta_bg;
        
        // dau-mom: t - gluon 
        TH1F * miniAOD_ISR_Eta_tg;

        // dau-mom : gluon -> gluon
        TH1F * miniAOD_ISR_Eta_gg;
};

// ------------------------------
// constructors and destructor
// ------------------------------
ISR_miniAOD::ISR_miniAOD(const edm::ParameterSet& iConfig) :
    genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("genParticles")))

{
    edm::Service<TFileService> fs;

    // ---------------------------------------
    // ISR definition - dau-mom: qg, gq, gg   
    // ---------------------------------------
    miniAOD_ISR_Eta = fs->make<TH1F>( "miniAOD_ISR_Eta", "ISR Eta", 100, -6, 6  );
    
    // -------------------------------
    // Intermediate ISR definitions   
    // ------------------------------- 
    // dau-mom: gluon - udsc 
    miniAOD_ISR_Eta_gudsc = fs->make<TH1F>( "miniAOD_ISR_Eta_gudsc", "ISR Eta", 100, -6, 6  );
    
    // dau-mom: gluon - b
    miniAOD_ISR_Eta_gb = fs->make<TH1F>( "miniAOD_ISR_Eta_gb", "ISR Eta", 100, -6, 6  );
    
    // dau-mom: gluon - t
    miniAOD_ISR_Eta_gt = fs->make<TH1F>( "miniAOD_ISR_Eta_gt", "ISR Eta", 100, -6, 6  );
    
    // dau-mom: udsc - gluon
    miniAOD_ISR_Eta_udscg = fs->make<TH1F>( "miniAOD_ISR_Eta_udscg", "ISR Eta", 100, -6, 6  );
    
    // dau-mom: b - gluon
    miniAOD_ISR_Eta_bg = fs->make<TH1F>( "miniAOD_ISR_Eta_bg", "ISR Eta", 100, -6, 6  );
    
    // dau-mom: t - gluon 
    miniAOD_ISR_Eta_tg = fs->make<TH1F>( "miniAOD_ISR_Eta_tg", "ISR Eta", 100, -6, 6  );
    
    // dau-mom : gluon -> gluon
    miniAOD_ISR_Eta_gg = fs->make<TH1F>( "miniAOD_ISR_Eta_gg", "ISR Eta", 100, -6, 6  );

}

ISR_miniAOD::~ISR_miniAOD()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

// -------------------------------
// Member functions
// method called for each event 
// -------------------------------
void ISR_miniAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;

    // ------------
    // GenParticles
    // ------------
    Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByToken(genParticlesToken_, genParticles);   
    for (reco::GenParticleCollection::const_iterator igGen = genParticles->begin();
        igGen != genParticles->end();
        ++igGen)
    {
        const reco::GenParticle &p = (*igGen);
        int grandDauNum = p.numberOfDaughters();
        int dauPdgId    = p.pdgId();
        int momNum      = p.numberOfMothers();
        int status      = p.status();

        if (momNum == 0) continue;        
        
        const reco::Candidate *mom = p.mother();
        int momPdgId  = mom->pdgId();
        int dauNum    = mom->numberOfDaughters();
        double momPt  = mom->pt();
        double dauPt  = p.p4().pt();
        double dauEta = p.p4().eta();

        // ---------------------------------------
        // ISR definition - dau-mom: qg, gq, gg 
        // ---------------------------------------
        if ( ( abs(dauPdgId) <= 6  && abs(momPdgId) == 21 && status == 23 ) ||
           (   abs(dauPdgId) == 21 && abs(momPdgId) <= 6  && status == 23 ) || 
           (   abs(dauPdgId) == 21 && abs(momPdgId) == 21 && status == 23 ) ) 
        {
            miniAOD_ISR_Eta->Fill(dauEta); 
        }
       
        // -------------------------------
        // Intermediate ISR definitions   
        // -------------------------------
        // dau-mom: gluon - udsc 
        if ( abs(dauPdgId) == 21  && abs(momPdgId) <= 4 && status == 23 )
        {
            miniAOD_ISR_Eta_gudsc->Fill(dauEta);   
        }

        // dau-mom: gluon - b
        if ( abs(dauPdgId) == 21  && abs(momPdgId) == 5 && status == 23 )
        {
            //std::cout << "dau: " << dauPdgId << " || " << "mom: " << momPdgId << std::endl;
            miniAOD_ISR_Eta_gb->Fill(dauEta);
        }

        // dau-mom: gluon - t
        if ( abs(dauPdgId) == 21  && abs(momPdgId) == 6 && status == 23 )
        {
            //std::cout << "dau: " << dauPdgId << " || " << "mom: " << momPdgId << std::endl;
            miniAOD_ISR_Eta_gt->Fill(dauEta);
        }
    
        // dau-mom: udsc - gluon 
        if ( abs(dauPdgId) <= 4  && abs(momPdgId) == 21 && status == 23 )
        {
            miniAOD_ISR_Eta_udscg->Fill(dauEta);
        }

        // dau-mom: b - gluon
        if ( abs(dauPdgId) == 5  && abs(momPdgId) == 21 && status == 23 )
        {
            miniAOD_ISR_Eta_bg->Fill(dauEta);
        }

        // dau-mom: t - gluon
        if ( abs(dauPdgId) == 6  && abs(momPdgId) == 21 && status == 23 )
        {
            //std::cout << "dau: " << dauPdgId << " || " << "mom: " << momPdgId << std::endl;
            miniAOD_ISR_Eta_tg->Fill(dauEta);
        }

        // dau-mom: gluon - gluon
        if ( abs(dauPdgId) == 21  && abs(momPdgId) == 21 && status == 23 )
        {
            miniAOD_ISR_Eta_gg->Fill(dauEta);
        }



    } // gen particles loop
    
} // event loop

// --------------------------------------------------------------
// Method called once each job just before starting event loop  
// --------------------------------------------------------------
void ISR_miniAOD::beginJob()
{
}

// ---------------------------------------------------------------
// Method called once each job just after ending the event loop  
// ---------------------------------------------------------------
void ISR_miniAOD::endJob()
{
}

// -------------------------------------------------------------------------
// Method fills 'descriptions' with the allowed parameters for the module 
// -------------------------------------------------------------------------
void ISR_miniAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
    // The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);

    // Specify that only 'tracks' is allowed
    // To use, remove the default given above and uncomment below
    //ParameterSetDescription desc;
    //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
    //descriptions.addDefault(desc);
}

// ---------------------------
// Define this as a plug-in
// ---------------------------
DEFINE_FWK_MODULE(ISR_miniAOD);




