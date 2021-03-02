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
        // ISR definition - mom-dau: gg, qg, gq  
        // --------------------------------------- 
        // before debug - first version - Semra
        TH1F * miniAOD_ISR_Eta;
        TH1F * miniAOD_ISR_Eta_gg;
        TH1F * miniAOD_ISR_Eta_qg;
        TH1F * miniAOD_ISR_Eta_gq;
        //TH1F * miniAOD_ISR_Pt;
        //TH1F * miniAOD_ISR_Pt_gg;
        //TH1F * miniAOD_ISR_Pt_qg;
        //TH1F * miniAOD_ISR_Pt_gq;
        
        // after debug - last version - Josh
        TH1F * miniAOD_ISR_Eta_g;
        TH1F * miniAOD_ISR_Eta_gg_g;
        TH1F * miniAOD_ISR_Eta_gg_q;
        TH1F * miniAOD_ISR_Eta_qq_g;
        TH1F * miniAOD_ISR_Eta_qq_q;
        TH1F * miniAOD_ISR_Eta_gq_gq;
        //TH1F * miniAOD_ISR_Pt_g;
        //TH1F * miniAOD_ISR_Pt_gg_g;
        //TH1F * miniAOD_ISR_Pt_gg_q;
        //TH1F * miniAOD_ISR_Pt_qq_g;
        //TH1F * miniAOD_ISR_Pt_qq_q;
        //TH1F * miniAOD_ISR_Pt_gq_gq;

        // events counters
        TH1F * totalEvents;
        TH1F * NumISR_Events;
        TH1F * NumISR_Events_gg;
        TH1F * NumISR_Events_qg;
        TH1F * NumISR_Events_gq;
        TH1F * NumISR_Events_g;
        TH1F * NumISR_Events_gg_g;
        TH1F * NumISR_Events_gg_q;
        TH1F * NumISR_Events_qq_g;
        TH1F * NumISR_Events_qq_q;
        TH1F * NumISR_Events_gq_g;
        TH1F * NumISR_Events_gq_q;

};

// ------------------------------
// constructors and destructor
// ------------------------------
ISR_miniAOD::ISR_miniAOD(const edm::ParameterSet& iConfig) :
    genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("genParticles")))

{
    edm::Service<TFileService> fs;

    // ---------------------------------------
    // ISR definition - mom-dau: gg, qg, gq 
    // ---------------------------------------
    // before debug - first version - Semra
    miniAOD_ISR_Eta      = fs->make<TH1F>( "miniAOD_ISR_Eta",    "ISR Eta", 100, -6, 6  );
    miniAOD_ISR_Eta_gg   = fs->make<TH1F>( "miniAOD_ISR_Eta_gg", "ISR Eta", 100, -6, 6  );
    miniAOD_ISR_Eta_qg   = fs->make<TH1F>( "miniAOD_ISR_Eta_qg", "ISR Eta", 100, -6, 6  );
    miniAOD_ISR_Eta_gq   = fs->make<TH1F>( "miniAOD_ISR_Eta_gq", "ISR Eta", 100, -6, 6  );
    //miniAOD_ISR_Pt       = fs->make<TH1F>( "miniAOD_ISR_Pt",    "ISR Pt", 100,  0, 1000 );
    //miniAOD_ISR_Pt_gg    = fs->make<TH1F>( "miniAOD_ISR_Pt_gg", "ISR Pt", 100,  0, 1000 );
    //miniAOD_ISR_Pt_qg    = fs->make<TH1F>( "miniAOD_ISR_Pt_qg", "ISR Pt", 100,  0, 1000 );
    //miniAOD_ISR_Pt_gq    = fs->make<TH1F>( "miniAOD_ISR_Pt_gq", "ISR Pt", 100,  0, 1000 );
    
    // after debug - last version - Josh
    miniAOD_ISR_Eta_g     = fs->make<TH1F>( "miniAOD_ISR_Eta_g",     "ISR Eta", 100, -6, 6  );
    miniAOD_ISR_Eta_gg_g  = fs->make<TH1F>( "miniAOD_ISR_Eta_gg_g",  "ISR Eta", 100, -6, 6  );
    miniAOD_ISR_Eta_gg_q  = fs->make<TH1F>( "miniAOD_ISR_Eta_gg_q",  "ISR Eta", 100, -6, 6  );
    miniAOD_ISR_Eta_qq_g  = fs->make<TH1F>( "miniAOD_ISR_Eta_qq_g",  "ISR Eta", 100, -6, 6  );
    miniAOD_ISR_Eta_qq_q  = fs->make<TH1F>( "miniAOD_ISR_Eta_qq_q",  "ISR Eta", 100, -6, 6  );
    miniAOD_ISR_Eta_gq_gq = fs->make<TH1F>( "miniAOD_ISR_Eta_gq_gq", "ISR Eta", 100, -6, 6  );
    //miniAOD_ISR_Pt_g      = fs->make<TH1F>( "miniAOD_ISR_Pt_g",      "ISR Pt", 100,  0, 1000 );
    //miniAOD_ISR_Pt_gg_g   = fs->make<TH1F>( "miniAOD_ISR_Pt_gg_g",   "ISR Pt", 100,  0, 1000 );
    //miniAOD_ISR_Pt_gg_q   = fs->make<TH1F>( "miniAOD_ISR_Pt_gg_q",   "ISR Pt", 100,  0, 1000 );
    //miniAOD_ISR_Pt_qq_g   = fs->make<TH1F>( "miniAOD_ISR_Pt_qq_g",   "ISR Pt", 100,  0, 1000 );
    //miniAOD_ISR_Pt_qq_q   = fs->make<TH1F>( "miniAOD_ISR_Pt_qq_q",   "ISR Pt", 100,  0, 1000 );
    //miniAOD_ISR_Pt_gq_gq  = fs->make<TH1F>( "miniAOD_ISR_Pt_gq_gq",  "ISR Pt", 100,  0, 1000 );   
   
    // event counters
    totalEvents         = fs->make<TH1F>( "totalEvents",        "Events",        5, 0, 5 );
    NumISR_Events       = fs->make<TH1F>( "NumISR_Events",      "NumISR_Events", 5, 0, 5 );
    NumISR_Events_gg    = fs->make<TH1F>( "NumISR_Events_gg",   "NumISR_Events", 5, 0, 5 ); 
    NumISR_Events_qg    = fs->make<TH1F>( "NumISR_Events_qg",   "NumISR_Events", 5, 0, 5 );
    NumISR_Events_gq    = fs->make<TH1F>( "NumISR_Events_gq",   "NumISR_Events", 5, 0, 5 );
    NumISR_Events_g     = fs->make<TH1F>( "NumISR_Events_g",    "NumISR_Events", 5, 0, 5 );
    NumISR_Events_gg_g  = fs->make<TH1F>( "NumISR_Events_gg_g", "NumISR_Events", 5, 0, 5 );
    NumISR_Events_gg_q  = fs->make<TH1F>( "NumISR_Events_gg_q", "NumISR_Events", 5, 0, 5 );
    NumISR_Events_qq_g  = fs->make<TH1F>( "NumISR_Events_qq_g", "NumISR_Events", 5, 0, 5 );
    NumISR_Events_qq_q  = fs->make<TH1F>( "NumISR_Events_qq_q", "NumISR_Events", 5, 0, 5 );
    NumISR_Events_gq_g = fs->make<TH1F>( "NumISR_Events_gq_g","NumISR_Events", 5, 0, 5 );
    NumISR_Events_gq_q = fs->make<TH1F>( "NumISR_Events_gq_q","NumISR_Events", 5, 0, 5 );

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

    int totEvents = 0;
    int nISR   = 0, nISR_gg   = 0, nISR_qg   = 0, nISR_gq   = 0;
    int nISR_g = 0, nISR_gg_g = 0, nISR_gg_q = 0, nISR_qq_g = 0, nISR_qq_q = 0, nISR_gq_g = 0, nISR_gq_q = 0;

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
        int momStatus = mom->status();
        double momPt  = mom->pt();
        double dauPt  = p.p4().pt();
        double dauEta = p.p4().eta();

        // ---------------------------------------
        // before debug - first version - Semra
        // -- ISR definition - mom-dau: gg, qg, gq
        // ---------------------------------------
        if ( ( abs(momPdgId) == 21 && abs(dauPdgId) == 21 && status == 23 ) ||
           (   abs(momPdgId) <= 6  && abs(dauPdgId) == 21 && status == 23 ) ||  
           (   abs(momPdgId) == 21 && abs(dauPdgId) <= 6  && status == 23 ) )
        {
            miniAOD_ISR_Eta->Fill(dauEta); 
            //miniAOD_ISR_Pt->Fill(dauPt);
            nISR++;
        }

        // ----------------
        // categorized ones
        // ----------------
        if ( abs(momPdgId) == 21 && abs(dauPdgId) == 21 && status == 23 )
        {
            miniAOD_ISR_Eta_gg->Fill(dauEta);
            //miniAOD_ISR_Pt_gg->Fill(dauPt);
            nISR_gg++;
        }

        if ( abs(momPdgId) <= 6 && abs(dauPdgId) == 21 && status == 23 )
        {
            miniAOD_ISR_Eta_qg->Fill(dauEta);
            //miniAOD_ISR_Pt_qg->Fill(dauPt);
            nISR_qg++;
        }

        if ( abs(momPdgId) == 21 && abs(dauPdgId) <= 6  && status == 23 )
        {
            miniAOD_ISR_Eta_gq->Fill(dauEta);
            //miniAOD_ISR_Pt_gq->Fill(dauPt);
            nISR_gq++;
        }

        // ---------------------------------
        // after debug - last version - Josh
        // -- mother status code 21 
        // -- number of mothers 2 
        // -- not double counting
        //    -- NO asymmetry 
        // ---------------------------------
        if (momNum == 2)
        {  
            if ( ( abs(p.mother(0)->pdgId()) <= 6 || p.mother(0)->pdgId()  == 21 ) &&
                 ( abs(p.mother(1)->pdgId()) <= 6 || p.mother(1)->pdgId()  == 21 ) &&
                 (    p.mother(0)->status() == 21 && p.mother(1)->status() == 21 ) &&
                 ( abs(dauPdgId) <= 6             || abs(dauPdgId) == 21 ) && status == 23 )
            {
                miniAOD_ISR_Eta_g->Fill(dauEta);
                //miniAOD_ISR_Pt_g->Fill(dauPt);
                nISR_g++;
            } 
       
            // ----------------
            // categorized ones
            // ----------------
            if ( ((p.mother(0)->pdgId()  == 21) && (p.mother(1)->pdgId() == 21)) && 
                   p.mother(0)->status() == 21  && p.mother(1)->status() == 21   &&
                   abs(dauPdgId)         == 21  && status                == 23)
            {
                miniAOD_ISR_Eta_gg_g->Fill(dauEta);
                //miniAOD_ISR_Pt_gg_g->Fill(dauPt);
                nISR_gg_g++;
            }

            if ( ((p.mother(0)->pdgId()  == 21) && (p.mother(1)->pdgId() == 21)) &&  
                   p.mother(0)->status() == 21  && p.mother(1)->status() == 21   &&
                   abs(dauPdgId)         <= 6   && status                == 23 )
            {
                miniAOD_ISR_Eta_gg_q->Fill(dauEta);
                //miniAOD_ISR_Pt_gg_q->Fill(dauPt);
                nISR_gg_q++;
            }

            if ( (abs(p.mother(0)->pdgId()) <= 6  && abs(p.mother(1)->pdgId()) <= 6) &&  
                  p.mother(0)->status()     == 21 && p.mother(1)->status()     == 21 &&
                  abs(dauPdgId)             == 21 && status                    == 23 )
            {
                miniAOD_ISR_Eta_qq_g->Fill(dauEta);
                //miniAOD_ISR_Pt_qq_g->Fill(dauPt);
                nISR_qq_g++;
            }

            if ( (abs(p.mother(0)->pdgId()) <= 6  && abs(p.mother(1)->pdgId()) <= 6) && 
                  p.mother(0)->status()     == 21 && p.mother(1)->status()     == 21 && 
                  abs(dauPdgId)             <= 6  && status                    == 23 )
            {
                miniAOD_ISR_Eta_qq_q->Fill(dauEta);
                //miniAOD_ISR_Pt_qq_q->Fill(dauPt);
                nISR_qq_q++;
            }

            if ( ( ((p.mother(0)->pdgId()  <= 6)  && (p.mother(1)->pdgId() == 21))   || 
                   ((p.mother(0)->pdgId()  == 21) && (p.mother(1)->pdgId() <= 6 )) ) &&
                   p.mother(0)->status()   == 21  && p.mother(1)->status() == 21     &&
                 ( abs(dauPdgId)         == 21 )   && status == 23 )
            {
                //miniAOD_ISR_Eta_gq_g->Fill(dauEta);
                nISR_gq_g++;    
            }

            if ( ( ((p.mother(0)->pdgId()  <= 6)  && (p.mother(1)->pdgId() == 21))   || 
                   ((p.mother(0)->pdgId()  == 21) && (p.mother(1)->pdgId() <= 6 )) ) &&
                   p.mother(0)->status()   == 21  && p.mother(1)->status() == 21     &&
                 ( abs(dauPdgId)           <= 6 )   && status == 23 )
            {
                //miniAOD_ISR_Eta_gq_q->Fill(dauEta);
                nISR_gq_q++;    
            }

        }
    } // gen particles loop
   
    // event counters
    totEvents++;
    totalEvents->Fill(totEvents);

    NumISR_Events->Fill(nISR);
    NumISR_Events_gg->Fill(nISR_gg);
    NumISR_Events_qg->Fill(nISR_qg);
    NumISR_Events_gq->Fill(nISR_gq);

    NumISR_Events_g->Fill(nISR_g);
    NumISR_Events_gg_g->Fill(nISR_gg_g);
    NumISR_Events_gg_q->Fill(nISR_gg_q);
    NumISR_Events_qq_g->Fill(nISR_qq_g);
    NumISR_Events_qq_q->Fill(nISR_qq_q);
    NumISR_Events_gq_g->Fill(nISR_gq_g);
    NumISR_Events_gq_q->Fill(nISR_gq_q);


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




