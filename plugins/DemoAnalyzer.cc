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

#include "TFile.h"
#include "TTree.h"

// --------------------
// class declaration
// --------------------
class DemoAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  
{
    public:
        explicit DemoAnalyzer(const edm::ParameterSet&);
        ~DemoAnalyzer();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        // -----------
        // Member data
        // ----------- 
        edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;        
        
        // mother:proton - daughter:quark - status:71,73
        TH1F * NumMom_pq_71;
        TH1F * NumDau_pq_71;
        TH1F * NumGrandDau_pq_71;
        TH1F * MomPt_pq_71;
        TH1F * DauPt_pq_71;        

        TH1F * NumMom_pq_73; 
        TH1F * NumDau_pq_73;
        TH1F * NumGrandDau_pq_73;
        TH1F * MomPt_pq_73; 
        TH1F * DauPt_pq_73;
        
        // mother:proton - daughter:gluon - status:71,73
        TH1F * NumMom_pg_71;
        TH1F * NumDau_pg_71;
        TH1F * NumGrandDau_pg_71;
        TH1F * MomPt_pg_71;
        TH1F * DauPt_pg_71;

        TH1F * NumMom_pg_73;
        TH1F * NumDau_pg_73; 
        TH1F * NumGrandDau_pg_73;
        TH1F * MomPt_pg_73; 
        TH1F * DauPt_pg_73;
 
        // mother:quark - daughter:quark - status:23,44,71,73
        TH1F * NumMom_qq_23;
        TH1F * NumDau_qq_23; 
        TH1F * NumGrandDau_qq_23;
        TH1F * MomPt_qq_23; 
        TH1F * DauPt_qq_23;

        TH1F * NumMom_qq_44;
        TH1F * NumDau_qq_44;
        TH1F * NumGrandDau_qq_44; 
        TH1F * MomPt_qq_44; 
        TH1F * DauPt_qq_44;

        TH1F * NumMom_qq_71;
        TH1F * NumDau_qq_71;
        TH1F * NumGrandDau_qq_71;
        TH1F * MomPt_qq_71;
        TH1F * DauPt_qq_71;

        TH1F * NumMom_qq_73;
        TH1F * NumDau_qq_73; 
        TH1F * NumGrandDau_qq_73;
        TH1F * MomPt_qq_73; 
        TH1F * DauPt_qq_73;

        // mother:quark - daughter:gluon - status:23,71,73
        TH1F * NumMom_qg_23;
        TH1F * NumDau_qg_23; 
        TH1F * NumGrandDau_qg_23;
        TH1F * MomPt_qg_23; 
        TH1F * DauPt_qg_23; 

        TH1F * NumMom_qg_71;
        TH1F * NumDau_qg_71;
        TH1F * NumGrandDau_qg_71;
        TH1F * MomPt_qg_71;
        TH1F * DauPt_qg_71;

        TH1F * NumMom_qg_73;
        TH1F * NumDau_qg_73; 
        TH1F * NumGrandDau_qg_73;
        TH1F * MomPt_qg_73; 
        TH1F * DauPt_qg_73;       

        // mother:gluon - daughter:gluon - status:23,44,71,73
        TH1F * NumMom_gg_23;
        TH1F * NumDau_gg_23; 
        TH1F * NumGrandDau_gg_23;
        TH1F * MomPt_gg_23; 
        TH1F * DauPt_gg_23; 

        TH1F * NumMom_gg_44;
        TH1F * NumDau_gg_44;
        TH1F * NumGrandDau_gg_44; 
        TH1F * MomPt_gg_44; 
        TH1F * DauPt_gg_44;

        TH1F * NumMom_gg_71;
        TH1F * NumDau_gg_71;
        TH1F * NumGrandDau_gg_71;
        TH1F * MomPt_gg_71;
        TH1F * DauPt_gg_71;

        TH1F * NumMom_gg_73;
        TH1F * NumDau_gg_73; 
        TH1F * NumGrandDau_gg_73;
        TH1F * MomPt_gg_73; 
        TH1F * DauPt_gg_73;

        // mother:gluon - daughter:quark - status:23,71,73
        TH1F * NumMom_gq_23;
        TH1F * NumDau_gq_23; 
        TH1F * NumGrandDau_gq_23;
        TH1F * MomPt_gq_23; 
        TH1F * DauPt_gq_23; 

        TH1F * NumMom_gq_71;
        TH1F * NumDau_gq_71;
        TH1F * NumGrandDau_gq_71;
        TH1F * MomPt_gq_71;
        TH1F * DauPt_gq_71;

        TH1F * NumMom_gq_73;
        TH1F * NumDau_gq_73; 
        TH1F * NumGrandDau_gq_73;
        TH1F * MomPt_gq_73; 
        TH1F * DauPt_gq_73;
};

// ------------------------------
// constructors and destructor
// ------------------------------
DemoAnalyzer::DemoAnalyzer(const edm::ParameterSet& iConfig) :
    genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("genParticles")))

{
    edm::Service<TFileService> fs;

    // mother:proton - daughter:quark - status:71,73
    NumMom_pq_71      = fs->make<TH1F>( "NumMom_pq_71", "Number of Mother", 50,  0, 50 );
    NumDau_pq_71      = fs->make<TH1F>( "NumDau_pq_71", "Number of Daughter", 50,  0, 50 );
    NumGrandDau_pq_71 = fs->make<TH1F>( "NumGrandDau_pq_71", "Number of GrandDaughter", 50,  0, 50 );
    MomPt_pq_71       = fs->make<TH1F>( "MomPt_pq_71", "Mom p_{T}", 100,  0, 1000 );
    DauPt_pq_71       = fs->make<TH1F>( "DauPt_pq_71", "Daughter p_{T}", 100,  0, 1000 );
 
    NumMom_pq_73      = fs->make<TH1F>( "NumMom_pq_73", "Number of Mother", 50,  0, 50 );
    NumDau_pq_73      = fs->make<TH1F>( "NumDau_pq_73", "Number of Daughter", 50,  0, 50 ); 
    NumGrandDau_pq_73 = fs->make<TH1F>( "NumGrandDau_pq_73", "Number of GrandDaughter", 50,  0, 50 );
    MomPt_pq_73       = fs->make<TH1F>( "MomPt_pq_73", "Mom p_{T}", 100,  0, 1000 );
    DauPt_pq_73       = fs->make<TH1F>( "DauPt_pq_73", "Daughter p_{T}", 100,  0, 1000 ); 
 
    // mother:proton - daughter:gluon - status:71,73
    NumMom_pg_71      = fs->make<TH1F>( "NumMom_pg_71", "Number of Mother", 50,  0, 50 );
    NumDau_pg_71      = fs->make<TH1F>( "NumDau_pg_71", "Number of Daughter", 50,  0, 50 );
    NumGrandDau_pg_71 = fs->make<TH1F>( "NumGrandDau_pg_71", "Number of GrandDaughter", 50,  0, 50 );
    MomPt_pg_71       = fs->make<TH1F>( "MomPt_pg_71", "Mom p_{T}", 100,  0, 1000 );
    DauPt_pg_71       = fs->make<TH1F>( "DauPt_pg_71", "Daughter p_{T}", 100,  0, 1000 );
    
    NumMom_pg_73      = fs->make<TH1F>( "NumMom_pg_73", "Number of Mother", 50,  0, 50 );
    NumDau_pg_73      = fs->make<TH1F>( "NumDau_pg_73", "Number of Daughter", 50,  0, 50 );
    NumGrandDau_pg_73 = fs->make<TH1F>( "NumGrandDau_pg_73", "Number of GrandDaughter", 50,  0, 50 );
    MomPt_pg_73       = fs->make<TH1F>( "MomPt_pg_73", "Mom p_{T}", 100,  0, 1000 );
    DauPt_pg_73       = fs->make<TH1F>( "DauPt_pg_73", "Daughter p_{T}", 100,  0, 1000 );

    // mother:quark - daughter:quark - status:23,44,71,73
    NumMom_qq_23      = fs->make<TH1F>( "NumMom_qq_23", "Number of Mother", 50,  0, 50 );
    NumDau_qq_23      = fs->make<TH1F>( "NumDau_qq_23", "Number of Daughter", 50,  0, 50 );
    NumGrandDau_qq_23 = fs->make<TH1F>( "NumGrandDau_qq_23", "Number of GrandDaughter", 50,  0, 50 );
    MomPt_qq_23       = fs->make<TH1F>( "MomPt_qq_23", "Mom p_{T}", 100,  0, 1000 ); 
    DauPt_qq_23       = fs->make<TH1F>( "DauPt_qq_23", "Daughter p_{T}", 100,  0, 1000 );
    
    NumMom_qq_44      = fs->make<TH1F>( "NumMom_qq_44", "Number of Mother", 50,  0, 50 );
    NumDau_qq_44      = fs->make<TH1F>( "NumDau_qq_44", "Number of Daughter", 50,  0, 50 ); 
    NumGrandDau_qq_44 = fs->make<TH1F>( "NumGrandDau_qq_44", "Number of GrandDaughter", 50,  0, 50 );
    MomPt_qq_44       = fs->make<TH1F>( "MomPt_qq_44", "Mom p_{T}", 100,  0, 1000 ); 
    DauPt_qq_44       = fs->make<TH1F>( "DauPt_qq_44", "Daughter p_{T}", 100,  0, 1000 );
    
    NumMom_qq_71      = fs->make<TH1F>( "NumMom_qq_71", "Number of Mother", 50,  0, 50 );
    NumDau_qq_71      = fs->make<TH1F>( "NumDau_qq_71", "Number of Daughter", 50,  0, 50 );
    NumGrandDau_qq_71 = fs->make<TH1F>( "NumGrandDau_qq_71", "Number of GrandDaughter", 50,  0, 50 );
    MomPt_qq_71       = fs->make<TH1F>( "MomPt_qq_71", "Mom p_{T}", 100,  0, 1000 );
    DauPt_qq_71       = fs->make<TH1F>( "DauPt_qq_71", "Daughter p_{T}", 100,  0, 1000 );
    
    NumMom_qq_73      = fs->make<TH1F>( "NumMom_qq_73", "Number of Mother", 50,  0, 50 );
    NumDau_qq_73      = fs->make<TH1F>( "NumDau_qq_73", "Number of Daughter", 50,  0, 50 ); 
    NumGrandDau_qq_73 = fs->make<TH1F>( "NumGrandDau_qq_73", "Number of GrandDaughter", 50,  0, 50 );
    MomPt_qq_73       = fs->make<TH1F>( "MomPt_qq_73", "Mom p_{T}", 100,  0, 1000 ); 
    DauPt_qq_73       = fs->make<TH1F>( "DauPt_qq_73", "Daughter p_{T}", 100,  0, 1000 );
 
    // mother:quark - daughter:gluon - status:23,71,73   
    NumMom_qg_23      = fs->make<TH1F>( "NumMom_qg_23", "Number of Mother", 50,  0, 50 );
    NumDau_qg_23      = fs->make<TH1F>( "NumDau_qg_23", "Number of Daughter", 50,  0, 50 );
    NumGrandDau_qg_23 = fs->make<TH1F>( "NumGrandDau_qg_23", "Number of GrandDaughter", 50,  0, 50 );
    MomPt_qg_23       = fs->make<TH1F>( "MomPt_qg_23", "Mom p_{T}", 100,  0, 1000 ); 
    DauPt_qg_23       = fs->make<TH1F>( "DauPt_qg_23", "Daughter p_{T}", 100,  0, 1000 ); 
    
    NumMom_qg_71      = fs->make<TH1F>( "NumMom_qg_71", "Number of Mother", 50,  0, 50 );
    NumDau_qg_71      = fs->make<TH1F>( "NumDau_qg_71", "Number of Daughter", 50,  0, 50 );
    NumGrandDau_qg_71 = fs->make<TH1F>( "NumGrandDau_qg_71", "Number of GrandDaughter", 50,  0, 50 );
    MomPt_qg_71       = fs->make<TH1F>( "MomPt_qg_71", "Mom p_{T}", 100,  0, 1000 );
    DauPt_qg_71       = fs->make<TH1F>( "DauPt_qg_71", "Daughter p_{T}", 100,  0, 1000 );
    
    NumMom_qg_73      = fs->make<TH1F>( "NumMom_qg_73", "Number of Mother", 50,  0, 50 );
    NumDau_qg_73      = fs->make<TH1F>( "NumDau_qg_73", "Number of Daughter", 50,  0, 50 ); 
    NumGrandDau_qg_73 = fs->make<TH1F>( "NumGrandDau_qg_73", "Number of GrandDaughter", 50,  0, 50 );
    MomPt_qg_73       = fs->make<TH1F>( "MomPt_qg_73", "Mom p_{T}", 100,  0, 1000 ); 
    DauPt_qg_73       = fs->make<TH1F>( "DauPt_qg_73", "Daughter p_{T}", 100,  0, 1000 ); 

    // mother:gluon - daughter:gluon - status:23,44,71,73   
    NumMom_gg_23      = fs->make<TH1F>( "NumMom_gg_23", "Number of Mother", 50,  0, 50 );
    NumDau_gg_23      = fs->make<TH1F>( "NumDau_gg_23", "Number of Daughter", 50,  0, 50 ); 
    NumGrandDau_gg_23 = fs->make<TH1F>( "NumGrandDau_gg_23", "Number of GrandDaughter", 50,  0, 50 );
    MomPt_gg_23       = fs->make<TH1F>( "MomPt_gg_23", "Mom p_{T}", 100,  0, 1000 ); 
    DauPt_gg_23       = fs->make<TH1F>( "DauPt_gg_23", "Daughter p_{T}", 100,  0, 1000 );
    
    NumMom_gg_44      = fs->make<TH1F>( "NumMom_gg_44", "Number of Mother", 50,  0, 50 );
    NumDau_gg_44      = fs->make<TH1F>( "NumDau_gg_44", "Number of Daughter", 50,  0, 50 );
    NumGrandDau_gg_44 = fs->make<TH1F>( "NumGrandDau_gg_44", "Number of GrandDaughter", 50,  0, 50 );
    MomPt_gg_44       = fs->make<TH1F>( "MomPt_gg_44", "Mom p_{T}", 100,  0, 1000 ); 
    DauPt_gg_44       = fs->make<TH1F>( "DauPt_gg_44", "Daughter p_{T}", 100,  0, 1000 ); 
    
    NumMom_gg_71      = fs->make<TH1F>( "NumMom_gg_71", "Number of Mother", 50,  0, 50 );
    NumDau_gg_71      = fs->make<TH1F>( "NumDau_gg_71", "Number of Daughter", 50,  0, 50 );
    NumGrandDau_gg_71 = fs->make<TH1F>( "NumGrandDau_gg_71", "Number of GrandDaughter", 50,  0, 50 );
    MomPt_gg_71       = fs->make<TH1F>( "MomPt_gg_71", "Mom p_{T}", 100,  0, 1000 );
    DauPt_gg_71       = fs->make<TH1F>( "DauPt_gg_71", "Daughter p_{T}", 100,  0, 1000 );
    
    NumMom_gg_73      = fs->make<TH1F>( "NumMom_gg_73", "Number of Mother", 50,  0, 50 );
    NumDau_gg_73      = fs->make<TH1F>( "NumDau_gg_73", "Number of Daughter", 50,  0, 50 ); 
    NumGrandDau_gg_73 = fs->make<TH1F>( "NumGrandDau_gg_73", "Number of GrandDaughter", 50,  0, 50 );
    MomPt_gg_73       = fs->make<TH1F>( "MomPt_gg_73", "Mom p_{T}", 100,  0, 1000 ); 
    DauPt_gg_73       = fs->make<TH1F>( "DauPt_gg_73", "Daughter p_{T}", 100,  0, 1000 );

    // mother:gluon - daughter:quark - status:23,71,73
    NumMom_gq_23      = fs->make<TH1F>( "NumMom_gq_23", "Number of Mother", 50,  0, 50 );
    NumDau_gq_23      = fs->make<TH1F>( "NumDau_gq_23", "Number of Daughter", 50,  0, 50 );
    NumGrandDau_gq_23 = fs->make<TH1F>( "NumGrandDau_gq_23", "Number of GrandDaughter", 50,  0, 50 );
    MomPt_gq_23       = fs->make<TH1F>( "MomPt_gq_23", "Mom p_{T}", 100,  0, 1000 ); 
    DauPt_gq_23       = fs->make<TH1F>( "DauPt_gq_23", "Daughter p_{T}", 100,  0, 1000 ); 
 
    NumMom_gq_71      = fs->make<TH1F>( "NumMom_gq_71", "Number of Mother", 50,  0, 50 );
    NumDau_gq_71      = fs->make<TH1F>( "NumDau_gq_71", "Number of Daughter", 50,  0, 50 );
    NumGrandDau_gq_71 = fs->make<TH1F>( "NumGrandDau_gq_71", "Number of GrandDaughter", 50,  0, 50 );
    MomPt_gq_71       = fs->make<TH1F>( "MomPt_gq_71", "Mom p_{T}", 100,  0, 1000 );
    DauPt_gq_71       = fs->make<TH1F>( "DauPt_gq_71", "Daughter p_{T}", 100,  0, 1000 );
    
    NumMom_gq_73      = fs->make<TH1F>( "NumMom_gq_73", "Number of Mother", 50,  0, 50 );
    NumDau_gq_73      = fs->make<TH1F>( "NumDau_gq_73", "Number of Daughter", 50,  0, 50 );
    NumGrandDau_gq_73 = fs->make<TH1F>( "NumGrandDau_gq_73", "Number of GrandDaughter", 50,  0, 50 );
    MomPt_gq_73       = fs->make<TH1F>( "MomPt_gq_73", "Mom p_{T}", 100,  0, 1000 ); 
    DauPt_gq_73       = fs->make<TH1F>( "DauPt_gq_73", "Daughter p_{T}", 100,  0, 1000 );

}

DemoAnalyzer::~DemoAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

// -------------------------------
// Member functions
// method called for each event 
// -------------------------------
void DemoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
        int momPdgId = mom->pdgId();
        int dauNum   = mom->numberOfDaughters();
        double momPt = mom->pt();
        double dauPt = p.p4().pt();

        //for (int d = 0; d < dauNum; d++)
        //{          
            //const reco::Candidate *dau = p.daughter(d);
            //int dauPdgId  = dau->pdgId();
            //int dauNumber = dau->numberOfDaughters();

            // ----------------
            // print out things
            // ----------------
            //if ( ( abs(dauPdgId) <= 6      || abs(dauPdgId) == 21 )
            //       && ( abs(momPdgId) <= 6 || abs(momPdgId) == 21            || abs(momPdgId) == 2212 )
            //       && ( status == 23       || (status <= 49 && status >= 41) || (status <= 79 && status >= 71) ) )
            //{
            //    std::cout << "momPdgId:  " << momPdgId  << " || " 
            //              << "dauPdgId:  " << dauPdgId  << " || "  
            //              << "status:    " << status    << std::endl;
            //}       
 
            // --------------------------------------------- 
            // mother:proton - daughter:quark - status:71,73
            // --------------------------------------------- 
            if ( abs(momPdgId) == 2212 && abs(dauPdgId) <= 6 && status == 71 )
            {
                NumMom_pq_71->Fill(momNum);
                NumDau_pq_71->Fill(dauNum);
                NumGrandDau_pq_71->Fill(grandDauNum); 
                MomPt_pq_71->Fill(momPt);
                DauPt_pq_71->Fill(dauPt);
            }            

            if ( abs(momPdgId) == 2212 && abs(dauPdgId) <= 6 && status == 73 )
            {
                NumMom_pq_73->Fill(momNum);
                NumDau_pq_73->Fill(dauNum);
                NumGrandDau_pq_73->Fill(grandDauNum);
                MomPt_pq_73->Fill(momPt);
                DauPt_pq_73->Fill(dauPt);
            }            
            
            // --------------------------------------------- 
            // mother:proton - daughter:gluon - status:71,73
            // --------------------------------------------- 
            if ( abs(momPdgId) == 2212 && abs(dauPdgId) == 21 && status == 71 )
            {
                NumMom_pg_71->Fill(momNum);
                NumDau_pg_71->Fill(dauNum);
                NumGrandDau_pg_71->Fill(grandDauNum);
                MomPt_pg_71->Fill(momPt);
                DauPt_pg_71->Fill(dauPt);
            }            

            if ( abs(momPdgId) == 2212 && abs(dauPdgId) == 21 && status == 73 )
            {
                NumMom_pg_73->Fill(momNum);
                NumDau_pg_73->Fill(dauNum);
                NumGrandDau_pg_73->Fill(grandDauNum);
                MomPt_pg_73->Fill(momPt);
                DauPt_pg_73->Fill(dauPt);
            }
    
            // --------------------------------------------------
            // mother:quark - daughter:quark - status:23,44,71,73
            // -------------------------------------------------- 
            if ( abs(momPdgId) <= 6 && abs(dauPdgId) <= 6 && status == 23 )
            {
                NumMom_qq_23->Fill(momNum);
                NumDau_qq_23->Fill(dauNum);
                NumGrandDau_qq_23->Fill(grandDauNum);
                MomPt_qq_23->Fill(momPt);
                DauPt_qq_23->Fill(dauPt);
            }
            
            if ( abs(momPdgId) <= 6 && abs(dauPdgId) <= 6 && status == 44 )
            {
                NumMom_qq_44->Fill(momNum);
                NumDau_qq_44->Fill(dauNum);
                NumGrandDau_qq_44->Fill(grandDauNum);
                MomPt_qq_44->Fill(momPt);
                DauPt_qq_44->Fill(dauPt);
            }

            if ( abs(momPdgId) <= 6 && abs(dauPdgId) <= 6 && status == 71 )
            {
                NumMom_qq_71->Fill(momNum);
                NumDau_qq_71->Fill(dauNum);
                NumGrandDau_qq_71->Fill(grandDauNum);
                MomPt_qq_71->Fill(momPt);
                DauPt_qq_71->Fill(dauPt);
            }

            if ( abs(momPdgId) <= 6 && abs(dauPdgId) <= 6 && status == 73 )
            {
                NumMom_qq_73->Fill(momNum);   
                NumDau_qq_73->Fill(dauNum);
                NumGrandDau_qq_73->Fill(grandDauNum);
                MomPt_qq_73->Fill(momPt);
                DauPt_qq_73->Fill(dauPt);    
            }

            // ----------------------------------------------- 
            // mother:quark - daughter:gluon - status:23,71,73          
            // ----------------------------------------------- 
            if ( abs(momPdgId) <= 6 && abs(dauPdgId) == 21 && status == 23 )
            {
                NumMom_qg_23->Fill(momNum);
                NumDau_qg_23->Fill(dauNum);
                NumGrandDau_qg_23->Fill(grandDauNum);
                MomPt_qg_23->Fill(momPt);
                DauPt_qg_23->Fill(dauPt);
            }
            
            if ( abs(momPdgId) <= 6 && abs(dauPdgId) == 21 && status == 71 )
            {
                NumMom_qg_71->Fill(momNum);
                NumDau_qg_71->Fill(dauNum);
                NumGrandDau_qg_71->Fill(grandDauNum);
                MomPt_qg_71->Fill(momPt);
                DauPt_qg_71->Fill(dauPt);
            }

            if ( abs(momPdgId) <= 6 && abs(dauPdgId) == 21 && status == 73 )
            {
                NumMom_qg_73->Fill(momNum);   
                NumDau_qg_73->Fill(dauNum);
                NumGrandDau_qg_73->Fill(grandDauNum);
                MomPt_qg_73->Fill(momPt);
                DauPt_qg_73->Fill(dauPt);       
            }

            // -------------------------------------------------- 
            // mother:gluon - daughter:gluon - status:23,44,71,73
            // -------------------------------------------------- 
            if ( abs(momPdgId) == 21 && abs(dauPdgId) == 21 && status == 23 )
            {
                NumMom_gg_23->Fill(momNum);
                NumDau_gg_23->Fill(dauNum);
                NumGrandDau_gg_23->Fill(grandDauNum);
                MomPt_gg_23->Fill(momPt);
                DauPt_gg_23->Fill(dauPt);
            }

            if ( abs(momPdgId) == 21 && abs(dauPdgId) == 21 && status == 44 )
            {
                NumMom_gg_44->Fill(momNum);
                NumDau_gg_44->Fill(dauNum);
                NumGrandDau_gg_44->Fill(grandDauNum);
                MomPt_gg_44->Fill(momPt);
                DauPt_gg_44->Fill(dauPt);                
            }

            if ( abs(momPdgId) == 21 && abs(dauPdgId) == 21 && status == 71 )
            {
                NumMom_gg_71->Fill(momNum);
                NumDau_gg_71->Fill(dauNum);
                NumGrandDau_gg_71->Fill(grandDauNum);
                MomPt_gg_71->Fill(momPt);
                DauPt_gg_71->Fill(dauPt);
            }

            if ( abs(momPdgId) == 21 && abs(dauPdgId) == 21 && status == 73 )
            {
                NumMom_gg_73->Fill(momNum);
                NumDau_gg_73->Fill(dauNum);
                NumGrandDau_gg_73->Fill(grandDauNum);
                MomPt_gg_73->Fill(momPt);
                DauPt_gg_73->Fill(dauPt);           
            }   

            // ----------------------------------------------- 
            // mother:gluon - daughter:quark - status:23,71,73
            // ----------------------------------------------- 
            if ( abs(momPdgId) == 21 && abs(dauPdgId) <= 6 && status == 23 )
            {
                NumMom_gq_23->Fill(momNum);
                NumDau_gq_23->Fill(dauNum);
                NumGrandDau_gq_23->Fill(grandDauNum);
                MomPt_gq_23->Fill(momPt);
                DauPt_gq_23->Fill(dauPt);                    
            }

            if ( abs(momPdgId) == 21 && abs(dauPdgId) <= 6 && status == 71 )
            {
                NumMom_gq_71->Fill(momNum);
                NumDau_gq_71->Fill(dauNum);
                NumGrandDau_gq_71->Fill(grandDauNum);
                MomPt_gq_71->Fill(momPt);
                DauPt_gq_71->Fill(dauPt);
            }

            if ( abs(momPdgId) == 21 && abs(dauPdgId) <= 6 && status == 73 )
            {
                NumMom_gq_73->Fill(momNum);
                NumDau_gq_73->Fill(dauNum);
                NumGrandDau_gq_73->Fill(grandDauNum);
                MomPt_gq_73->Fill(momPt);
                DauPt_gq_73->Fill(dauPt);
            }
    
        //} // daughher loop

    } // gen particles loop
    
} // event loop

// --------------------------------------------------------------
// Method called once each job just before starting event loop  
// --------------------------------------------------------------
void DemoAnalyzer::beginJob()
{
}

// ---------------------------------------------------------------
// Method called once each job just after ending the event loop  
// ---------------------------------------------------------------
void DemoAnalyzer::endJob()
{
}

// -------------------------------------------------------------------------
// Method fills 'descriptions' with the allowed parameters for the module 
// -------------------------------------------------------------------------
void DemoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
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
DEFINE_FWK_MODULE(DemoAnalyzer);




