// -*- C++ -*-
//
// Package:    FastPUPPI/SeedJetProducer
// Class:      SeedJetProducer
// 
/**\class SeedJetProducer SeedJetProducer.cc FastPUPPI/SeedJetProducer/plugins/SeedJetProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  zhenbin wu
//         Created:  Fri, 27 Oct 2017 20:47:24 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

// The input from the P2L1 are in l1t:PFCandidate type
#include "DataFormats/Phase2L1ParticleFlow/interface/PFCandidate.h"

// The output is expected to be Reco::Candidate type for storing
#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "RecoJets/JetProducers/interface/JetSpecific.h"

#include "TLorentzVector.h"
//
// class declaration
//


class SeedJetProducer : public edm::stream::EDProducer<> {
   public:
      explicit SeedJetProducer(const edm::ParameterSet&);
      ~SeedJetProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      bool ProdPar(edm::Handle<l1t::PFCandidateCollection> ParHdl);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;


      // ----------member data ---------------------------
      edm::EDGetTokenT<l1t::PFCandidateCollection> ParTok_;
      edm::Handle<l1t::PFCandidateCollection> ParHdl;
      std::unique_ptr<reco::PFJetCollection> jetCollection;
      // Clusting
      double seedPt;
      unsigned seedID;
      double cone;

      bool SeedCone(std::vector<unsigned> seeds, const edm::EventSetup& iSetup);
      std::vector<unsigned> GetSeed(const int pt, const int &id);




};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
SeedJetProducer::SeedJetProducer(const edm::ParameterSet& iConfig):
  ParTok_( consumes<l1t::PFCandidateCollection> (iConfig.getParameter<edm::InputTag>("src")) ),
  seedPt(iConfig.getParameter<double>       ("seedPtCut")),
  seedID(iConfig.getParameter<unsigned>       ("seedID")),
  cone(iConfig.getParameter<double>       ("cone"))
{
   produces<reco::PFJetCollection>();
}


SeedJetProducer::~SeedJetProducer()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
SeedJetProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

  iEvent.getByToken(ParTok_, ParHdl); 

  SeedCone(GetSeed(seedPt, seedID), iSetup);

  iEvent.put(std::move(jetCollection));

}


// ===  FUNCTION  ============================================================
//         Name:  SeedJetProducer::GetSeed
//  Description:  
// ===========================================================================
std::vector<unsigned> SeedJetProducer::GetSeed(const int pt, const int &id)
{
  unsigned n= 0;
  std::vector<unsigned> idx;
  std::map< double, unsigned > ptmap;

  for (const l1t::PFCandidate &g : *ParHdl) {
    if (g.pt() >= pt)
    {
        ptmap[g.pt()] = n;
    }
    n++;
  }


  for(std::map<double, unsigned>::reverse_iterator i=ptmap.rbegin(); i!=ptmap.rend(); ++i)
  {
    idx.push_back(i->second);
  }
  return idx;
}       // -----  end of function SeedJetProducer::GetSeed  -----


// ===  FUNCTION  ============================================================
//         Name:  SeedJetProducer::SeedCone
//  Description:  
// ===========================================================================
bool SeedJetProducer::SeedCone(std::vector<unsigned> seeds, const edm::EventSetup& iSetup)
{

  std::vector<bool> taken(ParHdl->size(), false);
  std::vector<TLorentzVector> rejets;

  for(auto i : seeds)
  {
    if (taken[i]) continue;
    TLorentzVector seed(0, 0, 0, 0);
    seed.SetPtEtaPhiM(ParHdl->at(i).pt(), ParHdl->at(i).eta(), ParHdl->at(i).phi(), 0);
    taken[i] = true;

    std::vector<TLorentzVector> incone;
    for(unsigned int j=0; j < ParHdl->size(); ++j)
    {
      if (taken[j]) continue;
      const l1t::PFCandidate &p = ParHdl->at(j);
      if ( (pow((seed.Eta() - p.eta()), 2) + pow((seed.Phi() - p.phi()), 2) ) < cone * cone)
      {
        TLorentzVector par(0, 0, 0, 0);
        par.SetPtEtaPhiM(p.pt(), p.eta(), p.phi(), 0);
        incone.push_back(par);
        taken[j] = true;
      }
    }

    for(auto jj : incone)
    {
      seed += jj;
    }
    rejets.push_back(seed);
  }


  jetCollection = std::make_unique<reco::PFJetCollection>();
  std::vector<math::XYZTLorentzVector> p4_Jets;
  for(auto g : rejets)
  {
    reco::PFJet jet;
    reco::Particle::Point           vertex_(0,0,0);
    std::vector<reco::CandidatePtr> constituents;
    reco::writeSpecific(jet,
                  reco::Particle::LorentzVector( g.Px(), g.Py(), g.Pz(), g.E()),
                  vertex_,
                  constituents,
                  iSetup);
    jetCollection->emplace_back( jet );
  }

  return true;
}       // -----  end of function SeedJetProducer::SeedCone  -----

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
SeedJetProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
SeedJetProducer::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
SeedJetProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
SeedJetProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
SeedJetProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
SeedJetProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SeedJetProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SeedJetProducer);
