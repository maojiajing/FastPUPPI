// -*- C++ -*-
//
// Package:    FastPUPPI/GenJetProducer
// Class:      GenJetProducer
// 
/**\class GenJetProducer GenJetProducer.cc FastPUPPI/GenJetProducer/plugins/GenJetProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  zhenbin wu
//         Created:  Sat, 30 Sep 2017 17:17:57 GMT
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

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "TTree.h"


//
// class declaration
//

class GenJetProducer : public edm::stream::EDProducer<> {
   public:
      explicit GenJetProducer(const edm::ParameterSet&);
      ~GenJetProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      bool ProdGenPar(edm::Handle<reco::GenParticleCollection> GenParHdl);
      bool ProdGenJet(edm::Handle<reco::GenJetCollection> GenJetHdl);
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      TTree *tree_;
      edm::EDGetTokenT<reco::GenJetCollection>      GenJetTok_;
      edm::EDGetTokenT<reco::GenParticleCollection> GenParTok_;
      std::vector<float> *genjetPt;
      std::vector<float> *genjetEta;
      std::vector<float> *genjetPhi;
      std::vector<float> *genjetMass;

      std::vector<float> *genparPt;
      std::vector<float> *genparEta;
      std::vector<float> *genparPhi;
      std::vector<float> *genparMass;
      std::vector<int> *genparID;
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
GenJetProducer::GenJetProducer(const edm::ParameterSet& iConfig):
  GenJetTok_( consumes<reco::GenJetCollection> (iConfig.getParameter<edm::InputTag>("GenJetTag")) ),
  GenParTok_( consumes<reco::GenParticleCollection> (iConfig.getParameter<edm::InputTag>("GenParTag")) )
{

  genjetPt   = new std::vector<float>();
  genjetEta  = new std::vector<float>();
  genjetPhi  = new std::vector<float>();
  genjetMass = new std::vector<float>();


  genparPt   = new std::vector<float>();
  genparEta  = new std::vector<float>();
  genparPhi  = new std::vector<float>();
  genparMass  = new std::vector<float>();
  genparID  = new std::vector<int>();

   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed
  
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree","tree");
  tree_->Branch("GenJet_pt", "vector<float>", genjetPt);
  tree_->Branch("GenJet_eta", "vector<float>", genjetEta);
  tree_->Branch("GenJet_phi", "vector<float>", genjetPhi);
  tree_->Branch("GenJet_mass", "vector<float>", genjetMass);
  tree_->Branch("GenPar_pt", "vector<float>", genparPt);
  tree_->Branch("GenPar_eta", "vector<float>", genparEta);
  tree_->Branch("GenPar_phi", "vector<float>", genparPhi);
  tree_->Branch("GenPar_mass", "vector<float>", genparMass);
  tree_->Branch("GenPar_ID", "vector<int>", genparID);
}


GenJetProducer::~GenJetProducer()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
GenJetProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<reco::GenJetCollection>  GenJetHdl;
  edm::Handle<reco::GenParticleCollection>  GenParHdl;
  iEvent.getByToken(GenJetTok_, GenJetHdl); 
  iEvent.getByToken(GenParTok_, GenParHdl); 

  ProdGenJet(GenJetHdl);
  ProdGenPar(GenParHdl);
  tree_->Fill();
}


// ===  FUNCTION  ============================================================
//         Name:  GenJetProducer::ProdGenJet
//  Description:  
// ===========================================================================
bool GenJetProducer::ProdGenJet(edm::Handle<reco::GenJetCollection>  GenJetHdl)
{

  genjetPt->clear();
  genjetEta->clear();
  genjetPhi->clear();
  genjetMass->clear();

  for(unsigned int i=0; i < GenJetHdl->size(); ++i)
  {
    reco::GenJet j = GenJetHdl->at(i);
    genjetPt->push_back(j.pt());
    genjetEta->push_back(j.eta());
    genjetPhi->push_back(j.phi());
    genjetMass->push_back(j.mass());
  }
  return true;
}       // -----  end of function GenJetProducer::ProdGenJet  -----


// ===  FUNCTION  ============================================================
//         Name:  GenJetProducer::ProdGenPar
//  Description:  
// ===========================================================================
bool GenJetProducer::ProdGenPar(edm::Handle<reco::GenParticleCollection> GenParHdl)
{
  genparPt->clear();
  genparEta->clear();
  genparPhi->clear();
  genparMass->clear();
  genparID->clear();

  for (const reco::GenParticle &g : *GenParHdl) {
    if (g.isHardProcess() ||g.fromHardProcessFinalState())
    {
      genparPt->push_back(g.pt());
      genparEta->push_back(g.eta());
      genparPhi->push_back(g.phi());
      genparMass->push_back(g.mass());
      genparID->push_back(g.pdgId());
    }
  }

    std::cout <<  "End of eventi"<< std::endl;

  return true;
}       // -----  end of function GenJetProducer::ProdGenPar  -----

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
GenJetProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
GenJetProducer::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
GenJetProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
GenJetProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
GenJetProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
GenJetProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenJetProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenJetProducer);
