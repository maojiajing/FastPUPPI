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
    static const Int_t Max = 1500;
    Int_t njets;
    Int_t npars;
    Float_t genjetPt[Max];
    Float_t genjetEta[Max];
    Float_t genjetPhi[Max];
    Float_t genjetMass[Max];
    Float_t genparPt[Max];
    Float_t genparEta[Max];
    Float_t genparPhi[Max];
    Float_t genparMass[Max];
    Int_t genparID[Max];
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
  tree_->Branch("njets",&njets,"njets/I");
  tree_->Branch("GenJet_pt",genjetPt,"genjetPt[njets]/F");
  tree_->Branch("GenJet_eta",genjetEta,"genjetEta[njets]/F");
  tree_->Branch("GenJet_phi",genjetPhi,"genjetPhi[njets]/F");
  tree_->Branch("GenJet_mass",genjetMass,"genjetMass[njets]/F");

  tree_->Branch("npars",&npars,"npars/I");
  tree_->Branch("GenPar_pt",genparPt,"genparPt[npars]/F");
  tree_->Branch("GenPar_eta",genparEta,"genparEta[npars]/F");
  tree_->Branch("GenPar_phi",genparPhi,"genparPhi[npars]/F");
  tree_->Branch("GenPar_mass",genparMass,"genparMass[npars]/F");
  tree_->Branch("GenPar_ID",genparID,"genparID[npars]/I");
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

  njets = GenJetHdl->size();
  for(unsigned int i=0; i < GenJetHdl->size(); ++i)
  {
    reco::GenJet j = GenJetHdl->at(i);
    genjetPt[i]  = j.pt();
    genjetEta[i]  = j.eta();
    genjetPhi[i]  = j.phi();
    genjetMass[i]  = j.mass();
  }
  return true;
}       // -----  end of function GenJetProducer::ProdGenJet  -----


// ===  FUNCTION  ============================================================
//         Name:  GenJetProducer::ProdGenPar
//  Description:  
// ===========================================================================
bool GenJetProducer::ProdGenPar(edm::Handle<reco::GenParticleCollection> GenParHdl)
{
  int i = 0;
  for (const reco::GenParticle &g : *GenParHdl) {
    if (g.isHardProcess() ||g.fromHardProcessFinalState())
    {
      genparPt[i] = g.pt();
      genparEta[i] = g.eta();
      genparPhi[i] = g.phi();
      genparMass[i] = g.mass();
      genparID[i] = g.pdgId();
      i++;
    }
  }
  npars = i;

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
