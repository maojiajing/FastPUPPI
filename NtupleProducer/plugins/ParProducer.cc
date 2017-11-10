// -*- C++ -*-
//
// Package:    FastPUPPI/ParProducer
// Class:      ParProducer
// 
/**\class ParProducer ParProducer.cc FastPUPPI/ParProducer/plugins/ParProducer.cc

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

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
//
// class declaration
//

typedef std::vector< reco::PFCandidate >               PFOutputCollection;
class ParProducer : public edm::stream::EDProducer<> {
   public:
      explicit ParProducer(const edm::ParameterSet&);
      ~ParProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      bool ProdPar(edm::Handle<PFOutputCollection> ParHdl);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      TTree *tree_;
      edm::EDGetTokenT<PFOutputCollection> ParTok_;

      static const Int_t Max = 1500;
      Int_t npar;
      Int_t id[Max];
      Int_t charge[Max];
      Float_t pt[Max];
      Float_t eta[Max];
      Float_t phi[Max];
      Float_t mass[Max];
      Float_t ecal[Max];
      Float_t hcal[Max];





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
ParProducer::ParProducer(const edm::ParameterSet& iConfig):
  ParTok_( consumes<reco::PFCandidateCollection> (iConfig.getParameter<edm::InputTag>("ParTag")) )
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
  tree_->Branch("npar",&npar,"npar/I");

  tree_->Branch("pt",pt,"pt[npar]/F");
  tree_->Branch("eta",eta,"eta[npar]/F");
  tree_->Branch("phi",phi,"phi[npar]/F");
  tree_->Branch("mass",mass,"mass[npar]/F");

  
}


ParProducer::~ParProducer()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
ParProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
/* This is an event example
   //Read 'ExampleData' from the Event
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);

   //Use the ExampleData to create an ExampleData2 which 
   // is put into the Event
   iEvent.put(std::make_unique<ExampleData2>(*pIn));
*/

/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/
  edm::Handle<PFOutputCollection>  ParHdl;
  iEvent.getByToken(ParTok_, ParHdl); 

  ProdPar(ParHdl);
  tree_->Fill();
 
}

// ===  FUNCTION  ============================================================
//         Name:  ParProducer::ProdPar
//  Description:  
// ===========================================================================
bool ParProducer::ProdPar(edm::Handle<PFOutputCollection> ParHdl)
{
  int n= 0;
  npar = ParHdl->size();
  
  for (const reco::PFCandidate &g : *ParHdl) {
    pt[n] = g.pt();
    eta[n] = g.eta();
    phi[n] = g.phi();
    mass[n] = g.mass();

    id[n] = g.particleId();
    charge[n]=g.charge();
    ecal[n] = g.ecalEnergy();
    hcal[n] = g.hcalEnergy();

    if (n++ > Max) break;
  }

  return true;
}       // -----  end of function ParProducer::ProdPar  -----

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
ParProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
ParProducer::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
ParProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
ParProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
ParProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
ParProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ParProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ParProducer);
