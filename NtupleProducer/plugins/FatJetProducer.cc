// -*- C++ -*-
//
// Package:    FastPUPPI/FatJet
// Class:      FatJet
// 
/**\class FatJet FatJet.cc FastPUPPI/FatJet/plugins/FatJet.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  zhenbin wu
//         Created:  Sat, 16 Sep 2017 00:51:12 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "TTree.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class FatJet : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit FatJet(const edm::ParameterSet&);
    ~FatJet();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    // ----------member data ---------------------------
    TTree *tree_;
    edm::EDGetTokenT<reco::PFJetCollection>   FatJetTok_;
    edm::EDGetTokenT<edm::ValueMap<float> >   Tau1Tok_;
    edm::EDGetTokenT<edm::ValueMap<float> >   Tau2Tok_;
    edm::EDGetTokenT<edm::ValueMap<float> >   Tau3Tok_;
    edm::EDGetTokenT<edm::ValueMap<float> >   Tau4Tok_;
    float pt;
    std::vector<float> *fatjetPt;
    std::vector<float> *fatjetEta;
    std::vector<float> *fatjetPhi;
    std::vector<float> *fatjetMass;
    std::vector<float> *tau1;
    std::vector<float> *tau2;
    std::vector<float> *tau3;
    std::vector<float> *tau4;
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
FatJet::FatJet(const edm::ParameterSet& iConfig):
  FatJetTok_( consumes<reco::PFJetCollection> (iConfig.getParameter<edm::InputTag>("FatJetTag")) ),
  Tau1Tok_( consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>("Tau1Tag") )),
  Tau2Tok_( consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>("Tau2Tag")) ),
  Tau3Tok_( consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>("Tau3Tag")) ),
  Tau4Tok_( consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>("Tau4Tag")) )
{
  fatjetPt   = new std::vector<float>();
  fatjetEta  = new std::vector<float>();
  fatjetPhi  = new std::vector<float>();
  fatjetMass = new std::vector<float>();
  tau1 = new std::vector<float>();
  tau2 = new std::vector<float>();
  tau3 = new std::vector<float>();
  tau4 = new std::vector<float>();


  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree","tree");
  tree_->Branch("mc_pt", &pt, "mc_pt/F");
  tree_->Branch("FatJet_pt", "vector<float>", fatjetPt);
  tree_->Branch("FatJet_eta", "vector<float>", fatjetEta);
  tree_->Branch("FatJet_phi", "vector<float>", fatjetPhi);
  tree_->Branch("FatJet_mass", "vector<float>", fatjetMass);
  tree_->Branch("Tau1", "vector<float>", tau1);
  tree_->Branch("Tau2", "vector<float>", tau2);
  tree_->Branch("Tau3", "vector<float>", tau3);
  tree_->Branch("Tau4", "vector<float>", tau4);
}


FatJet::~FatJet()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
FatJet::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<reco::PFJetCollection>  FatJetHdl;
  edm::Handle<edm::ValueMap<float> >  Tau1Hdl;
  edm::Handle<edm::ValueMap<float> >  Tau2Hdl;
  edm::Handle<edm::ValueMap<float> >  Tau3Hdl;
  edm::Handle<edm::ValueMap<float> >  Tau4Hdl;

  //pt = 10;
  iEvent.getByToken(FatJetTok_, FatJetHdl); 
  iEvent.getByToken(Tau1Tok_, Tau1Hdl); 
  iEvent.getByToken(Tau2Tok_, Tau2Hdl); 
  iEvent.getByToken(Tau3Tok_, Tau3Hdl); 
  iEvent.getByToken(Tau4Tok_, Tau4Hdl); 

  fatjetPt->clear();
  fatjetEta->clear();
  fatjetPhi->clear();
  fatjetMass->clear();
  tau1->clear();
  tau2->clear();
  tau3->clear();
  tau4->clear();
  for(unsigned int i=0; i < FatJetHdl->size(); ++i)
  {
    reco::PFJet j = FatJetHdl->at(i);
    reco::PFJetRef jref(FatJetHdl, i);
    //tau1 = j.userFloat(NJettiLabel_+":tau1");
    fatjetPt->push_back(j.pt());
    fatjetEta->push_back(j.eta());
    fatjetPhi->push_back(j.phi());
    fatjetMass->push_back(j.mass());
    tau1->push_back((*Tau1Hdl)[jref]);
    tau2->push_back((*Tau2Hdl)[jref]);
    tau3->push_back((*Tau3Hdl)[jref]);
    tau4->push_back((*Tau4Hdl)[jref]);
  }
  tree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
FatJet::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
FatJet::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FatJet::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FatJet);
