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
    static const Int_t Max = 1500;
    Int_t npar;
    Float_t pt[Max];
    Float_t eta[Max];
    Float_t phi[Max];
    Float_t mass[Max];
    Float_t tau1[Max];
    Float_t tau2[Max];
    Float_t tau3[Max];
    Float_t tau4[Max];
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
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree","tree");

  tree_->Branch("npar",&npar,"npar/I");
  tree_->Branch("pt",pt,"pt[npar]/F");
  tree_->Branch("eta",eta,"eta[npar]/F");
  tree_->Branch("phi",phi,"phi[npar]/F");
  tree_->Branch("mass",mass,"mass[npar]/F");
  tree_->Branch("Tau1",tau1,"Tau1[npar]/F");
  tree_->Branch("Tau2",tau2,"Tau2[npar]/F");
  tree_->Branch("Tau3",tau3,"Tau3[npar]/F");
  tree_->Branch("Tau4",tau4,"Tau4[npar]/F");
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

  npar = FatJetHdl->size();

  for(unsigned int i=0; i < FatJetHdl->size(); ++i)
  {
    reco::PFJet j = FatJetHdl->at(i);
    reco::PFJetRef jref(FatJetHdl, i);
    //tau1 = j.userFloat(NJettiLabel_+":tau1");
    pt[i] = j.pt();
    eta[i] = j.eta();
    phi[i] = j.phi();
    mass[i] = j.mass();
    tau1[i] = (*Tau1Hdl)[jref];
    tau2[i] = (*Tau2Hdl)[jref];
    tau3[i] = (*Tau3Hdl)[jref];
    tau4[i] = (*Tau4Hdl)[jref];
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
