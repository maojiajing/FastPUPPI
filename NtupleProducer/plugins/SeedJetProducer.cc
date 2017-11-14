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

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "TLorentzVector.h"
//
// class declaration
//

typedef std::vector< reco::PFCandidate >               PFOutputCollection;
class SeedJetProducer : public edm::stream::EDProducer<> {
   public:
      explicit SeedJetProducer(const edm::ParameterSet&);
      ~SeedJetProducer();

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
      bool SeedCone(std::vector<unsigned> seeds);
      std::vector<unsigned> GetSeed(const int pt, const int &id);
      TTree *tree_;
      edm::EDGetTokenT<PFOutputCollection> ParTok_;
      edm::Handle<PFOutputCollection> ParHdl;
      static const Int_t Max = 1500;
      Int_t npar;
      Int_t nleft;
      Float_t pt[Max];
      Float_t eta[Max];
      Float_t phi[Max];
      Float_t mass[Max];

      // Clusting
      double seedPt;
      unsigned seedID;
      double cone;



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
  ParTok_( consumes<reco::PFCandidateCollection> (iConfig.getParameter<edm::InputTag>("ParTag")) ),
  seedPt(iConfig.getParameter<double>       ("seedPtCut")),
  seedID(iConfig.getParameter<unsigned>       ("seedID")),
  cone(iConfig.getParameter<double>       ("cone"))
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
  tree_->Branch("nleft",&nleft,"nleft/I");
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


  iEvent.getByToken(ParTok_, ParHdl); 
  SeedCone(GetSeed(seedPt, seedID));


  //ProdPar(ParHdl);
  tree_->Fill();
 
}


// ===  FUNCTION  ============================================================
//         Name:  SeedJetProducer::GetSeed
//  Description:  
// ===========================================================================
std::vector<unsigned> SeedJetProducer::GetSeed(const int pt, const int &id)
{
  unsigned n= 0;
  std::vector<unsigned> idx;
  std::map< float, unsigned > ptmap;

  for (const reco::PFCandidate &g : *ParHdl) {
    if (g.pt() >= pt)
    {
      if (id == 999 || g.particleId() == id)
      {
        ptmap[g.pt()] = n;
      }
    }
    n++;
  }


  for(std::map<float, unsigned>::reverse_iterator i=ptmap.rbegin(); i!=ptmap.rend(); ++i)
  {
    idx.push_back(i->second);
  }
  return idx;
}       // -----  end of function SeedJetProducer::GetSeed  -----


// ===  FUNCTION  ============================================================
//         Name:  SeedJetProducer::SeedCone
//  Description:  
// ===========================================================================
bool SeedJetProducer::SeedCone(std::vector<unsigned> seeds)
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
      reco::PFCandidate p = ParHdl->at(j);
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


  int n= 0;
  npar = rejets.size();
  for(auto g : rejets)
  {
    pt[n] = g.Pt();
    eta[n] = g.Eta();
    phi[n] = g.Phi();
    mass[n] = g.M();
    if (n++ > Max) break;
  }
  nleft = count(taken.begin(), taken.end(), false);

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
