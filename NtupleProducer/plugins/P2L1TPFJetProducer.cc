// -*- C++ -*-
//
// Package:    FastPUPPI/P2L1TPFJetProducer
// Class:      P2L1TPFJetProducer
// 
/**\class P2L1TPFJetProducer P2L1TPFJetProducer.cc FastPUPPI/P2L1TPFJetProducer/plugins/P2L1TPFJetProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Zhenbin Wu
//         Created:  Wed, 13 Mar 2019 02:48:14 GMT
//
//

#include "P2L1TPFJetProducer.h"

//
// constructors and destructor
//
P2L1TPFJetProducer::P2L1TPFJetProducer(const edm::ParameterSet& iConfig):
  ParTok_( consumes<l1t::PFCandidateCollection> (iConfig.getParameter<edm::InputTag>("src")) ),
  rParam     ( iConfig.getParameter<double> ("rParam")),
  jetPtMin   ( iConfig.getParameter<double> ("jetPtMin")),
  inputEtMin ( iConfig.getParameter<double> ("inputEtMin")),
  ktsign     ( iConfig.getParameter<int>    ("ktsign")),
  metricKt   ( iConfig.getParameter<bool>   ("metricKt")),
  metricdR   ( iConfig.getParameter<bool>   ("metricdR")),
  mergingE   ( iConfig.getParameter<bool>   ("mergingE")),
  mergingWTA ( iConfig.getParameter<bool>   ("mergingWTA"))
{
  produces<reco::PFJetCollection>();
}


P2L1TPFJetProducer::~P2L1TPFJetProducer()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
P2L1TPFJetProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Getting the inputs ~~~~~
  iEvent.getByToken(ParTok_, ParHdl); 
  GetInput(iEvent, iSetup);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Run clustering ~~~~~
  Iterations();

  WriteJetCollection(iSetup);


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Store the output ~~~~~
  iEvent.put(std::move(jetCollection));
}


// ===  FUNCTION  ============================================================
//         Name:  P2L1TPFJetProducer::WriteJetCollection
//  Description:  
// ===========================================================================
void P2L1TPFJetProducer::WriteJetCollection(const edm::EventSetup& iSetup)
{
  jetCollection = std::make_unique<reco::PFJetCollection>();
  std::vector<math::XYZTLorentzVector> p4_Jets;
  for(auto g : fjJets_)
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
}       // -----  end of function P2L1TPFJetProducer::WriteJetCollection  -----

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
P2L1TPFJetProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
P2L1TPFJetProducer::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
P2L1TPFJetProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
P2L1TPFJetProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
P2L1TPFJetProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
P2L1TPFJetProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ===  FUNCTION  ============================================================
//         Name:  P2L1TPFJetProducer::GetInput
//  Description:  
// ===========================================================================
void P2L1TPFJetProducer::GetInput(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Convert to local TLV ~~~~~
  for (size_t i = 0; i < ParHdl->size(); ++i) {
    TLorentzVector temp( ParHdl->at(i).px(), ParHdl->at(i).py(),
        ParHdl->at(i).pz(), ParHdl->at(i).energy() );
    fjInputs_.push_back(temp);
  }

  //for (size_t i = 0; i < fjInputs_.size(); ++i)
  //{
    //std::cout << i <<  " pt " << fjInputs_.at(i).Pt() << std::endl;
  //}

}       // -----  end of function P2L1TPFJetProducer::GetInput  -----

// ===  FUNCTION  ============================================================
//         Name:  P2L1TPFJetProducer::GetCombinationsKt
//  Description:  
// ===========================================================================
bool P2L1TPFJetProducer::GetCombinationsKt()
{
  if(!metricKt) return false;
  
  for (unsigned i = 0; i < fjInputs_.size(); ++i)
  {
      for (unsigned j = i; j < fjInputs_.size(); ++j)
      {
        float val= -999;
        if (fjInputs_.at(i).DeltaR(fjInputs_.at(j)) > rParam) continue;
        val = Metric_AK(fjInputs_.at(i), fjInputs_.at(j), i==j);

        //combVal[std::make_pair(i, j)]= val;
        Valcomb[val] = std::make_pair(i, j);
      }
  }
  return true;
}       // -----  end of function P2L1TPFJetProducer::GetCombinationsKt  -----


// ===  FUNCTION  ============================================================
//         Name:  P2L1TPFJetProducer::GetCombinationsdR
//  Description:  
// ===========================================================================
bool P2L1TPFJetProducer::GetCombinationsdR()
{
  if (!metricdR) return false;

  std::vector<float> pts;
  for(auto i : fjInputs_)
  {
    pts.push_back(i.Pt());
  }

  unsigned maxptidx= std::distance(pts.begin(), std::max_element(pts.begin(), pts.end()));
  
  for (size_t i = 0; i < fjInputs_.size(); ++i)
  {
    float val= -999;
    if (fjInputs_.at(i).DeltaR(fjInputs_.at(maxptidx)) > rParam) continue;
    val = Metric_dR(fjInputs_.at(i), fjInputs_.at(maxptidx), i==maxptidx);
    Valcomb[val] = std::make_pair(i, maxptidx);
  };

  return true;
}       // -----  end of function P2L1TPFJetProducer::GetCombinationsdR  -----


// ===  FUNCTION  ============================================================
//         Name:  P2L1TPFJetProducer::Iterations
//  Description:  
// ===========================================================================
bool P2L1TPFJetProducer::Iterations()
{
  while (true)
  {
    Valcomb.clear();
    if (metricKt)
      GetCombinationsKt();
    if (metricdR)
      GetCombinationsdR();
    //float minval = Valcomb.begin()->first;
    std::pair<unsigned, unsigned > minpair = Valcomb.begin()->second;
    if (minpair.first == minpair.second) // Jet is ready
    {
      TLorentzVector njet =fjInputs_.at(minpair.second);
      fjJets_.push_back(njet);
      fjInputs_.erase(fjInputs_.begin()+minpair.second);
    }
    else{
      TLorentzVector newObj;
      if (mergingE)
        newObj = MergingE(fjInputs_.at(minpair.first), fjInputs_.at(minpair.second));
      if (mergingWTA)
        newObj = MergingWTA(fjInputs_.at(minpair.first), fjInputs_.at(minpair.second));
      fjInputs_.erase(fjInputs_.begin()+minpair.first);
      fjInputs_.erase(fjInputs_.begin()+minpair.second);
      fjInputs_.push_back(newObj);
    }
    if (fjInputs_.size() == 0) break;
  }
  
  return true;
}       // -----  end of function P2L1TPFJetProducer::Iterations  -----

// ===  FUNCTION  ============================================================
//         Name:  P2L1TPFJetProducer::Metric_AK
//  Description:  [> cursor <]
// ===========================================================================
inline float P2L1TPFJetProducer::Metric_AK(TLorentzVector p1, TLorentzVector p2, bool sameobj) const
{

  if(sameobj)
  {
    assert(p1.Pt() == p2.Pt());
    return pow(p1.Pt(), ktsign*2);
  }

  float ptterm = std::min(pow(p1.Pt(), ktsign*2), pow(p2.Pt(), ktsign*2));

  float dR = p1.DeltaR(p2);
  return ptterm* dR*dR/(rParam *rParam);
}       // -----  end of function P2L1TPFJetProducer::Metric_AK  -----

// ===  FUNCTION  ============================================================
//         Name:  P2L1TPFJetProducer::MergingE
//  Description:  [> cursor <]
// ===========================================================================
inline TLorentzVector P2L1TPFJetProducer::MergingE(TLorentzVector &p1, TLorentzVector &p2) const
{
  return p1+p2;
}       // -----  end of function P2L1TPFJetProducer::MergingE  -----

// ===  FUNCTION  ============================================================
//         Name:  P2L1TPFJetProducer::MergingWTA
//  Description:  
// ===========================================================================
inline TLorentzVector P2L1TPFJetProducer::MergingWTA(TLorentzVector &p1, TLorentzVector &p2) const
{
  TLorentzVector winner= p1.Pt() > p2.Pt() ? p1 : p2;
  winner.SetPerp(p1.Pt() + p2.Pt());
  return winner;
}       // -----  end of function P2L1TPFJetProducer::MergingWTA  -----


// ===  FUNCTION  ============================================================
//         Name:  P2L1TPFJetProducer::Metric_dR
//  Description:  
// ===========================================================================
inline float P2L1TPFJetProducer::Metric_dR(TLorentzVector p1, TLorentzVector p2, bool sameobj) const
{
  if (sameobj) return pow(rParam, 2);
  return pow(p1.DeltaR(p2),2);
}       // -----  end of function P2L1TPFJetProducer::Metric_dR  -----

 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
P2L1TPFJetProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(P2L1TPFJetProducer);
