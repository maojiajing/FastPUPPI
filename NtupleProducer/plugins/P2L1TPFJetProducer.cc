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
  mergingWTA ( iConfig.getParameter<bool>   ("mergingWTA")),
  N2Tile ( iConfig.getParameter<bool>   ("N2Tile")),
  N2Group ( iConfig.getParameter<bool>   ("N2Group"))
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
  jetCollection.reset();
  fjInputs.clear();
  fjJets.clear();
  gInputs.clear();
  gJets.clear();
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Getting the inputs ~~~~~
  iEvent.getByToken(ParTok_, ParHdl); 
  GetInput(iEvent, iSetup);

  if (N2Tile) SplitN2Tile(10, 10);
  else if (N2Group) SplitN2Groups();
  else
  {
    gInputs.push_back(fjInputs);
    gJets.push_back(fjJets);
  }
  std::cout << " Run with N2Tile " << N2Tile <<" N2Group "<< N2Group <<" MetricKt " << metricKt <<" metricdR " 
    << metricdR <<" mergingE " << mergingE <<" mergingWTA " << mergingWTA << " inputs " << ParHdl->size() << " Ngroups " << gInputs.size() << std::endl;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Run clustering ~~~~~
  for (unsigned i = 0; i < gInputs.size(); ++i)
  {
    //std::cout << i<< "---------Initial size " << gInputs.at(i).size()<< std::endl;
    //std::cout << gInputs.at(i).size()  << " " << gJets.at(i).size()<< std::endl;
    Iterations(gInputs.at(i), gJets.at(i));
  }


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Store the output ~~~~~
  RemergeJets();
  WriteJetCollection(iSetup);
  iEvent.put(std::move(jetCollection));
}


// ===  FUNCTION  ============================================================
//         Name:  P2L1TPFJetProducer::SplitN2Tile
//  Description:  /* cursor */
// ===========================================================================
void P2L1TPFJetProducer::SplitN2Tile(int Neta, int Nphi) 
{

  std::vector<float> etabounds, phibounds;
  float deta = 10/Neta;
  for (int i = 0; i < Neta; ++i)
  {
    etabounds.push_back(-5 + i*deta);
  }

  float dphi = 2*3.14/Neta;
  for (int i = 0; i < Neta; ++i)
  {
    phibounds.push_back(-3.14 + i*dphi);
  }

  for (int i = 0; i < Neta*Nphi; ++i)
  {
    gJets.push_back(std::vector<TLorentzVector>());
    gInputs.push_back(std::vector<TLorentzVector>());
  }

  for (unsigned i = 0; i < fjInputs.size(); i++)
  {
    int etaidx=0;
    for (unsigned j = 0; j < etabounds.size()-1; ++j)
    {
      if(fjInputs.at(i).Eta() > etabounds[j]  && fjInputs.at(i).Eta() < etabounds[j+1]  )
      {
        etaidx = j;
        break;
      }
    }

    int phiidx=0;
    for (unsigned k = 0; k < phibounds.size()-1; ++k)
    {
      if(fjInputs.at(i).Phi() > phibounds[k]  && fjInputs.at(i).Phi() < phibounds[k+1]  )
      {
        phiidx = k;
        break;
      }
    }
    gInputs[etaidx*Nphi + phiidx].push_back(fjInputs.at(i));
  }


}       // -----  end of function P2L1TPFJetProducer::SplitN2Tile  -----

// ===  FUNCTION  ============================================================
//         Name:  P2L1TPFJetProducer::SplitN2Groupsv2
//  Description:  
// ===========================================================================
void P2L1TPFJetProducer::SplitN2Groupsv2()
{
  size_t intsize = fjInputs.size();
  short neighours[intsize][intsize] = {{0}};

  for (unsigned i = 0; i < fjInputs.size(); i++)
  {
    // keep i=j for grouping later on
    for (unsigned j = i; j < fjInputs.size()-1; ++j)
    {
      double dEta = fabs( fjInputs[i].Eta() -  fjInputs[j].Eta());
      double dPhi = TVector2::Phi_mpi_pi( fjInputs[i].Phi() -  fjInputs[j].Phi());
      if (dEta < rParam && dPhi < rParam)
        neighours[i][j]=1;
    }
  }


  std::vector<std::set<int> > groupset;
  std::set<int> linkedidx = {0};
  std::vector<TLorentzVector> curgroup;

  do
  {
    for (std::set<int>::const_iterator it = linkedidx.begin(); it != linkedidx.end(); ++it)
    {
      unsigned i = *it;
      for (unsigned j = i; j < fjInputs.size(); ++j)
      {
        // Already read this row
        if (i==j && neighours[i][j] == 0)
          break;

        if ( neighours[i][j] == 1) 
        {
          linkedidx.insert(j);
          neighours[i][j] = 0;
        }
      }
    }

    groupset.push_back(linkedidx);
    linkedidx.clear();
    for (unsigned i = 0; i < fjInputs.size(); ++i)
    {
        if ( neighours[i][i] == 1) 
        {
          linkedidx.insert(i);
          break;
        }
    }
    if (linkedidx.empty()) break;
  }
  while (true);

  for(auto idx : groupset)
  {
    std::vector<TLorentzVector> temp;
    for(auto i : idx)
    {
      temp.push_back(fjInputs.at(i));
    }
    gInputs.push_back(temp);
    gJets.push_back(std::vector<TLorentzVector>());
  }


}       // -----  end of function P2L1TPFJetProducer::SplitN2Groupsv2  -----


// ===  FUNCTION  ============================================================
//         Name:  P2L1TPFJetProducer::SplitN2Groups
//  Description:  
// ===========================================================================
void P2L1TPFJetProducer::SplitN2Groups()
{
  size_t intsize = fjInputs.size();
  short neighours[intsize][intsize] = {{0}};

  for (unsigned i = 0; i < fjInputs.size(); i++)
  {
    // keep i=j for grouping later on
    for (unsigned j = i; j < fjInputs.size()-1; ++j)
    {
      double dEta = fabs( fjInputs[i].Eta() -  fjInputs[j].Eta());
      double dPhi = TVector2::Phi_mpi_pi( fjInputs[i].Phi() -  fjInputs[j].Phi());
      if (dEta < rParam && dPhi < rParam)
        neighours[i][j]=1;
    }
  }


  std::vector<std::set<int> > groupset;
  std::set<int> linkedidx = {0};
  std::vector<TLorentzVector> curgroup;

  do
  {
    for (std::set<int>::const_iterator it = linkedidx.begin(); it != linkedidx.end(); ++it)
    {
      unsigned i = *it;
      for (unsigned j = i; j < fjInputs.size(); ++j)
      {
        // Already read this row
        if (i==j && neighours[i][j] == 0)
          break;

        if ( neighours[i][j] == 1) 
        {
          linkedidx.insert(j);
          neighours[i][j] = 0;
        }
      }
    }

    groupset.push_back(linkedidx);
    linkedidx.clear();
    for (unsigned i = 0; i < fjInputs.size(); ++i)
    {
        if ( neighours[i][i] == 1) 
        {
          linkedidx.insert(i);
          break;
        }
    }
    if (linkedidx.empty()) break;
  }
  while (true);

  for(auto idx : groupset)
  {
    std::vector<TLorentzVector> temp;
    for(auto i : idx)
    {
      temp.push_back(fjInputs.at(i));
    }
    gInputs.push_back(temp);
    gJets.push_back(std::vector<TLorentzVector>());
  }


}       // -----  end of function P2L1TPFJetProducer::SplitN2Groups  -----


// ===  FUNCTION  ============================================================
//         Name:  P2L1TPFJetProducer::RemergeJets
//  Description:  
// ===========================================================================
bool P2L1TPFJetProducer::RemergeJets()
{
  for(auto gj : gJets)
  {
    fjJets.insert(fjJets.end(), gj.begin(), gj.end());
  }
  return true;
}       // -----  end of function P2L1TPFJetProducer::RemergeJets  -----

// ===  FUNCTION  ============================================================
//         Name:  P2L1TPFJetProducer::WriteJetCollection
//  Description:  
// ===========================================================================
void P2L1TPFJetProducer::WriteJetCollection(const edm::EventSetup& iSetup)
{
  jetCollection = std::make_unique<reco::PFJetCollection>();
  std::vector<math::XYZTLorentzVector> p4_Jets;
  for(auto g : fjJets)
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
    //if (ParHdl->at(i).pt() < 2)
      //continue;
    TLorentzVector temp( ParHdl->at(i).px(), ParHdl->at(i).py(),
        ParHdl->at(i).pz(), ParHdl->at(i).energy() );
    fjInputs.push_back(temp);
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
bool P2L1TPFJetProducer::GetCombinationsKt( std::vector<TLorentzVector> &fjInputs_, std::map<float, std::pair<unsigned, unsigned> > &Valcomb)
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
bool P2L1TPFJetProducer::GetCombinationsdR( std::vector<TLorentzVector> &fjInputs_, std::map<float, std::pair<unsigned, unsigned> > &Valcomb)
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
bool P2L1TPFJetProducer::Iterations( std::vector<TLorentzVector> &fjInputs_, std::vector<TLorentzVector> &fjJets_)
{
  if(fjInputs_.empty())
    return false;
  std::map<float, std::pair<unsigned, unsigned> > Valcomb;
  while (true)
  {
    Valcomb.clear();
    if (metricKt)
      GetCombinationsKt(fjInputs_, Valcomb);
    if (metricdR)
      GetCombinationsdR(fjInputs_, Valcomb);
    //float minval = Valcomb.begin()->first;
    if(Valcomb.empty()) break;
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
      {
        newObj = MergingWTA(fjInputs_.at(minpair.first), fjInputs_.at(minpair.second));
      }
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
