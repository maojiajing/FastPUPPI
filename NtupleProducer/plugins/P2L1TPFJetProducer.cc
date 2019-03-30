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
  ParTok_    ( consumes<l1t::PFCandidateCollection> ( iConfig.getParameter<edm::InputTag> ( "src")) ),
  rParam     ( iConfig.getParameter<double>         ( "rParam")),
  groupR     ( iConfig.getParameter<double>         ( "groupR")),
  jetPtMin   ( iConfig.getParameter<double>         ( "jetPtMin")),
  inputEtMin ( iConfig.getParameter<double>         ( "inputEtMin")),
  ktsign     ( iConfig.getParameter<int>            ( "ktsign")),
  metricKt   ( iConfig.getParameter<bool>           ( "metricKt")),
  metricdR   ( iConfig.getParameter<bool>           ( "metricdR")),
  mergingE   ( iConfig.getParameter<bool>           ( "mergingE")),
  mergingWTA ( iConfig.getParameter<bool>           ( "mergingWTA")),
  resumWTA   ( iConfig.getParameter<bool>           ( "resumWTA")),
  N2Tile     ( iConfig.getParameter<bool>           ( "N2Tile")),
  N2Group    ( iConfig.getParameter<bool>           ( "N2Group"))
{
  debug=true;
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
  else if (N2Group) SplitN2Groupsv2();
  else
  {
    gInputs.push_back(fjInputs);
    gJets.push_back(fjJets);
  }

  if (debug)
    std::cout << " Run with N2Tile " << N2Tile <<" N2Group "<< N2Group <<" MetricKt " << metricKt <<" metricdR " 
      << metricdR <<" mergingE " << mergingE <<" mergingWTA " << mergingWTA << " inputs " << ParHdl->size() << " Ngroups " << gInputs.size() 
      << " groupdR "<< groupR<< std::endl;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Run clustering ~~~~~
  for (unsigned i = 0; i < gInputs.size(); ++i)
  {
    //std::cout << i<< "---------Initial size " << gInputs.at(i).size()<< std::endl;
    //std::cout << gInputs.at(i).size()  << " " << gJets.at(i).size()<< std::endl;
    Iterations(gInputs.at(i), gJets.at(i));
    //std::cout<<"Run to \033[0;31m"<<__func__<<"\033[0m at \033[1;36m"<< __FILE__<<"\033[0m, line \033[0;34m"<< __LINE__<<"\033[0m"<< std::endl; 
  }


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Store the output ~~~~~
  RemergeJets();
  WriteJetCollection(iSetup);
  iEvent.put(std::move(jetCollection));

  //std::vector<unsigned> grpsizes;
  //for(auto i : gInputs)
  //{
    //grpsizes.push_back(i.size());
  //}
  //auto outg = std::make_unique<nanoaod::FlatTable>(grpsizes.size(), "GroupsInfo", true); 
  //outg->addColumn<unsigned>("GroupSize", grpsizes,  "pt of matched gen jet", nanoaod::FlatTable::UInt8Column);
  //outg->addColumnValue<unsigned>("InputSize",  fjInputs.size(), "genMet pt", nanoaod::FlatTable::UInt8Column);
  //iEvent.put(std::move(outg), "GroupsInfo");
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
//  Description:  2nd version, using arrays as in FPGA
// ===========================================================================
void P2L1TPFJetProducer::SplitN2Groupsv2()
{
#define MAXSIZE 300
  const size_t intsize = fjInputs.size();
  const size_t grpsize = intsize;
  if (intsize > MAXSIZE)
    std::cout << " WARNING: Input array is too big for  MAXSIZE " << MAXSIZE   << std::endl;
  std::bitset<MAXSIZE> neighours[intsize] = {0};
  std::bitset<MAXSIZE> groups[grpsize] = {0};
  grpInputMap.clear();


  for (unsigned i = 0; i < intsize; i++)
  {
    neighours[i] = 0;
    for (unsigned j = 0; j < intsize; ++j)
    {
      //if (fjInputs[i].DeltaR(fjInputs[j]) < rParam*2)
      if (fjInputs[i].DeltaR(fjInputs[j]) <= groupR)
      {
        neighours[i].set(j);
      }
    }
  }

  for (unsigned i = 0; i < intsize; i++)
  {
    //std::cout << "input "  << i <<"  bit " << neighours[i]<< std::endl;
    //if (neighours[i].none())
      //std::cout << " Empty ? " << i <<"  pt " << fjInputs.at(i).Pt() <<" eta "<< fjInputs.at(i).Eta() <<"  phi "<< fjInputs.at(i).Phi()<< std::endl;
    for (unsigned j = 0; j < grpsize; ++j)
    {
      if ((groups[j] & neighours[i]).any())
      {
        groups[j] |= neighours[i];
        break;
      }
      if(groups[j].none())
      {
        groups[j] = neighours[i];
        break;
      }
      // Run out of groups
      if (j == intsize)
        groups[j] |= neighours[i] ;
    }
  }

  for (unsigned i = 0; i < grpsize; i++)
  {
    //std::cout << " org group " << i << " " << groups[i]<< std::endl;
    for (unsigned j = i+1; j < grpsize; ++j)
    {
      if ((groups[j] & groups[i]).any())
      {
        groups[i] |= groups[j];
        groups[j].reset();
      }
    }
  }

  for (unsigned i = 0; i < grpsize; ++i)
  {
    grpInputMap[i]= std::vector<unsigned>();
    if (groups[i] == 0) continue;
    //std::cout << " final group " << i << " " << groups[i] << std::endl;
    std::vector<TLorentzVector> temp;
    for (unsigned j = 0; j < intsize; ++j)
    {
      if (groups[i].test(j))
      {
        //std::cout << i  << ", " << j <<", " << fjInputs.at(j).Pt() <<", " << fjInputs.at(j).Eta() <<", " << fjInputs.at(j).Phi() << std::endl;
        //std::cout << " Group " << i  << " input " << j <<" pt " << fjInputs.at(j).Pt() <<"  eta " << fjInputs.at(j).Eta() <<" Phi " << fjInputs.at(j).Phi() << std::endl;
        temp.push_back(fjInputs.at(j));
        grpInputMap[i].push_back(j);
      }
    }
    gInputs.push_back(temp);
    gJets.push_back(std::vector<TLorentzVector>());
  }

  unsigned total=0;
  for (unsigned i = 0; i < gInputs.size(); ++i)
  {
    //std::cout << "Group " <<i << " size " << gInputs.at(i).size()  << " index " ;
    //for(auto j : grpInputMap[i])
    //{
      //std::cout << j <<" ";
    //}

    //std::cout << " " << std::endl;
    total += gInputs.at(i).size() ;
  }

  if (total != intsize)
    std::cout << "Totoal size " << total <<" org input " << intsize << std::endl;

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
    for (unsigned j = i; j < fjInputs.size(); ++j)
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

}       // -----  end of function P2L1TPFJetProducer::GetInput  -----

// ===  FUNCTION  ============================================================
//         Name:  P2L1TPFJetProducer::GetCombinationsKt
//  Description:  
// ===========================================================================
bool P2L1TPFJetProducer::GetCombinationsKt( std::vector<TLorentzVector> &fjInputs_, std::map<float, std::pair<unsigned, unsigned> > &Valcomb, 
    std::vector<bool> &idx)
{
  if(!metricKt) return false;
  
  for (unsigned i = 0; i < idx.size(); ++i)
  {
    if(!idx[i]) continue;
      for (unsigned j = i; j < idx.size(); ++j)
      {
        if(!idx[j]) continue;
        float val= -999;
        if (fjInputs_.at(i).DeltaR(fjInputs_.at(j)) > rParam) continue;
        val = Metric_AK(fjInputs_.at(i), fjInputs_.at(j), i==j);

        Valcomb[val] = std::make_pair(i, j);
      }
  }
  return true;
}       // -----  end of function P2L1TPFJetProducer::GetCombinationsKt  -----


// ===  FUNCTION  ============================================================
//         Name:  P2L1TPFJetProducer::GetCombinationsdR
//  Description:  
// ===========================================================================
bool P2L1TPFJetProducer::GetCombinationsdR( std::vector<TLorentzVector> &fjInputs_, std::map<float, std::pair<unsigned, unsigned> > &Valcomb,
    std::vector<bool> &idx)
{
  if (!metricdR) return false;

  std::map<float, unsigned> ptmap;
  for (unsigned i = 0; i < idx.size(); ++i)
  {
    if(idx[i])
      ptmap[fjInputs_[i].Pt()] = i;
  }

  if (ptmap.empty()) return false;
  unsigned maxptidx= ptmap.rbegin()->second;

  //std::cout << " Maxptindx "<< maxptidx <<" max pt " << fjInputs_.at(maxptidx).Pt() << std::endl;
  for (size_t i = 0; i < idx.size(); ++i)
  {
    if (!idx[i]) continue;
    float val= -999;
    if (fjInputs_.at(i).DeltaR(fjInputs_.at(maxptidx)) > rParam) continue;
    val = Metric_dR(fjInputs_.at(i), fjInputs_.at(maxptidx), i==maxptidx);
    Valcomb[val] = std::make_pair(i, maxptidx);
    //std::cout << " valcomb " << val <<" pair " <<i <<" " << maxptidx << std::endl;
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
  std::vector<bool> idx(fjInputs_.size(), true);
  std::map<unsigned, std::pair<unsigned, unsigned> > history;
  std::vector<unsigned> jetidx;
  const unsigned orgsize = fjInputs_.size();

  while (true)
  {
    Valcomb.clear();
    if (metricKt)
      GetCombinationsKt(fjInputs_, Valcomb, idx);
    if (metricdR)
      GetCombinationsdR(fjInputs_, Valcomb, idx);
    if(Valcomb.empty()) break;
    std::pair<unsigned, unsigned > minpair = Valcomb.begin()->second;
    //for(auto v : Valcomb)
    //{
        //std::cout << v.first <<" ";
    //}
    //std::cout << "--------"<< Valcomb.begin()->first << "  "<<fjInputs_.size() <<" "
        //<<" first " << minpair.first <<"  "<< fjInputs_.at(minpair.first).Pt() <<"  "<< fjInputs_.at(minpair.first).Eta() <<"  "<< fjInputs_.at(minpair.first).Phi()
        //<<" second " << minpair.second <<"  "<< fjInputs_.at(minpair.second).Pt() <<"  "<< fjInputs_.at(minpair.second).Eta() <<"  "<< fjInputs_.at(minpair.second).Phi()
        //<< std::endl;
    if (minpair.first == minpair.second) // Jet is ready
    {
      TLorentzVector njet =fjInputs_.at(minpair.second);
      fjJets_.push_back(njet);
      jetidx.push_back(minpair.second);
      idx[minpair.second] = false;
      //std::cout << " newJet  "<< njet.Pt() <<" "<< njet.Eta() << " " << njet.Phi() << " pair index " << minpair.first <<" " << minpair.second<< std::endl;
      history[10000+minpair.first] = std::make_pair(minpair.first, minpair.second);
    }
    else{
      TLorentzVector newObj;
      if (mergingE)
        newObj = MergingE(fjInputs_.at(minpair.first), fjInputs_.at(minpair.second));
      if (mergingWTA)
        newObj = MergingWTA(fjInputs_.at(minpair.first), fjInputs_.at(minpair.second));

      //std::cout << " newObj  "<< newObj.Pt() <<" "<< newObj.Eta() << " " << newObj.Phi() << std::endl;
      fjInputs_.push_back(newObj);
      idx.push_back(true);
      idx[minpair.first] = false;
      idx[minpair.second] = false;
      //std::cout << " new idx size " << idx.size()-1 << " last " << idx.back() <<" dau" << minpair.first <<" " << minpair.second << std::endl;
      history[idx.size()-1] = std::make_pair(minpair.first, minpair.second);
     }
  }

  if (mergingWTA && resumWTA)
      ResumWTA(fjInputs_, fjJets_, orgsize, jetidx, history);
  
  return true;
}       // -----  end of function P2L1TPFJetProducer::Iterations  -----


// ===  FUNCTION  ============================================================
//         Name:  P2L1TPFJetProducer::SearchHistory
//  Description:  /* cursor */
// ===========================================================================
std::vector<unsigned> P2L1TPFJetProducer::SearchHistory(const unsigned &orgsize, const unsigned searchidx, std::map<unsigned, std::pair<unsigned, unsigned> > &history) const
{
  std::vector<unsigned> returnidx;

  std::pair<unsigned, unsigned> daughter(99999, 99999);
  if (history.find(searchidx) != history.end())
    daughter = history[searchidx];

  //std::cout << searchidx <<"  " << orgsize <<" in history " << (history.find(searchidx) != history.end()) <<"  " << daughter.first <<" " << daughter.second << std::endl;
  if(daughter.first < orgsize) 
  {
    returnidx.push_back(daughter.first);
  }
  else if (daughter.first != 99999)
  {
    std::vector<unsigned> tempfirst = SearchHistory(orgsize, daughter.first, history);
    returnidx.insert(returnidx.end(), tempfirst.begin(), tempfirst.end());
  }

  if (daughter.first== daughter.second) return returnidx;

  if(daughter.second < orgsize) 
    returnidx.push_back(daughter.second);
  else if (daughter.second != 99999)
  {
    std::vector<unsigned> tempsecond = SearchHistory(orgsize, daughter.second, history);
    returnidx.insert(returnidx.end(), tempsecond.begin(), tempsecond.end());
  }

  return returnidx;
}       // -----  end of function P2L1TPFJetProducer::SearchHistory  -----

// ===  FUNCTION  ============================================================
//         Name:  P2L1TPFJetProducer::ResumWTA
//  Description:  
// ===========================================================================
void P2L1TPFJetProducer::ResumWTA( std::vector<TLorentzVector> &fjInputs_, std::vector<TLorentzVector> &fjJets_, 
  const unsigned &orgsize, std::vector<unsigned> &jetidx, std::map<unsigned, std::pair<unsigned, unsigned> > &history) const
{

  for (unsigned i = 0; i < jetidx.size(); ++i)
  {
    TLorentzVector orgjet = fjJets_.at(i);
    std::vector<unsigned> considx= SearchHistory(orgsize, 10000+jetidx[i], history);
    TLorentzVector newjet(0, 0, 0, 0);
    //std::cout << " orgpt " << orgjet.Pt() <<" eta " << orgjet.Eta() <<" phi " << orgjet.Phi() <<" mass " << orgjet.M() << " from ";
    for(auto j : considx)
    {
      //std::cout << " " << j;

      //std::cout << "jets "<< i  << ", " << j <<", " << fjInputs_.at(j).Pt() <<", " << fjInputs_.at(j).Eta() <<", " << fjInputs_.at(j).Phi() << std::endl;
      newjet += fjInputs_[j];
    }
    //std::cout << " end of input " << std::endl;
    //std::cout << " orgpt " << orgjet.Pt() <<" eta " << orgjet.Eta() <<" phi " << orgjet.Phi() <<" mass " << orgjet.M()
      //<< " newpt " << newjet.Pt() <<" eta " << newjet.Eta() <<" phi " << newjet.Phi() <<" mass " << newjet.M() << std::endl;
    fjJets_.at(i) = newjet;
  }
}       // -----  end of function P2L1TPFJetProducer::ResumWTA  -----

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
    TLorentzVector newTLV(0, 0, 0, 0);
    //winner.SetPerp(p1.Pt() + p2.Pt()); //SetPerp changed jet axis
    newTLV.SetPtEtaPhiM(p1.Pt() + p2.Pt(), winner.Eta(), winner.Phi(), winner.M());
    //std::cout << " p1 " << p1.Pt() <<" p2 " << p2.Pt() <<" newPt " << newTLV.Pt() << std::endl;
    return newTLV;
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
