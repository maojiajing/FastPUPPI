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

#ifndef  MY_P2L1TPFJETPRODUCER_INC
#define  MY_P2L1TPFJETPRODUCER_INC


// system include files
#include <memory>
#include <cmath>
#include <cstdint>

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
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "RecoJets/JetProducers/interface/JetSpecific.h"

#include "TLorentzVector.h"
//
// class declaration
//

typedef uint16_t fixed_point_t;
#define FIXED_POINT_FRACTIONAL_BITS 6

inline double fixed_to_float(fixed_point_t input)
{
  return ((double)input / (double)(1 << FIXED_POINT_FRACTIONAL_BITS));
}

inline fixed_point_t float_to_fixed(double input)
{
  return (fixed_point_t)(input * (1 << FIXED_POINT_FRACTIONAL_BITS));
}


inline double rePrecision(double input)
{
  return fixed_to_float(float_to_fixed(input));
}

class P2L1TPFJetProducer : public edm::stream::EDProducer<> {
   public:
      explicit P2L1TPFJetProducer(const edm::ParameterSet&);
      ~P2L1TPFJetProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      edm::InputTag         src_;                       // input constituent source

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      inline float Metric_dR(TLorentzVector p1, TLorentzVector p2, bool sameobj) const;
      inline float Metric_AK(TLorentzVector p1, TLorentzVector p2, bool sameobj) const;
      inline TLorentzVector MergingE(TLorentzVector &p1, TLorentzVector &p2) const;
      inline TLorentzVector MergingWTA(TLorentzVector &p1, TLorentzVector &p2) const;
      bool Iterations( std::vector<TLorentzVector> &fjInputs_, std::vector<TLorentzVector> &fjJets_);
      bool GetCombinationsKt( std::vector<TLorentzVector> &fjInputs_, std::map<float, std::pair<unsigned, unsigned> > &Valcomb);
      bool GetCombinationsdR( std::vector<TLorentzVector> &fjInputs_, std::map<float, std::pair<unsigned, unsigned> > &Valcomb);
      void GetInput(edm::Event& iEvent, const edm::EventSetup& iSetup);
      void WriteJetCollection(const edm::EventSetup& iSetup);
      bool RemergeJets();
      void SplitN2Groups();
      void SplitN2Groupsv2();
      void SplitN2Tile(int Neta, int Nphi);
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<l1t::PFCandidateCollection> ParTok_;
      edm::Handle<l1t::PFCandidateCollection> ParHdl;
      std::unique_ptr<reco::PFJetCollection> jetCollection;

      std::vector<TLorentzVector> fjInputs;        // fastjet inputs
      std::vector<TLorentzVector> fjJets;          // fastjet jets

      std::vector<std::vector<TLorentzVector> > gInputs;
      std::vector<std::vector<TLorentzVector> > gJets;


      //std::map<std::pair<unsigned, unsigned>, float > combVal;
      //std::map<float, std::pair<unsigned, unsigned> > Valcomb;


      double rParam;
      double jetPtMin;
      double inputEtMin;
      int ktsign;
      bool metricKt;
      bool metricdR;
      bool mergingE;
      bool mergingWTA;
      bool N2Tile;
      bool N2Group;

};

#endif   // ----- #ifndef MY_P2L1TPFJETPRODUCER_INC  -----

