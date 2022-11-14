////////////////////////////////////////////////////////////////////////
// Class:       NuSelectionAnalysis
// Plugin Type: analyzer (Unknown Unknown)
// File:        NuSelectionAnalysis_module.cc
//
// Generated at Mon Nov 14 05:34:07 2022 by Isobel Mawby using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Services
#include "art_root_io/TFileService.h"

// LArSoft
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

// Utilities
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "protoduneana/Utilities/ProtoDUNEEmptyEventFinder.h"
#include "protoduneana/Utilities/ProtoDUNEBeamCuts.h"
#include "protoduneana/Utilities/ProtoDUNEBeamlineUtils.h"

// ROOT includes
#include "TTree.h"

namespace pduneana {
  class NuSelectionAnalysis;
}


class pduneana::NuSelectionAnalysis : public art::EDAnalyzer {
public:
  explicit NuSelectionAnalysis(fhicl::ParameterSet const& p);

  NuSelectionAnalysis(NuSelectionAnalysis const&) = delete;
  NuSelectionAnalysis(NuSelectionAnalysis&&) = delete;
  NuSelectionAnalysis& operator=(NuSelectionAnalysis const&) = delete;
  NuSelectionAnalysis& operator=(NuSelectionAnalysis&&) = delete;

  // Required functions.
  void analyze(art::Event const& evt) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  void reset();

private:

  const recob::PFParticle* GetRecoBeamParticle(art::Event const& evt);
  void FillBeamQualityCutsInfo(art::Event const& evt, const recob::PFParticle* beamPFParticle);
  void FillBeamParticlePIDInfo(art::Event const& evt, const recob::PFParticle* beamPFParticle);

  TTree *fTree;

  // Run information
  int run;
  int subrun;
  int event;
  int isMC;
  int isReconstructableBeamEvent;
  bool fIsGoodBeamlineTrigger;

  // Reco beam particle info
  int fRecoBP_passBeamQualityCuts; 
  bool fRecoBP_isTrack;
  bool fRecoBP_isShower;
  double fRecoBP_TOF;
  double fRecoBP_momentumByRangeMuonHyp;
  double fRecoBP_momentumByRangeProtonHyp;
  bool fRecoBP_passMuonCuts;
  bool fRecoBP_passPionCuts;
  bool fRecoBP_passElectronCuts;
  bool fRecoBP_passProtonCuts;

  // Config
  bool fCheckSlicesForBeam;

  // From FHICL file
  std::string fPFParticleLabel;
  std::string fTrackLabel;
  std::string fShowerLabel;
  protoana::ProtoDUNEEmptyEventFinder fEmptyEventFinder;
  protoana::ProtoDUNEBeamCuts fBeamCuts;
  protoana::ProtoDUNEBeamlineUtils fBeamlineUtils;
};


pduneana::NuSelectionAnalysis::NuSelectionAnalysis(fhicl::ParameterSet const& pset)
  : EDAnalyzer{pset},
    fCheckSlicesForBeam(pset.get< bool >("CheckSlicesForBeam")),
    fPFParticleLabel(pset.get< std::string >("PFParticleLabel")),
    fTrackLabel(pset.get< std::string >("TrackLabel")),
    fShowerLabel(pset.get< std::string >("ShowerLabel")),
    fEmptyEventFinder(pset.get< fhicl::ParameterSet >("EmptyEventFinder")),
    fBeamCuts(pset.get< fhicl::ParameterSet >("BeamCuts")),
    fBeamlineUtils(pset.get<fhicl::ParameterSet>("BeamlineUtils"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void pduneana::NuSelectionAnalysis::analyze(art::Event const& evt)
{
    reset();

    run = evt.run();
    subrun = evt.subRun();
    event = evt.id().event();
    std::cout<<"########## EvtNo."<<event<<std::endl;

    // Is MC or data?
    isMC = evt.isRealData() ? 0 : 1;

    // Is this a reconstructable beam event?
    isReconstructableBeamEvent = !fEmptyEventFinder.IsEmptyEvent(evt);

    // Is this a good trigger?
    fIsGoodBeamlineTrigger = fBeamlineUtils.IsGoodBeamlineTrigger(evt);

    // Get reconstructed beam particle
    const recob::PFParticle* beamPFParticle = GetRecoBeamParticle(evt);

    if (!beamPFParticle)
    {
        std::cout << "DID NOT FIND A BEAM PARTICLE. BOOOOOOOOO." << std::endl;
        return;
    }

    // Fill beam quality cut info
    FillBeamQualityCutsInfo(evt, beamPFParticle);

    // Fill beam particle ID info
    FillBeamParticlePIDInfo(evt, beamPFParticle);

    fTree->Fill();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void pduneana::NuSelectionAnalysis::beginJob()
{
    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("tree","tree");

    fTree->Branch("run", &run);
    fTree->Branch("subrun", &subrun);
    fTree->Branch("event", &event);
    fTree->Branch("isMC", &isMC);
    fTree->Branch("isReconstructableBeamEvent", &isReconstructableBeamEvent);
    fTree->Branch("fIsGoodBeamlineTrigger", &fIsGoodBeamlineTrigger);
    fTree->Branch("checkSlicesForBeam", &fCheckSlicesForBeam);
    fTree->Branch("RecoBP_passBeamQualityCuts", &fRecoBP_passBeamQualityCuts);
    fTree->Branch("RecoBP_isTrack", &fRecoBP_isTrack);
    fTree->Branch("RecoBP_isShower", &fRecoBP_isShower);
    fTree->Branch("RecoBP_TOF", &fRecoBP_TOF);
    fTree->Branch("RecoBP_momentumByRangeMuonHyp", &fRecoBP_momentumByRangeMuonHyp);
    fTree->Branch("RecoBP_momentumByRangeProtonHyp", &fRecoBP_momentumByRangeProtonHyp);
    fTree->Branch("RecoBP_passMuonCuts", &fRecoBP_passMuonCuts);
    fTree->Branch("RecoBP_passPionCuts", &fRecoBP_passPionCuts);
    fTree->Branch("RecoBP_passElectronCuts", &fRecoBP_passElectronCuts);
    fTree->Branch("RecoBP_passProtonCuts", &fRecoBP_passProtonCuts);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void pduneana::NuSelectionAnalysis::reset()
{
  run = -1;
  subrun = -1;
  event = -1;
  isMC = -1;
  isReconstructableBeamEvent = -1;
  fRecoBP_passBeamQualityCuts = -1;
  fRecoBP_isTrack = false;
  fRecoBP_isShower = false;
  fRecoBP_TOF = -999.0;
  fRecoBP_momentumByRangeMuonHyp = -999.0;
  fRecoBP_momentumByRangeProtonHyp = -999.0;
  fRecoBP_passMuonCuts = false;
  fRecoBP_passPionCuts = false;
  fRecoBP_passElectronCuts = false;
  fRecoBP_passProtonCuts = false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void pduneana::NuSelectionAnalysis::endJob()
{
  // Implementation of optional member function here.
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This is a bit overkill. Pandora normally returns a single beam particle
// There's an odd failure mode where the beam particle can be split and the downstream half is called the beam
// This is probably what Jake is trying to deal with

const recob::PFParticle* pduneana::NuSelectionAnalysis::GetRecoBeamParticle(art::Event const& evt)
{
  protoana::ProtoDUNEPFParticleUtils pfpUtil;

  // To set in this function
  const recob::PFParticle* beamPFParticle = nullptr;

  // If you want to accept Pandora's output
  if (!fCheckSlicesForBeam)
  {
      std::vector<const recob::PFParticle*> beamParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt, fPFParticleLabel);

      if (beamParticles.size() > 0)
      {
          std::cout << "FOUND BEAM PARTICLE! WAHOO!" << std::endl;
          beamPFParticle = beamParticles.at(0);
      }

      return beamPFParticle;
  }

  // If not
  // Find the slices that contain a beam particle (normally only 1)

  int nBeamSlices = 0;
  int nBeamParticles = 0;
  std::vector<std::vector<const recob::PFParticle*>> beamSlices;

  const std::map<unsigned int, std::vector<const recob::PFParticle*>> sliceMap = pfpUtil.GetPFParticleSliceMap(evt, fPFParticleLabel);

  for (auto slice : sliceMap)
  {
    std::vector<const recob::PFParticle*> slicePFParticles = slice.second;

    for (auto slicePFParticle : slicePFParticles) 
    {
      bool added = false;
      bool isBeamParticle = pfpUtil.IsBeamParticle(*slicePFParticle, evt, fPFParticleLabel);

      if (isBeamParticle)
      {
        if (!added) 
        {
           beamSlices.push_back(slicePFParticles);
           ++nBeamSlices;
           nBeamParticles += slicePFParticles.size();
           added = true;
        }
      }
      else
      {
        continue;
      }
    }
  }

  // Now find the best beam particle
  if (nBeamSlices == 0)
  {
      std::cout << "We found no beam particles for this event... moving on" << std::endl;
      return beamPFParticle;
  }

  // Honestly I have no idea why the cut value is needed? it shifts the vertex?
  int count = 0;
  int beamCandidateSliceIndex = 0;
  double cutValue = 9999.;

  for (std::vector<const recob::PFParticle*> beamSlice : beamSlices)
  {
      // Particles ordered in beam score order
      const recob::PFParticle* beamParticleCandidate = beamSlice.at(0); 
      const TVector3 vtx = pfpUtil.GetPFParticleVertex(*beamParticleCandidate, evt, fPFParticleLabel, fTrackLabel);

      // Figure of merit (to be compared to the cut value)
      double fom = abs(vtx.Z() - 30.0); 
          
      if (fom < cutValue) 
      {
        cutValue = fom;
        beamCandidateSliceIndex = count;
      }

      ++count;
  }
  
  // Get the associated reconstructed PFParticle 

  if (beamCandidateSliceIndex >= 0)
  {
      std::cout << "FOUND BEAM PARTICLE! WAHOO!" << std::endl;
      beamSlices.at(beamCandidateSliceIndex).at(0);
  }
  else
  {
      std::vector<const recob::PFParticle*> beamParticles = pfpUtil.GetPFParticlesFromBeamSlice(evt, fPFParticleLabel);

      if (beamParticles.size() > 0)
      {
          std::cout << "FOUND BEAM PARTICLE! WAHOO!" << std::endl;
          beamPFParticle = beamParticles.at(0);
      }
  }

  return beamPFParticle;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void pduneana::NuSelectionAnalysis::FillBeamQualityCutsInfo(art::Event const& evt, const recob::PFParticle* beamPFParticle)
{
  protoana::ProtoDUNEPFParticleUtils pfpUtil;

  // Determine if the beam particle is track-like or shower-like
  const recob::Track* thisTrack = pfpUtil.GetPFParticleTrack(*beamPFParticle, evt, fPFParticleLabel, fTrackLabel);
  const recob::Shower* thisShower = pfpUtil.GetPFParticleShower(*beamPFParticle, evt, fPFParticleLabel, fShowerLabel);

  if (thisTrack)
  {
    fRecoBP_passBeamQualityCuts = (fBeamCuts.IsBeamlike(*thisTrack, evt, "1") ? 1 : 0);
  }
  else if (thisShower)
  {
    fRecoBP_passBeamQualityCuts = (fBeamCuts.IsBeamlike(*thisShower, evt, "1") ? 1 : 0);
  }

  if (fRecoBP_passBeamQualityCuts  == 1)
  {
      std::cout << "PASSES THE BEAM CUTS" << std::endl;
  }
  else
  {
      std::cout << "DOES NOT PASS THE BEAM CUTS" << std::endl;
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void pduneana::NuSelectionAnalysis::FillBeamParticlePIDInfo(art::Event const& evt, const recob::PFParticle* beamPFParticle)
{
  protoana::ProtoDUNEPFParticleUtils pfpUtil;

  // Determine if the beam particle is track-like or shower-like
  const recob::Track* thisTrack = pfpUtil.GetPFParticleTrack(*beamPFParticle, evt, fPFParticleLabel, fTrackLabel);
  const recob::Shower* thisShower = pfpUtil.GetPFParticleShower(*beamPFParticle, evt, fPFParticleLabel, fShowerLabel);

  if (thisTrack)
    fRecoBP_isTrack = true;
  else if (thisShower)
    fRecoBP_isShower = true;

  // 1 here is for 1GeV
  const auto possibleParticles = fBeamlineUtils.GetPIDCandidates(evt, 1.0);

  fRecoBP_passMuonCuts = possibleParticles.muon;
  fRecoBP_passPionCuts = possibleParticles.pion;
  fRecoBP_passElectronCuts = possibleParticles.electron;
  fRecoBP_passProtonCuts = possibleParticles.proton;

  if (fRecoBP_isTrack)
  {
    // Calculate momentum according to CSDA range

    trkf::TrackMomentumCalculator trackMomCalc;

    fRecoBP_momentumByRangeMuonHyp = trackMomCalc.GetTrackMomentum(thisTrack->Length(), 13);
    fRecoBP_momentumByRangeProtonHyp = trackMomCalc.GetTrackMomentum(thisTrack->Length(), 2212);
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


DEFINE_ART_MODULE(pduneana::NuSelectionAnalysis)
