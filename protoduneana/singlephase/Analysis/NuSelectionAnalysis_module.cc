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
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Services
#include "art_root_io/TFileService.h"

// LArSoft
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/AnalysisBase/MVAPIDResult.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "dunecore/DuneObj/ProtoDUNEBeamEvent.h"

// Utilities
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "protoduneana/Utilities/ProtoDUNEEmptyEventFinder.h"
#include "protoduneana/Utilities/ProtoDUNEBeamCuts.h"
#include "protoduneana/Utilities/ProtoDUNEBeamlineUtils.h"
#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/Utilities/ProtoDUNEShowerUtils.h"
#include "protoduneana/Utilities/ProtoDUNECalibration.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"

// ROOT includes
#include "TTree.h"
#include "TGraph.h"
#include "TFile.h"
#include "TSpline.h"

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

  void FillTrueBeamParticleInfo(art::Event const& evt);
  const recob::PFParticle* GetRecoBeamParticle(art::Event const& evt);
  void FillBeamQualityCutsInfo(art::Event const& evt, const recob::PFParticle* beamPFParticle);
  void FillBeamParticlePIDInfo(art::Event const& evt, const recob::PFParticle* beamPFParticle);
  void FillRecoBeamInfo(art::Event const& evt);
  void FillMuonBDTInfo(art::Event const& evt, const recob::PFParticle* beamPFParticle);
  void FillWarwickMVAInfo(art::Event const& evt, const art::Ptr<recob::Track> thisTrack);
  void FilldEdxInfo(art::Event const& evt, const recob::Track* thisTrack);
  void FillTrackLengthInfo(art::Event const& evt, const recob::Track* thisTrack);
  void FillTrackDeflectionInfo(art::Event const& evt, const recob::Track* thisTrack);
  void FillMichelInfo(art::Event const& evt, const recob::PFParticle* beamPFParticle, const recob::Track* thisTrack);
  void FillCandidateMichelEnergyInfo(art::Event const& evt, const recob::PFParticle* michelCandidate);
  void FillCSDAScore(art::Event const& evt);

  bool IsMuonBDTSignal();
  bool IsMuonBDTBackground();
  TVector3 ApplyPositionSCECorrection(const geo::Point_t &inputPosition);

  TTree *fTree;
  TTree *fMuonBDTSignalTree;
  TTree *fMuonBDTBackgroundTree;

  TGraph * RvsKE;
  TSpline3 * RvsKESpline;

  // Run information
  int fRun;
  int fSubrun;
  int fEvent;
  int fIsMC;
  int fIsReconstructableBeamEvent;
  bool fIsGoodBeamlineTrigger;

  // Reco beam info
  double fRecoBeamMomentum;

  // True beam particle info
  int fTrueBP_PDG;
  double fTrueBeamMomentum;

  // Reco beam particle info
  int fRecoBP_passBeamQualityCuts; 
  bool fRecoBP_isTrack;
  bool fRecoBP_isShower;
  double fRecoBP_momentumByRangeMuonHyp;
  double fRecoBP_momentumByRangeProtonHyp;
  bool fRecoBP_passMuonCuts;
  bool fRecoBP_passPionCuts;
  bool fRecoBP_passElectronCuts;
  bool fRecoBP_passProtonCuts;
  double fRecoBP_CSDAScore;

  // Muon BDT variables (can't do a SCE version)
  double fRecoBP_EvalRatio_MuonBDT;
  double fRecoBP_Concentration_MuonBDT;
  double fRecoBP_CoreHaloRatio_MuonBDT;
  double fRecoBP_Conicalness_MuonBDT;

  double fRecoBP_dEdxStart_MuonBDT;
  double fRecoBP_dEdxEnd_MuonBDT;
  double fRecoBP_dEdxEndRatio_MuonBDT;
  double fRecoBP_TrackLength_MuonBDT;
  double fRecoBP_DeflecAngleSD_MuonBDT;
  double fRecoBP_candidateMichelNHits_MuonBDT;
  double fRecoBP_candidateMichelEnergy_MuonBDT;

  // PID
  double fRecoBP_protonChi2;
  double fRecoBP_muonChi2;

  // Config
  bool fRecalibrate;
  bool fCheckSlicesForBeam;
  protoana::ProtoDUNECalibration fCalibrationSCE;

  // From FHICL file
  std::string fGeneratorLabel;
  std::string fPFParticleLabel;
  std::string fHitLabel;
  std::string fTrackLabel;
  std::string fShowerLabel;
  std::string fPIDLabel; 
  std::string fSCECaloLabel;
  std::string fBeamLabel;
  protoana::ProtoDUNEEmptyEventFinder fEmptyEventFinder;
  protoana::ProtoDUNEBeamCuts fBeamCuts;
  protoana::ProtoDUNEBeamlineUtils fBeamlineUtils;
  protoana::ProtoDUNETrackUtils fTrackUtil;
  protoana::ProtoDUNEShowerUtils fShowerUtil;
  protoana::ProtoDUNETruthUtils fTruthUtil;
  protoana::ProtoDUNEPFParticleUtils fPFPUtil;
  art::ServiceHandle<geo::Geometry> fGeometryService;

  // Used in the Chi2 PID
  std::vector<double> fCalibrated_dEdX;
  std::vector<double> fResidualRange;
};


pduneana::NuSelectionAnalysis::NuSelectionAnalysis(fhicl::ParameterSet const& pset)
  : EDAnalyzer{pset},
    fRecalibrate(pset.get< bool >("Recalibrate")),
    fCheckSlicesForBeam(pset.get< bool >("CheckSlicesForBeam")),
    fCalibrationSCE(pset.get< fhicl::ParameterSet >("CalibrationParsSCE")),
    fGeneratorLabel(pset.get< std::string >("GeneratorLabel")),
    fPFParticleLabel(pset.get< std::string >("PFParticleLabel")),
    fHitLabel(pset.get< std::string >("HitLabel")),
    fTrackLabel(pset.get< std::string >("TrackLabel")),
    fShowerLabel(pset.get< std::string >("ShowerLabel")),
    fPIDLabel(pset.get< std::string >("PIDLabel")),
    fSCECaloLabel(pset.get< std::string >("SCECaloLabel")),
    fBeamLabel(pset.get< std::string >("BeamLabel")),
    fEmptyEventFinder(pset.get< fhicl::ParameterSet >("EmptyEventFinder")),
    fBeamCuts(pset.get< fhicl::ParameterSet >("BeamCuts")),
    fBeamlineUtils(pset.get<fhicl::ParameterSet>("BeamlineUtils"))
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void pduneana::NuSelectionAnalysis::analyze(art::Event const& evt)
{
    reset();

    fRun = evt.run();
    fSubrun = evt.subRun();
    fEvent = evt.id().event();
    std::cout<<"########## EvtNo."<< fEvent << std::endl;

    // Is MC or data?
    fIsMC = evt.isRealData() ? 0 : 1;

    // Is this a reconstructable beam event?
    fIsReconstructableBeamEvent = !fEmptyEventFinder.IsEmptyEvent(evt);

    // Is this a good trigger?
    fIsGoodBeamlineTrigger = fBeamlineUtils.IsGoodBeamlineTrigger(evt);

    FillRecoBeamInfo(evt);
    FillTrueBeamParticleInfo(evt);

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

    // Is track? Fill Pandizzle variables
    if (fRecoBP_isTrack)
    {
        FillMuonBDTInfo(evt, beamPFParticle);
        FillCSDAScore(evt);
    }
    // Is shower? Fill Pandrizzle variables

    // Fill trees
    fTree->Fill();

    if (IsMuonBDTSignal())
        fMuonBDTSignalTree->Fill();

    if (IsMuonBDTBackground())
        fMuonBDTBackgroundTree->Fill();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void pduneana::NuSelectionAnalysis::beginJob()
{
    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("tree","tree");
    fMuonBDTSignalTree = tfs->make<TTree>("MuonBDTSignalTree", "MuonBDTSignalTree");
    fMuonBDTBackgroundTree = tfs->make<TTree>("MuonBDTBackgroundTree", "MuonBDTBackgroundTree");

    double range_cm[29] =
            {9.833E-1, 1.786E0, 3.321E0, 6.598E0, 1.058E1, 3.084E1, 4.250E1, 6.732E1, 1.063E2, 1.725E2,
             2.385E2,  4.934E2, 6.163E2, 8.552E2, 1.202E3, 1.758E3, 2.297E3, 4.359E3, 5.354E3, 7.298E3,
             1.013E4,  1.469E4, 1.910E4, 3.558E4, 4.326E4, 5.768E4, 7.734E4, 1.060E5, 1.307E5};
    
    for (double& value : range_cm)
      value /= 1.396; // convert to cm

    double KE_MeV[29] = 
      {10, 14, 20, 30, 40, 80, 100, 140, 200, 300, 400, 800, 1000, 1400, 2000, 3000, 4000, 8000, 10000, 14000,
       20000, 30000, 40000, 80000, 100000, 140000, 200000, 300000, 400000};

    RvsKE = tfs->make<TGraph>(29, KE_MeV, range_cm);
    RvsKESpline = tfs->make<TSpline3>("RvsKES", RvsKE);

    fTree->Branch("Run", &fRun);
    fTree->Branch("Subrun", &fSubrun);
    fTree->Branch("Event", &fEvent);
    fTree->Branch("IsMC", &fIsMC);
    fTree->Branch("IsReconstructableBeamEvent", &fIsReconstructableBeamEvent);
    fTree->Branch("IsGoodBeamlineTrigger", &fIsGoodBeamlineTrigger);
    fTree->Branch("CheckSlicesForBeam", &fCheckSlicesForBeam);
    fTree->Branch("RecoBeamMomentum", &fRecoBeamMomentum);
    fTree->Branch("TrueBeamMomentum", &fTrueBeamMomentum);
    fTree->Branch("TrueBP_PDG", &fTrueBP_PDG);
    fTree->Branch("RecoBP_passBeamQualityCuts", &fRecoBP_passBeamQualityCuts);
    fTree->Branch("RecoBP_isTrack", &fRecoBP_isTrack);
    fTree->Branch("RecoBP_isShower", &fRecoBP_isShower);
    fTree->Branch("RecoBP_momentumByRangeMuonHyp", &fRecoBP_momentumByRangeMuonHyp);
    fTree->Branch("RecoBP_momentumByRangeProtonHyp", &fRecoBP_momentumByRangeProtonHyp);
    fTree->Branch("RecoBP_passMuonCuts", &fRecoBP_passMuonCuts);
    fTree->Branch("RecoBP_passPionCuts", &fRecoBP_passPionCuts);
    fTree->Branch("RecoBP_passElectronCuts", &fRecoBP_passElectronCuts);
    fTree->Branch("RecoBP_passProtonCuts", &fRecoBP_passProtonCuts);
    fTree->Branch("RecoBP_CSDAScore", &fRecoBP_CSDAScore);
    fTree->Branch("RecoBP_protonChi2", &fRecoBP_protonChi2);
    fTree->Branch("RecoBP_muonChi2", &fRecoBP_muonChi2);

    for (TTree *tree : {fTree, fMuonBDTSignalTree, fMuonBDTBackgroundTree})
    {
      tree->Branch("RecoBP_EvalRatio_MuonBDT", &fRecoBP_EvalRatio_MuonBDT);
      tree->Branch("RecoBP_Concentration_MuonBDT", &fRecoBP_Concentration_MuonBDT);
      tree->Branch("RecoBP_CoreHaloRatio_MuonBDT", &fRecoBP_CoreHaloRatio_MuonBDT);
      tree->Branch("RecoBP_Conicalness_MuonBDT", &fRecoBP_Conicalness_MuonBDT);
      tree->Branch("RecoBP_dEdxStart_MuonBDT", &fRecoBP_dEdxStart_MuonBDT);
      tree->Branch("RecoBP_dEdxEnd_MuonBDT", &fRecoBP_dEdxEnd_MuonBDT);
      tree->Branch("RecoBP_dEdxEndRatio_MuonBDT", &fRecoBP_dEdxEndRatio_MuonBDT);
      tree->Branch("RecoBP_TrackLength_MuonBDT", &fRecoBP_TrackLength_MuonBDT);
      tree->Branch("RecoBP_DeflecAngleSD_MuonBDT", &fRecoBP_DeflecAngleSD_MuonBDT);
      tree->Branch("RecoBP_candidateMichelNHits_MuonBDT", &fRecoBP_candidateMichelNHits_MuonBDT);
      tree->Branch("RecoBP_candidateMichelEnergy_MuonBDT", &fRecoBP_candidateMichelEnergy_MuonBDT);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void pduneana::NuSelectionAnalysis::reset()
{
  fRun = -1;
  fSubrun = -1;
  fEvent = -1;
  fIsMC = -1;
  fIsReconstructableBeamEvent = -1;
  fRecoBeamMomentum = -999.0;
  fTrueBP_PDG = -999;
  fTrueBeamMomentum = -999.0;
  fRecoBP_passBeamQualityCuts = -1;
  fRecoBP_isTrack = false;
  fRecoBP_isShower = false;
  fRecoBP_momentumByRangeMuonHyp = -999.0;
  fRecoBP_momentumByRangeProtonHyp = -999.0;
  fRecoBP_passMuonCuts = false;
  fRecoBP_passPionCuts = false;
  fRecoBP_passElectronCuts = false;
  fRecoBP_passProtonCuts = false;
  fRecoBP_CSDAScore = -999.0;
  fRecoBP_EvalRatio_MuonBDT = -999.0;
  fRecoBP_Concentration_MuonBDT = -999.0;
  fRecoBP_CoreHaloRatio_MuonBDT = -999.0;
  fRecoBP_Conicalness_MuonBDT = -999.0;
  fRecoBP_dEdxStart_MuonBDT = -999.0;
  fRecoBP_dEdxEnd_MuonBDT = -999.0;
  fRecoBP_dEdxEndRatio_MuonBDT = -999.0;
  fRecoBP_TrackLength_MuonBDT = -999.0;
  fRecoBP_DeflecAngleSD_MuonBDT = -999.0;
  fRecoBP_candidateMichelNHits_MuonBDT = -999.0;
  fRecoBP_candidateMichelEnergy_MuonBDT = -999.0;

  fRecoBP_protonChi2 = -999.0;
  fRecoBP_muonChi2 = -999.0;

  fCalibrated_dEdX.clear();
  fResidualRange.clear();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void pduneana::NuSelectionAnalysis::endJob()
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void pduneana::NuSelectionAnalysis::FillRecoBeamInfo(art::Event const& evt)
{
  // Get beam momentum
  std::vector<art::Ptr<beam::ProtoDUNEBeamEvent>> beamEvents;
  auto beamHandle = evt.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>(fBeamLabel);
  art::fill_ptr_vector(beamEvents, beamHandle);

  if (beamEvents.empty())
    return;

  art::Ptr<beam::ProtoDUNEBeamEvent> beamEvent = beamEvents.at(0);

  std::vector<double> momenta = beamEvent->GetRecoBeamMomenta(); 

  if (momenta.empty())
    return;

  fRecoBeamMomentum = momenta.at(0);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void pduneana::NuSelectionAnalysis::FillTrueBeamParticleInfo(art::Event const& evt)
{
  if (evt.isRealData())
    return;

  // This gets the true beam particle that generated the event
  const simb::MCParticle* trueBeamParticle = nullptr;

  auto mcTruths = evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorLabel);

  trueBeamParticle = fTruthUtil.GetGeantGoodParticle((*mcTruths)[0], evt);

  if (!trueBeamParticle)
    return;

  fTrueBP_PDG = trueBeamParticle->PdgCode();

  fTrueBeamMomentum = trueBeamParticle->P();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This is a bit overkill. Pandora normally returns a single beam particle
// There's an odd failure mode where the beam particle can be split and the downstream half is called the beam
// This is probably what Jake is trying to deal with

const recob::PFParticle* pduneana::NuSelectionAnalysis::GetRecoBeamParticle(art::Event const& evt)
{
  // To set in this function
  const recob::PFParticle* beamPFParticle = nullptr;

  // If you want to accept Pandora's output
  if (!fCheckSlicesForBeam)
  {
      std::vector<const recob::PFParticle*> beamParticles = fPFPUtil.GetPFParticlesFromBeamSlice(evt, fPFParticleLabel);

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

  const std::map<unsigned int, std::vector<const recob::PFParticle*>> sliceMap = fPFPUtil.GetPFParticleSliceMap(evt, fPFParticleLabel);

  for (auto slice : sliceMap)
  {
    std::vector<const recob::PFParticle*> slicePFParticles = slice.second;

    for (auto slicePFParticle : slicePFParticles) 
    {
      bool added = false;
      bool isBeamParticle = fPFPUtil.IsBeamParticle(*slicePFParticle, evt, fPFParticleLabel);

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
      const TVector3 vtx = fPFPUtil.GetPFParticleVertex(*beamParticleCandidate, evt, fPFParticleLabel, fTrackLabel);

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
      std::vector<const recob::PFParticle*> beamParticles = fPFPUtil.GetPFParticlesFromBeamSlice(evt, fPFParticleLabel);

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
  // Determine if the beam particle is track-like or shower-like
  const recob::Track* thisTrack = fPFPUtil.GetPFParticleTrack(*beamPFParticle, evt, fPFParticleLabel, fTrackLabel);
  const recob::Shower* thisShower = fPFPUtil.GetPFParticleShower(*beamPFParticle, evt, fPFParticleLabel, fShowerLabel);

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
  // Determine if the beam particle is track-like or shower-like
  const recob::Track* thisTrack = fPFPUtil.GetPFParticleTrack(*beamPFParticle, evt, fPFParticleLabel, fTrackLabel);
  const recob::Shower* thisShower = fPFPUtil.GetPFParticleShower(*beamPFParticle, evt, fPFParticleLabel, fShowerLabel);

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
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void pduneana::NuSelectionAnalysis::FillMuonBDTInfo(art::Event const& evt, const recob::PFParticle* beamPFParticle)
{
  const recob::Track* thisTrack = fPFPUtil.GetPFParticleTrack(*beamPFParticle, evt, fPFParticleLabel, fTrackLabel);

  if (!thisTrack)
      return;

  FilldEdxInfo(evt, thisTrack);
  FillTrackLengthInfo(evt, thisTrack);
  FillTrackDeflectionInfo(evt, thisTrack);
  FillMichelInfo(evt, beamPFParticle, thisTrack);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void pduneana::NuSelectionAnalysis::FillWarwickMVAInfo(art::Event const& evt, const art::Ptr<recob::Track> thisTrack)
{
  // Get the MVA result for the track
  art::ValidHandle< std::vector<recob::Track> > trackListHandle = evt.getValidHandle<std::vector<recob::Track>>(fTrackLabel);
  art::FindManyP<anab::MVAPIDResult> findMVA(trackListHandle, evt, fPIDLabel);
  std::vector<art::Ptr<anab::MVAPIDResult> > pids = findMVA.at(thisTrack.key());

  // Check that MVA exists
  if (pids.size() == 0)
      return;

  art::Ptr<anab::MVAPIDResult> pid = pids.at(0);

  fRecoBP_EvalRatio_MuonBDT = pid->evalRatio;
  fRecoBP_Concentration_MuonBDT = pid->concentration;
  fRecoBP_CoreHaloRatio_MuonBDT = pid->coreHaloRatio;
  fRecoBP_Conicalness_MuonBDT = pid->conicalness;
  fRecoBP_dEdxStart_MuonBDT = pid->dEdxStart;
  fRecoBP_dEdxEnd_MuonBDT = pid->dEdxEnd;
  fRecoBP_dEdxEndRatio_MuonBDT = pid->dEdxEndRatio;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void pduneana::NuSelectionAnalysis::FilldEdxInfo(art::Event const& evt, const recob::Track* thisTrack)
{
  // Primary Track Calorimetry
  // SCE-corrected
  auto trackCalo = fTrackUtil.GetRecoTrackCalorimetry(*thisTrack, evt, fTrackLabel, fSCECaloLabel);

  // Get the correct plane (collection) because it's not always the same
  bool foundCorrectPlane = false;
  size_t index = 0;

  for (size_t i = 0; i < trackCalo.size(); ++i) 
  {
      if (trackCalo[i].PlaneID().Plane == 2) 
      {
          foundCorrectPlane = true;
          index = i;
          break; 
      }
  }

  if (foundCorrectPlane) 
  {
      auto calo_dEdX = trackCalo[index].dEdx();
      auto calo_range = trackCalo[index].ResidualRange();

      if (fRecalibrate)
      {
        std::vector<float> new_dEdX = fCalibrationSCE.GetCalibratedCalorimetry(*thisTrack, evt, fTrackLabel, fSCECaloLabel, 2, -10.);

        for( size_t i = 0; i < new_dEdX.size(); ++i )
        {
          fCalibrated_dEdX.push_back(new_dEdX[i]);
          fResidualRange.push_back(calo_range[i]);
        }
      }
      else
      {
        for (size_t i = 0; i < calo_dEdX.size(); ++i)
        {
          fCalibrated_dEdX.push_back(calo_dEdX[i]);
          fResidualRange.push_back(calo_range[i]);
        }
      }
  }

  // First entry is at the end of the track (LArSoft hates you) but don't use end points
  if (fCalibrated_dEdX.size() > 2)
  {
    fRecoBP_dEdxStart_MuonBDT = fCalibrated_dEdX.at(fCalibrated_dEdX.size() - 2);
    fRecoBP_dEdxEnd_MuonBDT = fCalibrated_dEdX.at(1);
  }

  if (fCalibrated_dEdX.size() > 3)
    fRecoBP_dEdxEndRatio_MuonBDT = fCalibrated_dEdX.at(1) / fCalibrated_dEdX.at(2);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void pduneana::NuSelectionAnalysis::FillTrackLengthInfo(art::Event const& evt, const recob::Track* thisTrack)
{
  // To correct for the SCE
  auto sce = lar::providerFrom<spacecharge::SpaceChargeService>();

  geo::Point_t startPoint = thisTrack->Start();
  geo::Vector_t startOffset = sce->GetPosOffsets(startPoint);
  TVector3 startPointSCE = TVector3(startPoint.X() - startOffset.X(), startPoint.Y() + startOffset.Y(), startPoint.Z() + startOffset.Z());

  geo::Point_t endPoint = thisTrack->End();
  geo::Vector_t endOffset = sce->GetPosOffsets(endPoint);
  TVector3 endPointSCE = TVector3(endPoint.X() - endOffset.X(), endPoint.Y() + endOffset.Y(), endPoint.Z() + endOffset.Z());

  TVector3 lengthVectorSCE = endPointSCE - startPointSCE;
  double lengthSCE = sqrt((lengthVectorSCE.X() * lengthVectorSCE.X()) + (lengthVectorSCE.Y() * lengthVectorSCE.Y()) + (lengthVectorSCE.Z() * lengthVectorSCE.Z()));

  std::cout << "HERE HERE HERE" << std::endl;


  fRecoBP_TrackLength_MuonBDT = lengthSCE;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Totally copied from Dom, sass retained
void pduneana::NuSelectionAnalysis::FillTrackDeflectionInfo(art::Event const& evt, const recob::Track* thisTrack)
{
  // To correct for the SCE
  auto sce = lar::providerFrom<spacecharge::SpaceChargeService>();

  // Loop over the trajectory points in a track and compare adjacent pairs
  // Get the number of points
  size_t nPoints = thisTrack->NumberTrajectoryPoints();

  // Store the directions between adjacent points on a vector
  std::vector<TVector3> directions;
  for (size_t i_point = 0; i_point < (nPoints - 1); i_point++)
  {
    TVector3 position_i(thisTrack->TrajectoryPoint(i_point).position.X(), thisTrack->TrajectoryPoint(i_point).position.Y(), thisTrack->TrajectoryPoint(i_point).position.Z());
    auto offset_i = sce->GetPosOffsets({position_i.X(), position_i.Y(), position_i.Z()});
    TVector3 positionSCE_i(position_i.X() - offset_i.X(), position_i.Y() + offset_i.Y(), position_i.Z() + offset_i.Z());

    TVector3 position_iplus1(thisTrack->TrajectoryPoint(i_point + 1).position.X(), thisTrack->TrajectoryPoint(i_point + 1).position.Y(), thisTrack->TrajectoryPoint(i_point + 1).position.Z());
    auto offset_iplus1 = sce->GetPosOffsets({position_iplus1.X(), position_iplus1.Y(), position_iplus1.Z()});
    TVector3 positionSCE_iplus1(position_iplus1.X() - offset_iplus1.X(), position_iplus1.Y() + offset_iplus1.Y(), position_iplus1.Z() + offset_iplus1.Z());

    TVector3 direction = (positionSCE_iplus1 - positionSCE_i).Unit();
    directions.push_back(direction);
  }

  // Now let's loop through the direction and compare adjacent elements (wow!)
  std::vector<double> deflection_angles;
  for (size_t i_dir = 0; i_dir < (directions.size() - 1); i_dir++)
  {
    // Aim: rotate both direction so that the first direction is parallel to the z-axis.  Then take the x-projection of scattered track and calculate the angle between that and the first direction
    TVector3 z_axis(0, 0, 1);
    TVector3 direction_first = directions[i_dir];
    TVector3 direction_second = directions[i_dir + 1];

    // Ignore if either direction is 0 (i.e. not 1)
    if (direction_first.Mag() < 0.999 || direction_second.Mag() < 0.999)
          continue;
      
    double angle_dir_first_z_axis = direction_first.Angle(z_axis);
    TVector3 orthogonal_vector = direction_first.Cross(z_axis);

    if (orthogonal_vector.Unit().Mag() < 0.999)
        continue;

    direction_first.Rotate(angle_dir_first_z_axis, orthogonal_vector);
    direction_second.Rotate(angle_dir_first_z_axis, orthogonal_vector);
   
    // Now work out the angle between the vectors in th x-z plane
    direction_first.SetY(0);
    direction_second.SetY(0);

    double dot_product = direction_first.Dot(direction_second);
    dot_product = std::min(std::max(dot_product,-1.),1.);

    double angle = acos(dot_product) * 180 / 3.142;

    if (direction_second.X() < 0) 
      angle *= -1;

    deflection_angles.push_back(angle);
  }

  double angle_mean = 0;

  for (size_t i_angle = 0; i_angle < deflection_angles.size(); i_angle++)
      angle_mean += deflection_angles[i_angle];
  
  if (deflection_angles.size() > 0) 
    angle_mean /= deflection_angles.size();
  else 
    angle_mean = -100;

  double angle_var = 0;

  for (size_t i_angle = 0; i_angle < deflection_angles.size(); i_angle++)
      angle_var = (deflection_angles[i_angle] - angle_mean) * (deflection_angles[i_angle] - angle_mean);
  
  if (deflection_angles.size() > 1) 
    angle_var /= (deflection_angles.size() - 1);
  else 
    angle_var = -2.;

  fRecoBP_DeflecAngleSD_MuonBDT = (angle_var > 0.0) ? sqrt(angle_var) : -2.f;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void pduneana::NuSelectionAnalysis::FillMichelInfo(art::Event const& evt, const recob::PFParticle* beamPFParticle, const recob::Track* thisTrack)
{
  // Need to find a Michel candidate
  auto eventPFParticles = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);

  // Why is particle physics so gendered -.-
  std::vector<const recob::PFParticle*> childrenPFParticles;
  for (size_t childID : beamPFParticle->Daughters())
    childrenPFParticles.push_back(&eventPFParticles->at(childID));

  if (childrenPFParticles.empty())
    return;

  // Need to find SCE corrected track start/end point
  TVector3 trackEndpointSCE = ApplyPositionSCECorrection(thisTrack->End());

  const recob::PFParticle* michelCandidate = nullptr;
  double candidateSeparation = std::numeric_limits<double>::max();

  for (const recob::PFParticle* childPFParticle : childrenPFParticles)
  {
    double separation = std::numeric_limits<double>::max();

    const recob::Track* childTrack = fPFPUtil.GetPFParticleTrack(*childPFParticle, evt, fPFParticleLabel, fTrackLabel);

    if (childTrack)
    {
      TVector3 childStartpointSCE = ApplyPositionSCECorrection(childTrack->Start());
      separation = (childStartpointSCE - trackEndpointSCE).Mag();
    }

    const recob::Shower* childShower = fPFPUtil.GetPFParticleShower(*childPFParticle, evt, fPFParticleLabel, fShowerLabel);

    if (childShower)
    {
      geo::Point_t showerStart(childShower->ShowerStart().X(), childShower->ShowerStart().Y(), childShower->ShowerStart().Z());
      TVector3 childStartpointSCE = ApplyPositionSCECorrection(showerStart);
      separation = (childStartpointSCE - trackEndpointSCE).Mag();
    }

    if ((separation < 2.f) && (separation < candidateSeparation))
    {
      candidateSeparation = separation;
      michelCandidate = childPFParticle;
    }
  }

  if (!michelCandidate)
    return;

  // Fill hit info
  const std::vector<const recob::Hit*> michelHitVector = fPFPUtil.GetPFParticleHits(*michelCandidate, evt, fPFParticleLabel);
  fRecoBP_candidateMichelNHits_MuonBDT = michelHitVector.size();

  //Get shower energy
  FillCandidateMichelEnergyInfo(evt, michelCandidate);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void pduneana::NuSelectionAnalysis::FillCandidateMichelEnergyInfo(art::Event const& evt, const recob::PFParticle* michelCandidate)
{
  const std::vector< art::Ptr< recob::Hit > > candidateMichelHits = fPFPUtil.GetPFParticleHits_Ptrs(*michelCandidate, evt, fPFParticleLabel);

  // This is completely stolen from Jake
  // For hit energy calculation we need position of each hit
  // For some reason we have to find an average y, because sometimes y is non existant?
  art::FindManyP<recob::SpacePoint> spFromHits(candidateMichelHits, evt, fHitLabel);

  std::vector<double> xPositions, yPositions, zPositions;
  double totalY = 0.;
  int nGoodY = 0;

  std::vector<art::Ptr<recob::Hit>> collectionHits;

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
  auto const detProp =  art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(evt, clockData);

  for (size_t iHit = 0; iHit < candidateMichelHits.size(); ++iHit) 
  {
      auto theHit = candidateMichelHits[iHit];

      if (theHit->View() != 2) 
        continue;

      collectionHits.push_back(theHit);

      double xPosition = detProp.ConvertTicksToX(theHit->PeakTime(), theHit->WireID().Plane, theHit->WireID().TPC, 0);
      double zPosition = fGeometryService->Wire(theHit->WireID()).GetCenter().Z();

      xPositions.push_back(xPosition);
      zPositions.push_back(zPosition);

      std::vector<art::Ptr<recob::SpacePoint>> sps = spFromHits.at(iHit);

      if (!sps.empty()) 
      {
          yPositions.push_back(sps[0]->XYZ()[1]);
          totalY += yPositions.back();
          ++nGoodY;
      }
      else 
      {
          yPositions.push_back(-999.);
      }
  }

  if (nGoodY == 0)
    return;

  double totalEnergy = 0.;

  for (size_t iHit = 0; iHit < collectionHits.size(); ++iHit) 
  {
      auto theHit = collectionHits[iHit];

      if (theHit->View() != 2) 
        continue;

      if (yPositions[iHit] < -100.)
          yPositions[iHit] = totalY / nGoodY;

      totalEnergy += fCalibrationSCE.HitToEnergy(collectionHits[iHit], xPositions[iHit], yPositions[iHit], zPositions[iHit]);
  }

  fRecoBP_candidateMichelEnergy_MuonBDT = totalEnergy;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void pduneana::NuSelectionAnalysis::FillCSDAScore(art::Event const& evt)
{
  // Convert to KE
  double muonMass = 105.7;
  double beamMomentum =  evt.isRealData() ? fRecoBeamMomentum : fTrueBeamMomentum;
  double kineticEnergy = sqrt((1.e6 * beamMomentum * beamMomentum) + (muonMass * muonMass)) - muonMass;
  double CSDARange = RvsKESpline->Eval(kineticEnergy);


  std::cout << "beamMomentum: " << beamMomentum << std::endl;
  std::cout << "kineticEnergy: " << kineticEnergy << std::endl;
  std::cout << "CSDARange: " << CSDARange << std::endl;
  std::cout << "fRecoBP_TrackLength_MuonBDT: " << fRecoBP_TrackLength_MuonBDT << std::endl;

  // then divide SCE corrected track length by the CSDA range
  fRecoBP_CSDAScore = fRecoBP_TrackLength_MuonBDT / CSDARange;

  std::cout << "fRecoBP_CSDAScore: " << fRecoBP_CSDAScore << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool pduneana::NuSelectionAnalysis::IsMuonBDTSignal()
{
  if (fIsReconstructableBeamEvent != 1)
      return false;

  if (!fIsGoodBeamlineTrigger)
      return false;

  if (!fRecoBP_passBeamQualityCuts)
      return false;

  if (!fRecoBP_isTrack)
    return false;

  if (!fRecoBP_passMuonCuts)
    return false;

  if ((fRecoBP_muonChi2 < 0.9) || (fRecoBP_muonChi2 > 1.1))
    return false;

  return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool pduneana::NuSelectionAnalysis::IsMuonBDTBackground()
{
  if (fIsReconstructableBeamEvent != 1)
      return false;

  if (!fIsGoodBeamlineTrigger)
      return false;

  if (!fRecoBP_passBeamQualityCuts)
      return false;

  if (!fRecoBP_isTrack)
    return false;

  if (!fRecoBP_passMuonCuts)
    return true;

  if ((fRecoBP_muonChi2 < 0.9) || (fRecoBP_muonChi2 > 1.1))
    return true;

  return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TVector3 pduneana::NuSelectionAnalysis::ApplyPositionSCECorrection(const geo::Point_t &inputPosition)
{ 
  // To correct for the SCE
  auto sce = lar::providerFrom<spacecharge::SpaceChargeService>();

  geo::Vector_t offset = sce->GetPosOffsets(inputPosition);
  return TVector3(inputPosition.X() - offset.X(), inputPosition.Y() + offset.Y(), inputPosition.Z() + offset.Z());
}


DEFINE_ART_MODULE(pduneana::NuSelectionAnalysis)
