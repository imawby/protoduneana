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

// Utilities
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "protoduneana/Utilities/ProtoDUNEEmptyEventFinder.h"
#include "protoduneana/Utilities/ProtoDUNEBeamCuts.h"
#include "protoduneana/Utilities/ProtoDUNEBeamlineUtils.h"
#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/Utilities/ProtoDUNECalibration.h"

// ROOT includes
#include "TTree.h"
#include "TFile.h"

// File Handling
#include "cetlib/search_path.h"
#include "cetlib/filesystem.h"

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
  void FillMuonBDTInfo(art::Event const& evt, const recob::PFParticle* beamPFParticle);
  void FillWarwickMVAInfo(art::Event const& evt, const art::Ptr<recob::Track> thisTrack);
  void FillTrackDeflectionInfo(art::Event const& evt, const art::Ptr<recob::Track> thisTrack);
  void FillChi2BDTInfo();
  bool IsMuonBDTSignal();
  bool IsMuonBDTBackground();
  TFile * OpenFile(const std::string filename);

  TTree *fTree;
  TTree *fMuonBDTSignalTree;
  TTree *fMuonBDTBackgroundTree;

  // Run information
  int fRun;
  int fSubrun;
  int fEvent;
  int fIsMC;
  int fIsReconstructableBeamEvent;
  bool fIsGoodBeamlineTrigger;

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

  // Muon BDT variables (can't do a SCE version)
  double fRecoBP_EvalRatio_MuonBDT;
  double fRecoBP_Concentration_MuonBDT;
  double fRecoBP_CoreHaloRatio_MuonBDT;
  double fRecoBP_Conicalness_MuonBDT;

  double fRecoBP_dEdxStart_MuonBDT;
  double fRecoBP_dEdxEnd_MuonBDT;
  double fRecoBP_dEdxEndRatio_MuonBDT;
  double fRecoBP_TrackLength_MuonBDT;
  double fRecoBP_DeflecAngleSD_MuonBDT; //TODO
 
  // PID
  double fRecoBP_protonChi2;
  double fRecoBP_muonChi2;

  // Config
  bool fRecalibrate;
  bool fCheckSlicesForBeam;
  protoana::ProtoDUNECalibration fCalibrationSCE;

  // I hate larsoft
  std::map< int, TProfile* > dEdxTemplates;
  std::string dEdxTemplateName;
  TFile * dEdxTemplateFile;

  // From FHICL file
  std::string fPFParticleLabel;
  std::string fTrackLabel;
  std::string fShowerLabel;
  std::string fPIDLabel; 
  std::string fSCECaloLabel;
  protoana::ProtoDUNEEmptyEventFinder fEmptyEventFinder;
  protoana::ProtoDUNEBeamCuts fBeamCuts;
  protoana::ProtoDUNEBeamlineUtils fBeamlineUtils;
  protoana::ProtoDUNETrackUtils fTrackUtil;

  // Used in the Chi2 PID
  std::vector<double> fCalibrated_dEdX;
  std::vector<double> fResidualRange;
};


pduneana::NuSelectionAnalysis::NuSelectionAnalysis(fhicl::ParameterSet const& pset)
  : EDAnalyzer{pset},
    fRecalibrate(pset.get< bool >("Recalibrate")),
    fCheckSlicesForBeam(pset.get< bool >("CheckSlicesForBeam")),
    fCalibrationSCE(pset.get< fhicl::ParameterSet >("CalibrationParsSCE")),
    dEdxTemplateName(pset.get<std::string>("dEdX_template_name")),
    fPFParticleLabel(pset.get< std::string >("PFParticleLabel")),
    fTrackLabel(pset.get< std::string >("TrackLabel")),
    fShowerLabel(pset.get< std::string >("ShowerLabel")),
    fPIDLabel(pset.get< std::string >("PIDLabel")),
    fSCECaloLabel(pset.get< std::string >("SCECaloLabel")),
    fEmptyEventFinder(pset.get< fhicl::ParameterSet >("EmptyEventFinder")),
    fBeamCuts(pset.get< fhicl::ParameterSet >("BeamCuts")),
    fBeamlineUtils(pset.get<fhicl::ParameterSet>("BeamlineUtils"))
{
    dEdxTemplateFile = OpenFile(dEdxTemplateName);
    dEdxTemplates[ 211 ]  = (TProfile*)dEdxTemplateFile->Get( "dedx_range_pi"  );
    dEdxTemplates[ 321 ]  = (TProfile*)dEdxTemplateFile->Get( "dedx_range_ka"  );
    dEdxTemplates[ 13 ]   = (TProfile*)dEdxTemplateFile->Get( "dedx_range_mu"  );
    dEdxTemplates[ 2212 ] = (TProfile*)dEdxTemplateFile->Get( "dedx_range_pro" );

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
        FillMuonBDTInfo(evt, beamPFParticle);

    // Is shower? Fill Pandrizzle variables

    FillChi2BDTInfo();

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

    fTree->Branch("Run", &fRun);
    fTree->Branch("Subrun", &fSubrun);
    fTree->Branch("Event", &fEvent);
    fTree->Branch("IsMC", &fIsMC);
    fTree->Branch("IsReconstructableBeamEvent", &fIsReconstructableBeamEvent);
    fTree->Branch("IsGoodBeamlineTrigger", &fIsGoodBeamlineTrigger);
    fTree->Branch("checkSlicesForBeam", &fCheckSlicesForBeam);
    fTree->Branch("RecoBP_passBeamQualityCuts", &fRecoBP_passBeamQualityCuts);
    fTree->Branch("RecoBP_isTrack", &fRecoBP_isTrack);
    fTree->Branch("RecoBP_isShower", &fRecoBP_isShower);
    fTree->Branch("RecoBP_momentumByRangeMuonHyp", &fRecoBP_momentumByRangeMuonHyp);
    fTree->Branch("RecoBP_momentumByRangeProtonHyp", &fRecoBP_momentumByRangeProtonHyp);
    fTree->Branch("RecoBP_passMuonCuts", &fRecoBP_passMuonCuts);
    fTree->Branch("RecoBP_passPionCuts", &fRecoBP_passPionCuts);
    fTree->Branch("RecoBP_passElectronCuts", &fRecoBP_passElectronCuts);
    fTree->Branch("RecoBP_passProtonCuts", &fRecoBP_passProtonCuts);
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
  fRecoBP_passBeamQualityCuts = -1;
  fRecoBP_isTrack = false;
  fRecoBP_isShower = false;
  fRecoBP_momentumByRangeMuonHyp = -999.0;
  fRecoBP_momentumByRangeProtonHyp = -999.0;
  fRecoBP_passMuonCuts = false;
  fRecoBP_passPionCuts = false;
  fRecoBP_passElectronCuts = false;
  fRecoBP_passProtonCuts = false;
  fRecoBP_EvalRatio_MuonBDT = -999.0;
  fRecoBP_Concentration_MuonBDT = -999.0;
  fRecoBP_CoreHaloRatio_MuonBDT = -999.0;
  fRecoBP_Conicalness_MuonBDT = -999.0;
  fRecoBP_dEdxStart_MuonBDT = -999.0;
  fRecoBP_dEdxEnd_MuonBDT = -999.0;
  fRecoBP_dEdxEndRatio_MuonBDT = -999.0;
  fRecoBP_TrackLength_MuonBDT = -999.0;
  fRecoBP_DeflecAngleSD_MuonBDT = -999.0;

  fRecoBP_protonChi2 = -999.0;
  fRecoBP_muonChi2 = -999.0;

  fCalibrated_dEdX.clear();
  fResidualRange.clear();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void pduneana::NuSelectionAnalysis::endJob()
{
  dEdxTemplateFile->Close();
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

void pduneana::NuSelectionAnalysis::FillMuonBDTInfo(art::Event const& evt, const recob::PFParticle* beamPFParticle)
{
  art::ValidHandle< std::vector<recob::PFParticle> > pfpListHandle = evt.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);
  const art::FindManyP<recob::Track> findTracks(pfpListHandle, evt, fTrackLabel);
  const std::vector<art::Ptr<recob::Track>> pfpTracks = findTracks.at(beamPFParticle->Self());

  // Check that it is a track
  if (pfpTracks.size() == 0)
    return;

  // Determine if the beam particle is track-like or shower-like
  art::Ptr<recob::Track> artTrack = pfpTracks.at(0);

  //Primary Track Calorimetry
  //SCE-corrected
  protoana::ProtoDUNEPFParticleUtils pfpUtil;
  const recob::Track* thisTrack = pfpUtil.GetPFParticleTrack(*beamPFParticle, evt, fPFParticleLabel, fTrackLabel);
  auto trackCalo = fTrackUtil.GetRecoTrackCalorimetry(*thisTrack, evt, fTrackLabel, fSCECaloLabel);

  //Get the correct plane (collection) because it's not always the same
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
      fRecoBP_TrackLength_MuonBDT = trackCalo[index].Range();

      auto calo_dEdX = trackCalo[index].dEdx();
      auto calo_range = trackCalo[index].ResidualRange();

      // First entry is at the end of the track (LArSoft hates you)
      fRecoBP_dEdxStart_MuonBDT = calo_dEdX.back();
      fRecoBP_dEdxEnd_MuonBDT = calo_dEdX.front();
      fRecoBP_dEdxEndRatio_MuonBDT = calo_dEdX.at(0) / calo_dEdX.at(1);

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

  // Doesn't work with SCE
  //fRecoBP_TrackLength_MuonBDT = thisTrack->Length(); 
  //FillWarwickMVAInfo(evt, artTrack);

  FillTrackDeflectionInfo(evt, artTrack);
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

// Totally copied from Dom, sass retained
void pduneana::NuSelectionAnalysis::FillTrackDeflectionInfo(art::Event const& evt, const art::Ptr<recob::Track> thisTrack)
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
// Get the chi2-based PID for the SCE-corrected beam track
void pduneana::NuSelectionAnalysis::FillChi2BDTInfo()
{
  // First do proton
  std::pair<double, int> protonChi2 = fTrackUtil.Chi2PID(fCalibrated_dEdX, fResidualRange, dEdxTemplates[2212]);
  fRecoBP_protonChi2 = protonChi2.first;
  std::cout << "fRecoBP_protonChi2: " << fRecoBP_protonChi2 << std::endl;

  // Now do muon
  std::pair<double, int> muonChi2 = fTrackUtil.Chi2PID(fCalibrated_dEdX, fResidualRange, dEdxTemplates[13]);
  fRecoBP_muonChi2 = muonChi2.first;
  std::cout << "fRecoBP_muonChi2: " << fRecoBP_muonChi2 << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool pduneana::NuSelectionAnalysis::IsMuonBDTSignal()
{
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
  if (!fRecoBP_isTrack)
    return false;

  if (!fRecoBP_passMuonCuts)
    return true;

  if ((fRecoBP_muonChi2 < 0.9) || (fRecoBP_muonChi2 > 1.1))
    return true;

  return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//To make auto-finding files easier for the dEdX template file
TFile * pduneana::NuSelectionAnalysis::OpenFile(const std::string filename) 
{
  TFile * theFile = 0x0;
  mf::LogInfo("pduneana::OpenFile") << "Searching for " << filename;

  if (cet::file_exists(filename)) 
  {
    mf::LogInfo("pduneana::OpenFile") << "File exists. Opening " << filename;
    theFile = new TFile(filename.c_str());

    if (!theFile ||theFile->IsZombie() || !theFile->IsOpen()) 
    {
      delete theFile;
      theFile = 0x0;
      throw cet::exception("PDSPAnalyzer_module.cc") << "Could not open " << filename;
    }
  }
  else 
  {
    mf::LogInfo("pduneana::OpenFile") << "File does not exist here. Searching FW_SEARCH_PATH";
    cet::search_path sp{"FW_SEARCH_PATH"};
    std::string found_filename;
    auto found = sp.find_file(filename, found_filename);

    if (!found) 
      throw cet::exception("PDSPAnalyzer_module.cc") << "Could not find " << filename;

    mf::LogInfo("pduneana::OpenFile") << "Found file " << found_filename;
    theFile = new TFile(found_filename.c_str());

    if (!theFile ||theFile->IsZombie() || !theFile->IsOpen()) 
    {
      delete theFile;
      theFile = 0x0;
      throw cet::exception("PDSPAnalyzer_module.cc") << "Could not open " << found_filename;
    }
  }

  return theFile;
}

DEFINE_ART_MODULE(pduneana::NuSelectionAnalysis)
