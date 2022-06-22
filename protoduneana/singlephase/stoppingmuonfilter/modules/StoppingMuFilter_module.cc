////////////////////////////////////////////////////////////////////////
// Class:       StoppingMuFilter
// Plugin Type: filter (art v3_01_02)
// File:        StoppingMuFilter_module.cc
//
// Generated at Mon Apr 15 12:50:14 2019 by Adam Lister using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

// This code takes in protoDUNE simulated or physics data and filters
// out events with stopping muons.
//
// For events which pass, a set of information is saved for the stopping muons:
// -- recob::PFParticles
// -- recob::Tracks
// -- recob::Showers (for daughters of the muon)
// -- recob::Hits
// -- anab::T0s

// framework
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Persistency/Common/PtrMaker.h"

// larsoft
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/TrackingTypes.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// root
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TH3.h"
#include "TH1.h"

// cpp
#include <memory>

// local
#include "protoduneana/singlephase/stoppingmuonfilter/algorithms/FiducialVolume.h"
#include "protoduneana/singlephase/stoppingmuonfilter/algorithms/GeometryHelper.h"
#include "protoduneana/singlephase/stoppingmuonfilter/algorithms/SelectionCuts.h"
#include "protoduneana/singlephase/stoppingmuonfilter/algorithms/Structs.h"
#include "protoduneana/Utilities/ProtoDUNECalibration.h"

namespace pdsp {
  class StoppingMuFilter;
}


class pdsp::StoppingMuFilter : public art::EDFilter {
  public:
    explicit StoppingMuFilter(fhicl::ParameterSet const& p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    StoppingMuFilter(StoppingMuFilter const&) = delete;
    StoppingMuFilter(StoppingMuFilter&&) = delete;
    StoppingMuFilter& operator=(StoppingMuFilter const&) = delete;
    StoppingMuFilter& operator=(StoppingMuFilter&&) = delete;

    // Required functions.
    bool filter(art::Event& e) override;

    // Selected optional functions.
    void beginJob() override;

    bool beginRun(art::Run& r) override;

    /// function to clear all vectors
    void clearVectors();

    /// apply modified box model
    double getModBoxdEdX(float dqdx, float efield);

    /// get e field
    //double calcEField(float xval, float yval, float zval);

  private:

    // services 
    art::ServiceHandle< art::TFileService > tfs;
    art::ServiceHandle<geo::Geometry> geo;

    // fhicl parameters
    std::string fPfpLabel;
    std::string fTrackLabel;
    std::string fShowerLabel;
    std::string fHitLabel;
    std::string fSpacePointLabel;
    std::string fT0Label;
    std::string fCalorimetryLabel;
    bool fIsUseT0;
    bool fIsThetaXZCutInPlane;
    bool fIsThetaYZCutInPlane;
    bool fThetaXZPlaneToCutOn;
    bool fThetaYZPlaneToCutOn;
    bool fIsUseMinDistCut;
    bool fIsUseMinHitPeakTimeCut;

    // variables for tree
    TH2D* trackStartXZAllTracks;
    TH2D* trackEndXZAllTracks;
    TH2D* trackStartYZAllTracks;
    TH2D* trackEndYZAllTracks;
    TTree* anaTree;
    int run;
    int subRun;
    int event;
    std::vector<int>*    trackNHitsPlane0                          = nullptr;
    std::vector<int>*    trackNHitsPlane1                          = nullptr;
    std::vector<int>*    trackNHitsPlane2                          = nullptr;
    std::vector<double>* trackStartX                               = nullptr;
    std::vector<double>* trackStartY                               = nullptr;
    std::vector<double>* trackStartZ                               = nullptr;
    std::vector<double>* trackEndX                                 = nullptr;
    std::vector<double>* trackEndY                                 = nullptr;
    std::vector<double>* trackEndZ                                 = nullptr;
    std::vector<double>* trackStartXSCECorr                        = nullptr;
    std::vector<double>* trackStartYSCECorr                        = nullptr;
    std::vector<double>* trackStartZSCECorr                        = nullptr;
    std::vector<double>* trackEndXSCECorr                          = nullptr;
    std::vector<double>* trackEndYSCECorr                          = nullptr;
    std::vector<double>* trackEndZSCECorr                          = nullptr;
    std::vector<double>* trackLength                               = nullptr;
    std::vector<double>* trackTheta                                = nullptr;
    std::vector<double>* trackPhi                                  = nullptr;
    std::vector<double>* trackAzimuthal                            = nullptr;
    std::vector<double>* trackZenith                               = nullptr;
    std::vector<double>* trackThetaXZ                              = nullptr;
    std::vector<double>* trackThetaYZ                              = nullptr;
    std::vector<int>* trackHitMinPeakTime                          = nullptr;
    std::vector< std::vector< double > >* trackMinDistLen          = nullptr;
    std::vector< std::vector< double > >* trackMinDistAngle        = nullptr;
    std::vector< std::vector< float > >* trackdEdxByHitPlane0      = nullptr;
    std::vector< std::vector< float > >* trackdQdxByHitPlane0uncal = nullptr;
    std::vector< std::vector< float > >* trackdQdxByHitPlane0      = nullptr;
    std::vector< std::vector< float > >* trackResRangeByHitPlane0  = nullptr;
    std::vector< std::vector< float > >* trackdEdxByHitPlane1      = nullptr;
    std::vector< std::vector< float > >* trackdQdxByHitPlane1uncal = nullptr;
    std::vector< std::vector< float > >* trackdQdxByHitPlane1      = nullptr;
    std::vector< std::vector< float > >* trackResRangeByHitPlane1  = nullptr;
    std::vector< std::vector< float > >* trackdEdxByHitPlane2      = nullptr;
    std::vector< std::vector< float > >* trackdQdxByHitPlane2uncal = nullptr;
    std::vector< std::vector< float > >* trackdQdxByHitPlane2      = nullptr;
    std::vector< std::vector< float > >* trackResRangeByHitPlane2  = nullptr;

    // other vairables
    ::util::GeometryHelper _geomHelper;
    ::util::FiducialVolume _fidVolOuter;
    ::util::FiducialVolume _fidVolInner;
    ::util::SelectionCuts  _selCuts;
    int  isData;
    bool isPass;
    int  hitItP0;
    int  hitItP1;
    int  hitItP2;

};


pdsp::StoppingMuFilter::StoppingMuFilter(fhicl::ParameterSet const& p)
  : EDFilter{p}//,
{
  fhicl::ParameterSet const pLabels  = p.get< fhicl::ParameterSet >("ProducerLabels");
  fhicl::ParameterSet const pCuts    = p.get< fhicl::ParameterSet >("CutValues");
  fhicl::ParameterSet const pOuterFV = p.get< fhicl::ParameterSet >("OuterFV");
  fhicl::ParameterSet const pTpcFV   = p.get< fhicl::ParameterSet >("TpcFV");
  fhicl::ParameterSet const pInnerFV = p.get< fhicl::ParameterSet >("InnerFV");

  fPfpLabel               = pLabels.get< std::string >("PFPLabel");
  fTrackLabel             = pLabels.get< std::string >("TrackLabel");
  fShowerLabel            = pLabels.get< std::string >("ShowerLabel");
  fHitLabel               = pLabels.get< std::string >("HitLabel");
  fSpacePointLabel        = pLabels.get< std::string >("SpacePointLabel");
  fT0Label                = pLabels.get< std::string >("T0Label");
  fCalorimetryLabel       = pLabels.get< std::string >("CalorimetryLabel");
  fIsUseT0                = pCuts.get<bool>("IsUseT0");
  fIsThetaXZCutInPlane    = pCuts.get<bool>("IsThetaXZCutInPlane");
  fIsThetaYZCutInPlane    = pCuts.get<bool>("IsThetaYZCutInPlane");
  fThetaXZPlaneToCutOn    = pCuts.get<int>("ThetaXZPlaneToCutOn");
  fThetaYZPlaneToCutOn    = pCuts.get<int>("ThetaYZPlaneToCutOn");
  fIsUseMinDistCut        = pCuts.get<bool>("IsUseMinDistCut");
  fIsUseMinHitPeakTimeCut = pCuts.get<bool>("IsUseMinHitPeakTimeCut");

  produces< std::vector< recob::PFParticle > >();
  produces< std::vector< recob::Track > >();
  produces< std::vector< recob::Hit > >();
  produces< std::vector< recob::SpacePoint > >();
  produces< std::vector< anab::T0 > >();
  produces< std::vector< anab::Calorimetry > >();
  produces< art::Assns< recob::PFParticle, recob::Track > >();
  produces< art::Assns< recob::Track, recob::Hit > >();
  produces< art::Assns< recob::Hit, recob::SpacePoint > >();
  produces< art::Assns< recob::Track, anab::Calorimetry > >();
  produces< art::Assns< recob::Track, anab::T0 > >();

  // get detector geometry
  pdsp::GeometryBoundaries thisFullDetActiveVol = _geomHelper.GetFullDetectorActiveVolumeBoundaries(geo);
  pdsp::ActiveVolumeBoundaries thisActiveVol    = _geomHelper.GetActiveVolumeBoundaries(geo);

  MF_LOG_VERBATIM("StoppingMuFilter")
    << "Full Detector geometry: "
    << "\n -- lowx : "            << thisFullDetActiveVol.lowx
    << "\n -- highx: "            << thisFullDetActiveVol.highx
    << "\n -- lowy : "            << thisFullDetActiveVol.lowy
    << "\n -- highy: "            << thisFullDetActiveVol.highy
    << "\n -- lowz : "            << thisFullDetActiveVol.lowz
    << "\n -- highz: "            << thisFullDetActiveVol.highz
    << "\nFull Detector has "     << thisActiveVol.numVolumes
    << "subvolumes";

  auto const& thisActiveVolMap = thisActiveVol.activeVolumesMap;

  for (auto tAV_it=thisActiveVolMap.begin(); tAV_it!=thisActiveVolMap.end(); tAV_it++){
    pdsp::GeometryBoundaries thisTpcVol = tAV_it->second;

    MF_LOG_VERBATIM("StoppingMuFilter")
      << "Example TPC geometry for TPC: " << tAV_it->first
      << "\n -- lowx : " << thisTpcVol.lowx
      << "\n -- highx: " << thisTpcVol.highx
      << "\n -- lowy : " << thisTpcVol.lowy
      << "\n -- highy: " << thisTpcVol.highy
      << "\n -- lowz : " << thisTpcVol.lowz
      << "\n -- highz: " << thisTpcVol.highz;
  }

  // define fiducial volumes
  _fidVolOuter.Configure(pOuterFV, 
      thisFullDetActiveVol);
  _fidVolOuter.Configure(pTpcFV, 
      thisActiveVol);
  _fidVolInner.Configure(pInnerFV, 
      thisFullDetActiveVol);

  _fidVolOuter.PrintConfiguration();
  _fidVolInner.PrintConfiguration();

  // configure selection cuts 
  _selCuts.Configure(pCuts);
  _selCuts.PrintConfiguration();

}

bool pdsp::StoppingMuFilter::filter(art::Event& e)
{

  this->clearVectors();

  // sce
  auto const* sce = lar::providerFrom<spacecharge::SpaceChargeService>();

  // get auxiliary information
  run    = e.run();
  subRun = e.subRun();
  event  = e.event();
  isData = e.isRealData();
  isPass = false;

  MF_LOG_VERBATIM("pdsp::StoppingMuFilter::filter")
    << "event: " << run << "." << subRun << "." << event;

  // get handles to information we're interested in
  art::Handle< std::vector< recob::PFParticle > > pfpHandle;
  e.getByLabel(fPfpLabel, pfpHandle);
  std::vector< art::Ptr< recob::PFParticle > > pfpPtrVector;
  art::fill_ptr_vector(pfpPtrVector, pfpHandle);

  art::Handle< std::vector< recob::Track > > trackHandle;
  e.getByLabel(fTrackLabel, trackHandle);
  std::vector< art::Ptr< recob::Track > > trackPtrVector;
  art::fill_ptr_vector(trackPtrVector, trackHandle);

  art::Handle< std::vector< recob::Hit > > hitHandle;
  e.getByLabel(fHitLabel, hitHandle);

  // associations
  art::FindManyP< recob::Track >      tracksFromPfps(pfpHandle, e, fTrackLabel);
  art::FindManyP< anab::T0 >          t0sFromPfps(pfpHandle, e, fT0Label);
  art::FindManyP< recob::Hit >        hitsFromTracks(trackHandle, e, fTrackLabel);
  art::FindManyP< recob::SpacePoint > spacePointFromHits(hitHandle, e, fSpacePointLabel);
  art::FindManyP< anab::Calorimetry > caloFromTracks(trackHandle, e, fCalorimetryLabel);

  // define the collections of objects to be put in the event
  std::unique_ptr< std::vector< recob::PFParticle > >              pfpCollection( new std::vector< recob::PFParticle > );
  std::unique_ptr< std::vector< recob::Track > >                   trackCollection( new std::vector< recob::Track > );
  std::unique_ptr< std::vector< recob::Hit > >                     hitCollection( new std::vector< recob::Hit > );
  std::unique_ptr< std::vector< recob::SpacePoint > >              spacePointCollection( new std::vector< recob::SpacePoint > );
  std::unique_ptr< std::vector< anab::Calorimetry > >              calorimetryCollection( new std::vector< anab::Calorimetry > );
  std::unique_ptr< std::vector< anab::T0 > >                       t0Collection( new std::vector< anab::T0 > );
  std::unique_ptr< art::Assns< recob::PFParticle, recob::Track > > pfpTrackAssn( new art::Assns< recob::PFParticle, recob::Track > );
  std::unique_ptr< art::Assns< recob::Track, recob::Hit > >        trackHitAssn( new art::Assns< recob::Track, recob::Hit> );
  std::unique_ptr< art::Assns< recob::Track, anab::Calorimetry > > trackCalorimetryAssn( new art::Assns< recob::Track, anab::Calorimetry > );
  std::unique_ptr< art::Assns< recob::Hit, recob::SpacePoint > >   hitSpacePointAssn( new art::Assns< recob::Hit, recob::SpacePoint> );
  std::unique_ptr< art::Assns< recob::Track, anab::T0 > >          trackT0Assn( new art::Assns< recob::Track, anab::T0 > );

  // define for creation of ptrs
  art::PtrMaker<recob::PFParticle> makePfpPtr(e);
  art::PtrMaker<recob::Track>      makeTrackPtr(e);
  art::PtrMaker<anab::T0>          makeT0Ptr(e);
  art::PtrMaker<recob::Hit>        makeHitPtr(e);
  art::PtrMaker<recob::SpacePoint> makeSpacePointPtr(e);
  art::PtrMaker<anab::Calorimetry> makeCalorimetryPtr(e);

  MF_LOG_VERBATIM("pdsp::StoppingMuFilter::filter")
    << "looping PFP vector which has size " << pfpPtrVector.size();

  // check to make sure run number is calibrated

    for ( size_t i = 0; i < pfpPtrVector.size(); i++){

      art::Ptr< recob::PFParticle > thisPfp = pfpPtrVector.at(i);

      if ((tracksFromPfps.at(thisPfp.key())).size() == 0) continue;

      // find particles which are identified as track-like and have a single 
      // daughter
      if (thisPfp->PdgCode() == 13){

        art::Ptr< recob::Track > thisTrack = (tracksFromPfps.at(thisPfp.key())).at(0);

        trackStartXZAllTracks ->Fill(thisTrack->Start().Z(), thisTrack->Start().X());
        trackEndXZAllTracks   ->Fill(thisTrack->End().Z()  , thisTrack->End().X());
        trackStartYZAllTracks ->Fill(thisTrack->Start().Z(), thisTrack->Start().Y());
        trackEndYZAllTracks   ->Fill(thisTrack->End().Z()  , thisTrack->End().Y());

        // associated t0 must exist
        if (t0sFromPfps.at(thisPfp.key()).size() == 0 && fIsUseT0 == true) continue;

        MF_LOG_VERBATIM("pdsp::StoppingMuFilter::filter")
          << "-- track has associated t0!";

        // if t0s are being used then get the associated T0, 
        // if not, then create a dud to put in the file. 
        // Not the best handling but for now it'll do.

        art::Ptr< anab::T0 > thisT0;
        if (fIsUseT0)
          thisT0 = (t0sFromPfps.at(thisPfp.key())).at(0);

        // track must meet fiducial requirements (enter the detector and stop)
        // this goes through the space charge service, the calirbation module doesn't deal with spatial offsets
        if (  !_fidVolOuter.IsInDetectorFV(thisTrack -> Start().X() - sce -> GetPosOffsets(thisTrack -> Start()).X(),
                                           thisTrack -> Start().Y() - sce -> GetPosOffsets(thisTrack -> Start()).Y(),
                                           thisTrack -> Start().Z() - sce -> GetPosOffsets(thisTrack -> Start()).Z())
            && _fidVolInner.IsInDetectorFV(thisTrack -> End().X()   - sce -> GetPosOffsets(thisTrack -> End()).X() ,
                                           thisTrack -> End().Y()   - sce -> GetPosOffsets(thisTrack -> End()).Y() ,
                                           thisTrack -> End().Z()   - sce -> GetPosOffsets(thisTrack -> End()).Z())
            && _fidVolOuter.IsInTpcFV(     thisTrack -> Start().X() - sce -> GetPosOffsets(thisTrack -> Start()).X(),
                                           thisTrack -> Start().Y() - sce -> GetPosOffsets(thisTrack -> Start()).Y(),
                                           thisTrack -> Start().Z() - sce -> GetPosOffsets(thisTrack -> Start()).Z())
            && _fidVolOuter.IsInTpcFV(     thisTrack -> End().X()   - sce -> GetPosOffsets(thisTrack -> End()).X()  ,
                                           thisTrack -> End().Y()   - sce -> GetPosOffsets(thisTrack -> End()).Y()  ,
                                           thisTrack -> End().Z()   - sce -> GetPosOffsets(thisTrack -> End()).Z())){

          MF_LOG_VERBATIM("pdsp::StoppingMuFilter::filter")
            << "  -- Track meets fiducial requirements";

          std::vector< pdsp::CutResult > isPassThetaXZ = _selCuts.IsPassesThetaXZSelection(thisTrack);
          std::vector< pdsp::CutResult > isPassThetaYZ = _selCuts.IsPassesThetaYZSelection(thisTrack);


          // if there's a demand that the track pass the theta_xz 
          // cut in a specific plane then apply that
          if (fIsThetaXZCutInPlane){
            if (!(isPassThetaXZ.at(fThetaXZPlaneToCutOn).result)){
              MF_LOG_VERBATIM("pdsp::StoppingMuFilter::filter")
                << "  -- Track Fails Theta XZ";
              continue;
            }
          }

          // if there's a demand that the track pass the theta_yz 
          // cut in a specific plane then apply that
          if (fIsThetaYZCutInPlane){
            if (!(isPassThetaYZ.at(fThetaYZPlaneToCutOn).result)){
              MF_LOG_VERBATIM("pdsp::StoppingMuFilter::filter")
                << "  -- Track Fails Theta YZ";
              continue;
            }
          }

          // this demands that the minimum distance and angle cut be met
          // i.e. this removes broken tracks
          pdsp::MinDistCutResult thisMinDistCutResult = _selCuts.IsPassesMinimumDistanceCut(thisTrack, trackPtrVector); 

          if (fIsUseMinDistCut){
            if (!thisMinDistCutResult.result){
              MF_LOG_VERBATIM("pdsp::StoppingMuFilter::filter")
                << "  -- Track Fails min dist cut";
              continue;
            }
          }

          // get the hits associated with this track and make the 
          // demand that the minimum time be greater than some value
          std::vector< art::Ptr< recob::Hit > > theseHits = hitsFromTracks.at(thisTrack.key());

          pdsp::CutResult thisMinHitCutResult = _selCuts.IsPassesMinHitPeakTime(theseHits);

          if (fIsUseMinHitPeakTimeCut){
            if(!thisMinHitCutResult.result){
              MF_LOG_VERBATIM("pdsp::StoppingMuFilter::filter")
                << "  -- Track Fails min hit peak time";
              continue;
            }
          }

          // If we reach this point, the event is selected, and we're going 
          // to go ahead and save the information to the tree
          MF_LOG_VERBATIM("StoppingMuFilter")
            << "  -- Event selected";

          trackStartX         -> push_back(thisTrack -> Start().X());
          trackStartY         -> push_back(thisTrack -> Start().Y());
          trackStartZ         -> push_back(thisTrack -> Start().Z());
          trackEndX           -> push_back(thisTrack -> End().X());
          trackEndY           -> push_back(thisTrack -> End().Y());
          trackEndZ           -> push_back(thisTrack -> End().Z());
          trackStartXSCECorr  -> push_back(thisTrack -> Start().X() - sce -> GetPosOffsets(thisTrack -> Start()).X());
          trackStartYSCECorr  -> push_back(thisTrack -> Start().Y() - sce -> GetPosOffsets(thisTrack -> Start()).Y());
          trackStartZSCECorr  -> push_back(thisTrack -> Start().Z() - sce -> GetPosOffsets(thisTrack -> Start()).Z());
          trackEndXSCECorr    -> push_back(thisTrack -> End().X()   - sce -> GetPosOffsets(thisTrack -> End()).X());
          trackEndYSCECorr    -> push_back(thisTrack -> End().Y()   - sce -> GetPosOffsets(thisTrack -> End()).Y());
          trackEndZSCECorr    -> push_back(thisTrack -> End().Z()   - sce -> GetPosOffsets(thisTrack -> End()).Z());
          trackLength         -> push_back(thisTrack -> Length());
          trackTheta          -> push_back(thisTrack -> Theta());
          trackPhi            -> push_back(thisTrack -> Phi());
          trackAzimuthal      -> push_back(thisTrack -> AzimuthAngle());
          trackZenith         -> push_back(thisTrack -> ZenithAngle());
          trackThetaXZ        -> push_back((isPassThetaXZ.at(2)).value);
          trackThetaYZ        -> push_back((isPassThetaYZ.at(2)).value);
          trackHitMinPeakTime -> push_back(thisMinHitCutResult.value);
          trackMinDistLen     -> push_back(thisMinDistCutResult.lenValue);
          trackMinDistAngle   -> push_back(thisMinDistCutResult.angleValue);

          // now deal with getting information for associations
          recob::PFParticle pfpForCollection = *(thisPfp.get());
          recob::Track trackForCollection    = *(thisTrack.get());
          anab::T0 t0ForCollection(-1, -1, -1);
          if (fIsUseT0)
            t0ForCollection = *(thisT0.get());

          pfpCollection  ->push_back(pfpForCollection);
          trackCollection->push_back(trackForCollection);
          t0Collection   ->push_back(t0ForCollection);

          art::Ptr<recob::PFParticle> pfpForCollectionPtr   = makePfpPtr(pfpCollection->size()-1);
          art::Ptr<recob::Track>      trackForCollectionPtr = makeTrackPtr(trackCollection->size()-1);
          art::Ptr<anab::T0>          t0ForCollectionPtr    = makeT0Ptr(t0Collection->size() -1);

          // get calorimetry information
          std::vector< art::Ptr< anab::Calorimetry > > theseCalo = caloFromTracks.at(thisTrack.key());
          std::vector< art::Ptr< anab::Calorimetry > > calorimetryPtrCollection;

          for (size_t i_cal = 0; i_cal < theseCalo.size(); i_cal++){

            art::Ptr< anab::Calorimetry > thisCalo     = theseCalo.at(i_cal);
            anab::Calorimetry calorimetryForCollection = *(thisCalo.get());

            // calibrate dqdx points
            std::vector< float > trackdQdx = thisCalo->dQdx();
            std::vector< float > trackdEdx = thisCalo->dEdx();

            // first deal with getting the stuff for association creation
            calorimetryForCollection.fdQdx = trackdQdx;
            calorimetryForCollection.fdEdx = trackdEdx;

            calorimetryCollection->push_back(calorimetryForCollection);
            art::Ptr< anab::Calorimetry > calorimetryForCollectionPtr = makeCalorimetryPtr(calorimetryCollection->size() -1);
            calorimetryPtrCollection.push_back(calorimetryForCollectionPtr);

            // and now fill stuff for the output tree
            if (thisCalo->PlaneID().Plane == 0){
              trackdEdxByHitPlane0     ->push_back(calorimetryForCollection.dEdx());
              trackdQdxByHitPlane0uncal->push_back(thisCalo->dQdx());
              trackdQdxByHitPlane0     ->push_back(calorimetryForCollection.dQdx());
              trackResRangeByHitPlane0 ->push_back(calorimetryForCollection.ResidualRange());
            }
            if (thisCalo->PlaneID().Plane == 1){
              trackdEdxByHitPlane1     ->push_back(calorimetryForCollection.dEdx());
              trackdQdxByHitPlane1uncal->push_back(thisCalo->dQdx());
              trackdQdxByHitPlane1     ->push_back(calorimetryForCollection.dQdx());
              trackResRangeByHitPlane1 ->push_back(calorimetryForCollection.ResidualRange());
            }
            if (thisCalo->PlaneID().Plane == 2){
              trackdEdxByHitPlane2     ->push_back(calorimetryForCollection.dEdx());
              trackdQdxByHitPlane2uncal->push_back(thisCalo->dQdx());
              trackdQdxByHitPlane2     ->push_back(calorimetryForCollection.dQdx());
              trackResRangeByHitPlane2 ->push_back(calorimetryForCollection.ResidualRange());
            }

          }


          util::CreateAssn(e, trackForCollectionPtr, calorimetryPtrCollection, *trackCalorimetryAssn); 
          std::vector< art::Ptr<recob::Hit> >hitPtrCollection;

          // now get hit information for creating the associations

          hitItP0 = 0;
          hitItP1 = 0;
          hitItP2 = 0;
          for (size_t iH = 0; iH < theseHits.size(); iH++){

            if (theseHits[iH]->View() == 0)
              hitItP0++; 
            if (theseHits[iH]->View() == 1)
              hitItP1++; 
            if (theseHits[iH]->View() == 2)
              hitItP2++; 

            hitCollection->push_back(*((theseHits.at(iH)).get()));
            art::Ptr<recob::Hit> hitForCollectionPtr = makeHitPtr(hitCollection->size() -1);
            hitPtrCollection.push_back(hitForCollectionPtr);

            if (spacePointFromHits.at((theseHits.at(iH)).key()).size() != 0){
              art::Ptr<recob::SpacePoint> thisSpacePoint = spacePointFromHits.at((theseHits.at(iH)).key()).at(0);

              spacePointCollection->push_back(*(thisSpacePoint.get()));
              art::Ptr<recob::SpacePoint> spacePointForCollectionPtr = makeSpacePointPtr(spacePointCollection->size() -1);

              // need to create a 1-to-1 assciation for spacepoints to hits
              util::CreateAssn(e, spacePointForCollectionPtr, hitForCollectionPtr, *hitSpacePointAssn); 
            }
          }

          trackNHitsPlane0->push_back(hitItP0);
          trackNHitsPlane1->push_back(hitItP1);
          trackNHitsPlane2->push_back(hitItP2);

          util::CreateAssn(e, trackForCollectionPtr, pfpForCollectionPtr  , *pfpTrackAssn);
          util::CreateAssn(e, trackForCollectionPtr, hitPtrCollection     , *trackHitAssn);
          util::CreateAssn(e, t0ForCollectionPtr   , trackForCollectionPtr, *trackT0Assn);


          anaTree->Fill();

          isPass = true;
        }
      }
    }
  //}
  //else MF_LOG_VERBATIM("pdsp::stoppingmufilter::filter") << "Run not calibrated";

  // now put data products in the event
  e.put(std::move(pfpCollection));
  e.put(std::move(trackCollection));
  e.put(std::move(hitCollection));
  e.put(std::move(spacePointCollection));
  e.put(std::move(calorimetryCollection));
  e.put(std::move(t0Collection));
  e.put(std::move(pfpTrackAssn));
  e.put(std::move(trackHitAssn));
  e.put(std::move(hitSpacePointAssn));
  e.put(std::move(trackCalorimetryAssn));
  e.put(std::move(trackT0Assn));
  return isPass;

}

void pdsp::StoppingMuFilter::beginJob()
{  
  trackStartXZAllTracks = tfs->make<TH2D>("trackStartXZAllTracks", 
      ";Track Start Z (cm);Track Start X (cm)", 
      100, -100, 800, 100, -450, 450);

  trackEndXZAllTracks   = tfs->make<TH2D>("trackEndXZAllTracks", 
      ";Track End Z (cm);Track End X (cm)", 
      100, -100, 800, 100, -450, 450);

  trackStartYZAllTracks = tfs->make<TH2D>("trackStartYZAllTracks", 
      ";Track Start Z (cm);Track Start Y (cm)", 
      100, -100, 800, 100, -100, 800);

  trackEndYZAllTracks   = tfs->make<TH2D>("trackEndYZAllTracks", 
      ";Track End Z (cm);Track End Y (cm)", 
      100, -100, 800, 100, -100, 800);

  anaTree = tfs->make<TTree>("analysis_tree" , "analysis tree");
  anaTree->Branch("run"                      , &run);
  anaTree->Branch("subRun"                   , &subRun);
  anaTree->Branch("event"                    , &event);
  anaTree->Branch("trackNHitsPlane0"         , "std::vector<int>"                     , &trackNHitsPlane0);
  anaTree->Branch("trackNHitsPlane1"         , "std::vector<int>"                     , &trackNHitsPlane1);
  anaTree->Branch("trackNHitsPlane2"         , "std::vector<int>"                     , &trackNHitsPlane2);
  anaTree->Branch("trackStartX"              , "std::vector<double>"                  , &trackStartX);
  anaTree->Branch("trackStartY"              , "std::vector<double>"                  , &trackStartY);
  anaTree->Branch("trackStartZ"              , "std::vector<double>"                  , &trackStartZ);
  anaTree->Branch("trackEndX"                , "std::vector<double>"                  , &trackEndX);
  anaTree->Branch("trackEndY"                , "std::vector<double>"                  , &trackEndY);
  anaTree->Branch("trackEndZ"                , "std::vector<double>"                  , &trackEndZ);
  anaTree->Branch("trackStartXSCECorr"       , "std::vector<double>"                  , &trackStartXSCECorr);
  anaTree->Branch("trackStartYSCECorr"       , "std::vector<double>"                  , &trackStartYSCECorr);
  anaTree->Branch("trackStartZSCECorr"       , "std::vector<double>"                  , &trackStartZSCECorr);
  anaTree->Branch("trackEndXSCECorr"         , "std::vector<double>"                  , &trackEndXSCECorr);
  anaTree->Branch("trackEndYSCECorr"         , "std::vector<double>"                  , &trackEndYSCECorr);
  anaTree->Branch("trackEndZSCECorr"         , "std::vector<double>"                  , &trackEndZSCECorr);
  anaTree->Branch("trackLength"              , "std::vector<double>"                  , &trackLength);
  anaTree->Branch("trackTheta"               , "std::vector<double>"                  , &trackTheta);
  anaTree->Branch("trackPhi"                 , "std::vector<double>"                  , &trackPhi);
  anaTree->Branch("trackAzimuthal"           , "std::vector<double>"                  , &trackAzimuthal);
  anaTree->Branch("trackZenith"              , "std::vector<double>"                  , &trackZenith);
  anaTree->Branch("trackThetaXZ"             , "std::vector<double>"                  , &trackThetaXZ);
  anaTree->Branch("trackThetaYZ"             , "std::vector<double>"                  , &trackThetaYZ);
  anaTree->Branch("trackMinDistLen"          , "std::vector< std::vector< double > >" , &trackMinDistLen);
  anaTree->Branch("trackMinDistAngle"        , "std::vector< std::vector< double > >" , &trackMinDistAngle);

  anaTree->Branch("trackHitMinPeakTime"      , "std::vector< int >"                   , &trackHitMinPeakTime);
  anaTree->Branch("trackdEdxByHitPlane0"     , "std::vector< std::vector< float > >"  , &trackdEdxByHitPlane0);
  anaTree->Branch("trackdQdxByHitPlane0uncal", "std::vector< std::vector< float > >"  , &trackdQdxByHitPlane0uncal);
  anaTree->Branch("trackdQdxByHitPlane0"     , "std::vector< std::vector< float > >"  , &trackdQdxByHitPlane0);
  anaTree->Branch("trackResRangeByHitPlane0" , "std::vector< std::vector< float > >"  , &trackResRangeByHitPlane0);
  anaTree->Branch("trackdEdxByHitPlane1"     , "std::vector< std::vector< float > >"  , &trackdEdxByHitPlane1);
  anaTree->Branch("trackdQdxByHitPlane1uncal", "std::vector< std::vector< float > >"  , &trackdQdxByHitPlane1uncal);
  anaTree->Branch("trackdQdxByHitPlane1"     , "std::vector< std::vector< float > >"  , &trackdQdxByHitPlane1);
  anaTree->Branch("trackResRangeByHitPlane1" , "std::vector< std::vector< float > >"  , &trackResRangeByHitPlane1);
  anaTree->Branch("trackdEdxByHitPlane2"     , "std::vector< std::vector< float > >"  , &trackdEdxByHitPlane2);
  anaTree->Branch("trackdQdxByHitPlane2uncal", "std::vector< std::vector< float > >"  , &trackdQdxByHitPlane2uncal);
  anaTree->Branch("trackdQdxByHitPlane2"     , "std::vector< std::vector< float > >"  , &trackdQdxByHitPlane2);
  anaTree->Branch("trackResRangeByHitPlane2" , "std::vector< std::vector< float > >"  , &trackResRangeByHitPlane2);

}

void pdsp::StoppingMuFilter::clearVectors()
{
  trackNHitsPlane0          -> resize(0);
  trackNHitsPlane1          -> resize(0); 
  trackNHitsPlane2          -> resize(0);
  trackStartX               -> resize(0);
  trackStartY               -> resize(0);
  trackStartZ               -> resize(0);
  trackEndX                 -> resize(0);
  trackEndY                 -> resize(0);
  trackEndZ                 -> resize(0);
  trackStartXSCECorr        -> resize(0);
  trackStartYSCECorr        -> resize(0);
  trackStartZSCECorr        -> resize(0);
  trackEndXSCECorr          -> resize(0);
  trackEndYSCECorr          -> resize(0);
  trackEndZSCECorr          -> resize(0);
  trackLength               -> resize(0);
  trackTheta                -> resize(0);
  trackPhi                  -> resize(0);
  trackAzimuthal            -> resize(0);
  trackZenith               -> resize(0);
  trackThetaXZ              -> resize(0);
  trackThetaYZ              -> resize(0);
  trackHitMinPeakTime       -> resize(0);
  trackMinDistLen           -> resize(0);
  trackMinDistAngle         -> resize(0);
  trackdEdxByHitPlane0      -> resize(0);
  trackdQdxByHitPlane0uncal -> resize(0);
  trackdQdxByHitPlane0      -> resize(0);
  trackResRangeByHitPlane0  -> resize(0);
  trackdEdxByHitPlane1      -> resize(0);
  trackdQdxByHitPlane1      -> resize(0);
  trackResRangeByHitPlane1  -> resize(0);
  trackdEdxByHitPlane2      -> resize(0);
  trackdQdxByHitPlane2      -> resize(0);
  trackResRangeByHitPlane2  -> resize(0);
}

bool pdsp::StoppingMuFilter::beginRun(art::Run& r)
{
  return true;
}


DEFINE_ART_MODULE(pdsp::StoppingMuFilter)
