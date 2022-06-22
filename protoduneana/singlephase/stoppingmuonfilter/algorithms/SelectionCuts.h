// Class header for SelectionCuts class

#ifndef SELECTIONCUTS_H
#define SELECTIONCUTS_H

// art
#include "fhiclcpp/ParameterSet.h"

// larsoft
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"

// root
#include "TVector3.h"

// cpp
#include <iostream>

// local
#include "protoduneana/singlephase/stoppingmuonfilter/algorithms/Structs.h"

namespace util{

  class SelectionCuts{

    public:

      // default constructor
      SelectionCuts();

      // default destructor
      ~SelectionCuts(){}

      // configure SelectionCuts with fhicl parameters
      void Configure(fhicl::ParameterSet const& pSet);

      // print the configuration of this SelectionCuts
      void PrintConfiguration();

      // does the track pass the thetaXZ selection?
      // this is done per plane, although it's currently
      // the same cut values for each plane
      std::vector<pdsp::CutResult> IsPassesThetaXZSelection( art::Ptr<recob::Track> thisTrack);

      // does the track pass the thetaYZ selection in
      // each plane?
      std::vector<pdsp::CutResult> IsPassesThetaYZSelection( art::Ptr<recob::Track> thisTrack);

      // does the track pass the minimum distance cut?
      pdsp::MinDistCutResult IsPassesMinimumDistanceCut( art::Ptr< recob::Track > thisTrack, std::vector< art::Ptr< recob::Track > > allTracks );

      // is the minimum hit peak time greater than the 
      // value defined in fhicl?
      pdsp::CutResult IsPassesMinHitPeakTime( std::vector< art::Ptr< recob::Hit > > theseHits );

    private:
      pdsp::AngularLimits _plane0AngularLimits;
      pdsp::AngularLimits _plane1AngularLimits;
      pdsp::AngularLimits _plane2AngularLimits;
      pdsp::BrokenTrackCuts _brokenTracksTightCuts;
      pdsp::BrokenTrackCuts _brokenTracksLooseCuts;
      int _minHitPeakTime;
  
  };

}

#endif
