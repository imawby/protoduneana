// Class header for GeometryHelper

#ifndef GEOMETRYHELPER_H
#define GEOMETRYHELPER_H

// framework includes
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// larsoft
#include "larcore/Geometry/Geometry.h"

// cpp
#include <iostream>

// local
#include "protoduneana/singlephase/stoppingmuonfilter/algorithms/Structs.h"

namespace util{

  class GeometryHelper{

    public:

      // default constructor
      GeometryHelper();

      // default destructor
      ~GeometryHelper(){}

      // get total active volume bounded within the detector
      pdsp::GeometryBoundaries GetFullDetectorActiveVolumeBoundaries(art::ServiceHandle<geo::Geometry> geom);

      // get active sub-boundaries taking into account spacing
      // between TPCs
      pdsp::ActiveVolumeBoundaries GetActiveVolumeBoundaries(art::ServiceHandle<geo::Geometry> geom);

  };

}

#endif
