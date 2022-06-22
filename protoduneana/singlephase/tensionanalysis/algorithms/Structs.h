#ifndef TENSIONSTRUCTS_H
#define TENSIONSTRUCTS_H

// larsoft
#include "larcore/Geometry/Geometry.h"

namespace pdsp{

  enum APASide
  {
    kA,
    kB,
    kUnknown
  };

  // base struct for a simple geometry
  struct APAAndSide
  {
    int APABuildNumber;
    int APAInstallationNumber;
    int APAAnalysisNumber;
    APASide sideToUse = kUnknown;
  };

  // struct for containing calorimetry information
  struct CalorimetryInformation
  {
    float distanceFromSPToCalo;
    float energyDeposition;
    float chargeDeposition;
    float residualRange;
  };

  struct TensionInformation
  {
    float wireSegmentYStart  = -1;
    float wireSegmentYEnd    = -1;
    float wireSegmentZStart  = -1;
    float wireSegmentZEnd    = -1;
    float wireSegmentLength  = -1;
    float wireGeometryYStart = -1;
    float wireGeometryYEnd   = -1;
    float wireGeometryZStart = -1;
    float wireGeometryZEnd   = -1;
    float wireGeometryLength = -1;
    float distSPGeometry     = -1;
    float tension            = -1;
  };

}

#endif
