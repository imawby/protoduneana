#ifndef STOPMUSTRUCTS_H
#define STOPMUSTRUCTS_H

// larsoft
#include "larcore/Geometry/Geometry.h"

// cpp
#include <map>

namespace pdsp{

  // base struct for a simple geometry
  struct GeometryBoundaries
  {
    double lowx;
    double highx;
    double lowy;
    double highy;
    double lowz;
    double highz;
  };

  // allows multiple active volumes within detector
  struct ActiveVolumeBoundaries
  {
    int numVolumes=0;
    std::map<int, pdsp::GeometryBoundaries> activeVolumesMap;
  };

  // struct to hold angular cuts for low x (x<0) and highx (x>0)
  struct AngularLimits
  {
    std::pair<double,double> lowxThetaXZ  = std::make_pair(0.,0.);
    std::pair<double,double> lowxThetaYZ  = std::make_pair(0.,0.);
    std::pair<double,double> highxThetaYZ = std::make_pair(0.,0.);
    std::pair<double,double> highxThetaXZ = std::make_pair(0.,0.);
  };

  struct BrokenTrackCuts
  {
    double dist  = 0;
    double angle = 0;
  };

  struct CutResult
  {
    bool result = false;
    double value = 0.;
  };

  struct MinDistCutResult
  {
    bool result = true;
    std::vector<double> lenValue;
    std::vector<double> angleValue;
  };

  struct CalibrationParameters
  {
    double calibFactor          = -1;
    double normalisationFactor  = -1;
    std::string yzCorrFactorLoc = "noLocation";
    std::string xCorrFactorLoc  = "noLocation";
  };

}

#endif
