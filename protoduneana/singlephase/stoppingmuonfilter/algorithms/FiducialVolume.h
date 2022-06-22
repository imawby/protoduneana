// Class header for FiducialVolume class

#ifndef FIDUCIALVOLUME_H
#define FIDUCIALVOLUME_H

#include <iostream>
#include "fhiclcpp/ParameterSet.h"
#include "protoduneana/singlephase/stoppingmuonfilter/algorithms/Structs.h"

namespace util{

  class FiducialVolume{

    public:

      // default constructor
      FiducialVolume();

      // default destructor
      ~FiducialVolume(){}

      // configure FiducialVolume with fhicl parameters
      void Configure(fhicl::ParameterSet const& pSet, pdsp::GeometryBoundaries thisGeom);

      // configure multiple fiducial volumes 
      void Configure(fhicl::ParameterSet const& pSet, pdsp::ActiveVolumeBoundaries theseGeom);

      // print the configuration of this FiducialVolume
      void PrintConfiguration();

      // is position in detector fiducial volume
      bool IsInDetectorFV(double x, double y, double z);

      // is position in tpc fiducial volume
      // ths differs from above as it takes into account the 
      // space between tpc active volumes
      bool IsInTpcFV(double x, double y, double z);

    private:
      double _detectorActiveVolumeBorderXLow  = 0;
      double _detectorActiveVolumeBorderXHigh = 0;
      double _detectorActiveVolumeBorderYLow  = 0;
      double _detectorActiveVolumeBorderYHigh = 0;
      double _detectorActiveVolumeBorderZLow  = 0;
      double _detectorActiveVolumeBorderZHigh = 0;

      bool _isTPCFV = false;
      double _tpcActiveVolumeBorderXLow  = 0;
      double _tpcActiveVolumeBorderXHigh = 0;
      double _tpcActiveVolumeBorderYLow  = 0;
      double _tpcActiveVolumeBorderYHigh = 0;
      double _tpcActiveVolumeBorderZLow  = 0;
      double _tpcActiveVolumeBorderZHigh = 0;

      pdsp::GeometryBoundaries _thisGeom;
      pdsp::ActiveVolumeBoundaries _theseGeom;

  };

}

#endif
