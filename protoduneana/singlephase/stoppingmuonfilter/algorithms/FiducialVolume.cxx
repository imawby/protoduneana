#ifndef FIDUCIALVOLUME_CXX
#define FIDUCIALVOLUME_CXX

#include "FiducialVolume.h"

namespace util{

  FiducialVolume::FiducialVolume()
  {
  }

  void FiducialVolume::Configure(fhicl::ParameterSet const &pSet, pdsp::GeometryBoundaries thisGeom)
  {

    _detectorActiveVolumeBorderXLow  = pSet.get<double>("BorderXLow");
    _detectorActiveVolumeBorderXHigh = pSet.get<double>("BorderXHigh");
    _detectorActiveVolumeBorderYLow  = pSet.get<double>("BorderYLow");
    _detectorActiveVolumeBorderYHigh = pSet.get<double>("BorderYHigh");
    _detectorActiveVolumeBorderZLow  = pSet.get<double>("BorderZLow");
    _detectorActiveVolumeBorderZHigh = pSet.get<double>("BorderZHigh");

    _thisGeom = thisGeom;

  }

  void FiducialVolume::Configure(fhicl::ParameterSet const &pSet, pdsp::ActiveVolumeBoundaries theseGeom)
  {

    _tpcActiveVolumeBorderXLow  = pSet.get<double>("BorderXLow");
    _tpcActiveVolumeBorderXHigh = pSet.get<double>("BorderXHigh");
    _tpcActiveVolumeBorderYLow  = pSet.get<double>("BorderYLow");
    _tpcActiveVolumeBorderYHigh = pSet.get<double>("BorderYHigh");
    _tpcActiveVolumeBorderZLow  = pSet.get<double>("BorderZLow");
    _tpcActiveVolumeBorderZHigh = pSet.get<double>("BorderZHigh");

    _theseGeom = theseGeom;

    _isTPCFV = true;
  }


  void FiducialVolume::PrintConfiguration(){

    std::cout << "Printing Detector Fiducial Volume:" << std::endl;
    std::cout << " -- X, Low:  " << _detectorActiveVolumeBorderXLow  << std::endl;
    std::cout << " -- X, High: " << _detectorActiveVolumeBorderXHigh << std::endl;
    std::cout << " -- Y, Low:  " << _detectorActiveVolumeBorderYLow  << std::endl;
    std::cout << " -- Y, High: " << _detectorActiveVolumeBorderYHigh << std::endl;
    std::cout << " -- Z, Low:  " << _detectorActiveVolumeBorderZLow  << std::endl;
    std::cout << " -- Z, High: " << _detectorActiveVolumeBorderZHigh << std::endl;

    std::cout << "Printing bounding box:" << std::endl;
    std::cout << " -- X, Low:  " << _thisGeom.lowx  + _detectorActiveVolumeBorderXLow  << std::endl;
    std::cout << " -- X, High: " << _thisGeom.highx - _detectorActiveVolumeBorderXHigh << std::endl;
    std::cout << " -- Y, Low:  " << _thisGeom.lowy  + _detectorActiveVolumeBorderYLow  << std::endl;
    std::cout << " -- Y, High: " << _thisGeom.highy - _detectorActiveVolumeBorderYHigh << std::endl;
    std::cout << " -- Z, Low:  " << _thisGeom.lowz  + _detectorActiveVolumeBorderZLow  << std::endl;
    std::cout << " -- Z, High: " << _thisGeom.highz - _detectorActiveVolumeBorderZHigh << std::endl;

    if (_isTPCFV){
      std::cout << "Printing TPC Fiducial Volume:" << std::endl;
      std::cout << " -- X, Low:  " << _tpcActiveVolumeBorderXLow  << std::endl;
      std::cout << " -- X, High: " << _tpcActiveVolumeBorderXHigh << std::endl;
      std::cout << " -- Y, Low:  " << _tpcActiveVolumeBorderYLow  << std::endl;
      std::cout << " -- Y, High: " << _tpcActiveVolumeBorderYHigh << std::endl;
      std::cout << " -- Z, Low:  " << _tpcActiveVolumeBorderZLow  << std::endl;
      std::cout << " -- Z, High: " << _tpcActiveVolumeBorderZHigh << std::endl;
    }
  }

  bool FiducialVolume::IsInDetectorFV(double x, double y, double z){

    if (   x > (_thisGeom.lowx  + _detectorActiveVolumeBorderXLow) 
        && x < (_thisGeom.highx - _detectorActiveVolumeBorderXHigh)
        && y > (_thisGeom.lowy  + _detectorActiveVolumeBorderYLow) 
        && y < (_thisGeom.highy - _detectorActiveVolumeBorderYHigh)
        && z > (_thisGeom.lowz  + _detectorActiveVolumeBorderZLow) 
        && z < (_thisGeom.highz - _detectorActiveVolumeBorderZHigh))
    {
      return true;
    }
    else {
      return false;
    }
  }

  bool FiducialVolume::IsInTpcFV(double x, double y, double z){

    std::map<int, pdsp::GeometryBoundaries> thisMap = _theseGeom.activeVolumesMap;

    for (std::map<int, pdsp::GeometryBoundaries>::iterator it = thisMap.begin(); it != thisMap.end(); ++it){

      pdsp::GeometryBoundaries thisGeom = it->second;
      
      if ( x > (thisGeom.lowx  + _tpcActiveVolumeBorderXLow) 
        && x < (thisGeom.highx - _tpcActiveVolumeBorderXHigh)
        && y > (thisGeom.lowy  + _tpcActiveVolumeBorderYLow) 
        && y < (thisGeom.highy - _tpcActiveVolumeBorderYHigh)
        && z > (thisGeom.lowz  + _tpcActiveVolumeBorderZLow) 
        && z < (thisGeom.highz - _tpcActiveVolumeBorderZHigh))
      {
        return true;
      }
      else continue;

    }

    return false;

  }


}

#endif
