#ifndef GEOMETRYHELPER_CXX
#define GEOMETRYHELPER_CXX

#include "GeometryHelper.h"

namespace util{

  GeometryHelper::GeometryHelper()
  {
  }

  pdsp::GeometryBoundaries GeometryHelper::GetFullDetectorActiveVolumeBoundaries(art::ServiceHandle<geo::Geometry> geom){

    // note on protodune geometry:
    // x = 0 is the center cathode, the anodes are at +/- ~360
    // y = 0 is at the bottom of the TPC frame
    // z = 0 is at the upstream edge of the upstream APA
    //
    // in order to get min/max, initialise each variable to be the
    // middle of the detector

    pdsp::GeometryBoundaries detectorGeometry;
    detectorGeometry.lowx = 0;
    detectorGeometry.highx = 0;
    detectorGeometry.lowy = 300;
    detectorGeometry.highy = 300;
    detectorGeometry.lowz = 350;
    detectorGeometry.highz = 350;

    for (geo::TPCID const& tID: geom->IterateTPCIDs()) {

      geo::TPCGeo const& TPC = geom->TPC(tID);

      geo::BoxBoundedGeo const activeTpcBox = TPC.ActiveBoundingBox();

      double boxXLow  = activeTpcBox.MinX();
      double boxXHigh = activeTpcBox.MaxX();
      double boxYLow  = activeTpcBox.MinY();
      double boxYHigh = activeTpcBox.MaxY();
      double boxZLow  = activeTpcBox.MinZ();
      double boxZHigh = activeTpcBox.MaxZ();

      if (std::abs(boxXHigh - boxXLow) < 100) 
        continue;

      if (boxXLow  < detectorGeometry.lowx)
        detectorGeometry.lowx = boxXLow;
      if (boxXHigh > detectorGeometry.highx)
        detectorGeometry.highx = boxXHigh;
      if (boxYLow  < detectorGeometry.lowy)
        detectorGeometry.lowy = boxYLow;
      if (boxYHigh > detectorGeometry.highy)
        detectorGeometry.highy = boxYHigh;
      if (boxZLow  < detectorGeometry.lowz)
        detectorGeometry.lowz = boxZLow;
      if (boxZHigh > detectorGeometry.highz)
        detectorGeometry.highz = boxZHigh;

    }

    return detectorGeometry;

  }

  pdsp::ActiveVolumeBoundaries GeometryHelper::GetActiveVolumeBoundaries(art::ServiceHandle<geo::Geometry> geom){

    pdsp::ActiveVolumeBoundaries activeVolBoundaries;
    
    for (geo::TPCID const& tID: geom->IterateTPCIDs()) {

      geo::TPCGeo const& TPC = geom->TPC(tID);

      geo::BoxBoundedGeo const activeTpcBox = TPC.ActiveBoundingBox();
      
      // get GeometryBoundaries for this TPC
      pdsp::GeometryBoundaries tpcGeometry;
      tpcGeometry.lowx = activeTpcBox.MinX();
      tpcGeometry.highx = activeTpcBox.MaxX();
      tpcGeometry.lowy = activeTpcBox.MinY();
      tpcGeometry.highy = activeTpcBox.MaxY();
      tpcGeometry.lowz = activeTpcBox.MinZ();
      tpcGeometry.highz = activeTpcBox.MaxZ();

      // There are seprate TPCs which are simulated behind the primary TPCs
      // which are characterised by short (~20 cm) dimensions in the x direction
      //
      // I don't understand this but for now catch these and remove them
      if (std::abs(tpcGeometry.highx - tpcGeometry.lowx) < 100) 
        continue;

      activeVolBoundaries.activeVolumesMap.insert( std::make_pair((int)tID.TPC, tpcGeometry) );
      activeVolBoundaries.numVolumes++;

    }

    return activeVolBoundaries;

  }

}

#endif
