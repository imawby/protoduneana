#ifndef SELECTIONCUTS_CXX
#define SELECTIONCUTS_CXX

#include "SelectionCuts.h"

namespace util{

  SelectionCuts::SelectionCuts()
  {
  }

  void SelectionCuts::Configure(fhicl::ParameterSet const &pSet)
  {

    fhicl::ParameterSet const &pSetAngles = pSet.get< fhicl::ParameterSet >("AngularCuts");
    fhicl::ParameterSet const &pSetLowX   = pSetAngles.get< fhicl::ParameterSet >("LowX");
    fhicl::ParameterSet const &pSetHighX  = pSetAngles.get< fhicl::ParameterSet >("HighX");

    std::vector< std::vector<double> > angularCutsLowXThetaXZ = pSetLowX.get< std::vector<std::vector< double > > >("ThetaXZ");
    std::vector< std::vector<double> > angularCutsLowXThetaYZ = pSetLowX.get< std::vector<std::vector< double > > >("ThetaYZ");
    std::vector< std::vector<double> > angularCutsHighXThetaXZ = pSetHighX.get< std::vector<std::vector< double > > >("ThetaXZ");
    std::vector< std::vector<double> > angularCutsHighXThetaYZ = pSetHighX.get< std::vector<std::vector< double > > >("ThetaYZ");

    _plane0AngularLimits.lowxThetaXZ.first   = angularCutsLowXThetaXZ.at(0).at(0);
    _plane0AngularLimits.lowxThetaXZ.second  = angularCutsLowXThetaXZ.at(0).at(1);
    _plane0AngularLimits.lowxThetaYZ.first   = angularCutsLowXThetaYZ.at(0).at(0);
    _plane0AngularLimits.lowxThetaYZ.second  = angularCutsLowXThetaYZ.at(0).at(1);
    _plane0AngularLimits.highxThetaXZ.first  = angularCutsHighXThetaXZ.at(0).at(0);
    _plane0AngularLimits.highxThetaXZ.second = angularCutsHighXThetaXZ.at(0).at(1);
    _plane0AngularLimits.highxThetaYZ.first  = angularCutsHighXThetaYZ.at(0).at(0);
    _plane0AngularLimits.highxThetaYZ.second = angularCutsHighXThetaYZ.at(0).at(1);

    _plane1AngularLimits.lowxThetaXZ.first   = angularCutsLowXThetaXZ.at(1).at(0);
    _plane1AngularLimits.lowxThetaXZ.second  = angularCutsLowXThetaXZ.at(1).at(1);
    _plane1AngularLimits.lowxThetaYZ.first   = angularCutsLowXThetaYZ.at(1).at(0);
    _plane1AngularLimits.lowxThetaYZ.second  = angularCutsLowXThetaYZ.at(1).at(1);
    _plane1AngularLimits.highxThetaXZ.first  = angularCutsHighXThetaXZ.at(1).at(0);
    _plane1AngularLimits.highxThetaXZ.second = angularCutsHighXThetaXZ.at(1).at(1);
    _plane1AngularLimits.highxThetaYZ.first  = angularCutsHighXThetaYZ.at(1).at(0);
    _plane1AngularLimits.highxThetaYZ.second = angularCutsHighXThetaYZ.at(1).at(1);

    _plane2AngularLimits.lowxThetaXZ.first   = angularCutsLowXThetaXZ.at(2).at(0);
    _plane2AngularLimits.lowxThetaXZ.second  = angularCutsLowXThetaXZ.at(2).at(1);
    _plane2AngularLimits.lowxThetaYZ.first   = angularCutsLowXThetaYZ.at(2).at(0);
    _plane2AngularLimits.lowxThetaYZ.second  = angularCutsLowXThetaYZ.at(2).at(1);
    _plane2AngularLimits.highxThetaXZ.first  = angularCutsHighXThetaXZ.at(2).at(0);
    _plane2AngularLimits.highxThetaXZ.second = angularCutsHighXThetaXZ.at(2).at(1);
    _plane2AngularLimits.highxThetaYZ.first  = angularCutsHighXThetaYZ.at(2).at(0);
    _plane2AngularLimits.highxThetaYZ.second = angularCutsHighXThetaYZ.at(2).at(1);

    fhicl::ParameterSet const &pSetBrokenTracks = pSet.get< fhicl::ParameterSet >("BrokenTracksCuts");
    fhicl::ParameterSet const &pSetTightCuts = pSetBrokenTracks.get< fhicl::ParameterSet >("TightCuts");
    fhicl::ParameterSet const &pSetLooseCuts = pSetBrokenTracks.get< fhicl::ParameterSet >("LooseCuts");
    _brokenTracksLooseCuts.dist  = pSetLooseCuts.get< double >("Distance");
    _brokenTracksLooseCuts.angle = pSetLooseCuts.get< double >("Angle");
    _brokenTracksTightCuts.dist  = pSetTightCuts.get< double >("Distance");
    _brokenTracksTightCuts.angle = pSetTightCuts.get< double >("Angle");

    _minHitPeakTime = pSet.get< int >("HitMinPeakTime");

  }

  void SelectionCuts::PrintConfiguration(){

    std::cout << "Printing angular selection:"     << std::endl;
    std::cout << " -- plane 0, low x, theta_xz:  " << _plane0AngularLimits.lowxThetaXZ.first  << ", " << _plane0AngularLimits.lowxThetaXZ.second  << std::endl;
    std::cout << " -- plane 0, low x, theta_yz:  " << _plane0AngularLimits.lowxThetaYZ.first  << ", " << _plane0AngularLimits.lowxThetaYZ.second  << std::endl;
    std::cout << " -- plane 0, high x, theta_xz: " << _plane0AngularLimits.highxThetaXZ.first << ", " << _plane0AngularLimits.highxThetaXZ.second << std::endl;
    std::cout << " -- plane 0, high x, theta_yz: " << _plane0AngularLimits.highxThetaYZ.first << ", " << _plane0AngularLimits.highxThetaYZ.second << std::endl;
    std::cout << " -- plane 1, low x, theta_xz:  " << _plane1AngularLimits.lowxThetaXZ.first  << ", " << _plane1AngularLimits.lowxThetaXZ.second  << std::endl;
    std::cout << " -- plane 1, low x, theta_yz:  " << _plane1AngularLimits.lowxThetaYZ.first  << ", " << _plane1AngularLimits.lowxThetaYZ.second  << std::endl;
    std::cout << " -- plane 1, high x, theta_xz: " << _plane1AngularLimits.highxThetaXZ.first << ", " << _plane1AngularLimits.highxThetaXZ.second << std::endl;
    std::cout << " -- plane 1, high x, theta_yz: " << _plane1AngularLimits.highxThetaYZ.first << ", " << _plane1AngularLimits.highxThetaYZ.second << std::endl;
    std::cout << " -- plane 2, low x, theta_xz:  " << _plane2AngularLimits.lowxThetaXZ.first  << ", " << _plane2AngularLimits.lowxThetaXZ.second  << std::endl;
    std::cout << " -- plane 2, low x, theta_yz:  " << _plane2AngularLimits.lowxThetaYZ.first  << ", " << _plane2AngularLimits.lowxThetaYZ.second  << std::endl;
    std::cout << " -- plane 2, high x, theta_xz: " << _plane2AngularLimits.highxThetaXZ.first << ", " << _plane2AngularLimits.highxThetaXZ.second << std::endl;
    std::cout << " -- plane 2, high x, theta_yz: " << _plane2AngularLimits.highxThetaYZ.first << ", " << _plane2AngularLimits.highxThetaYZ.second << std::endl;

    std::cout << "Printing broken track cuts:"     << std::endl;
    std::cout << " -- Tight cuts (distance, angle):  " << _brokenTracksTightCuts.dist  << ", " << _brokenTracksTightCuts.angle  << std::endl;
    std::cout << " -- Loose cuts (distance, angle):  " << _brokenTracksLooseCuts.dist  << ", " << _brokenTracksLooseCuts.angle  << std::endl;
  
    std::cout << "Printing minimum hit peak time cut: " << _minHitPeakTime << std::endl;

  }

  std::vector< pdsp::CutResult > SelectionCuts::IsPassesThetaXZSelection( art::Ptr<recob::Track> thisTrack){
  
    std::vector< pdsp::CutResult > passesSelection;
    pdsp::CutResult cutResultPlane0;
    pdsp::CutResult cutResultPlane1;
    pdsp::CutResult cutResultPlane2;

    TVector3 startPosition  = thisTrack->Vertex<TVector3>();
    TVector3 startDirection = thisTrack->VertexDirection<TVector3>();

    double thetaXZ = std::atan2(startDirection.X(), startDirection.Z());
   
    // low x
    if (startPosition.X() < 0){
      
      // plane 0
      if (   std::abs(thetaXZ*180./TMath::Pi()) < _plane0AngularLimits.lowxThetaXZ.first
          || std::abs(thetaXZ*180./TMath::Pi()) > _plane0AngularLimits.lowxThetaXZ.second){
        cutResultPlane0.result = true;
        cutResultPlane0.value = thetaXZ;
        passesSelection.push_back(cutResultPlane0);
      }
      else {
        cutResultPlane0.result = false;
        cutResultPlane0.value = thetaXZ;
        passesSelection.push_back(cutResultPlane0);
      }

      // plane 1
      if (   std::abs(thetaXZ*180./TMath::Pi()) < _plane1AngularLimits.lowxThetaXZ.first
          || std::abs(thetaXZ*180./TMath::Pi()) > _plane1AngularLimits.lowxThetaXZ.second){
        cutResultPlane1.result = true;
        cutResultPlane1.value = thetaXZ;
        passesSelection.push_back(cutResultPlane1);
      }
      else {
        cutResultPlane1.result = false;
        cutResultPlane1.value = thetaXZ;
        passesSelection.push_back(cutResultPlane1);
      }

      // plane 2
      if (   std::abs(thetaXZ*180./TMath::Pi()) < _plane2AngularLimits.lowxThetaXZ.first
          || std::abs(thetaXZ*180./TMath::Pi()) > _plane2AngularLimits.lowxThetaXZ.second){
        cutResultPlane2.result = true;
        cutResultPlane2.value = thetaXZ;
        passesSelection.push_back(cutResultPlane2);
      }
      else {
        cutResultPlane2.result = false;
        cutResultPlane2.value = thetaXZ;
        passesSelection.push_back(cutResultPlane2);
      }

    }
    // high x
    else {

      // plane 0
      if (   std::abs(thetaXZ*180./TMath::Pi()) < _plane0AngularLimits.highxThetaXZ.first
          || std::abs(thetaXZ*180./TMath::Pi()) > _plane0AngularLimits.highxThetaXZ.second){
        cutResultPlane0.result = true;
        cutResultPlane0.value = thetaXZ;
        passesSelection.push_back(cutResultPlane0);
      }
      else {
        cutResultPlane0.result = false;
        cutResultPlane0.value = thetaXZ;
        passesSelection.push_back(cutResultPlane0);
      }

      // plane 1
      if (   std::abs(thetaXZ*180./TMath::Pi()) < _plane1AngularLimits.highxThetaXZ.first
          || std::abs(thetaXZ*180./TMath::Pi()) > _plane1AngularLimits.highxThetaXZ.second){
        cutResultPlane1.result = true;
        cutResultPlane1.value = thetaXZ;
        passesSelection.push_back(cutResultPlane1);
      }
      else {
        cutResultPlane1.result = false;
        cutResultPlane1.value = thetaXZ;
        passesSelection.push_back(cutResultPlane1);
      }

      // plane 2
      if (   std::abs(thetaXZ*180./TMath::Pi()) < _plane2AngularLimits.highxThetaXZ.first
          || std::abs(thetaXZ*180./TMath::Pi()) > _plane2AngularLimits.highxThetaXZ.second){
        cutResultPlane2.result = true;
        cutResultPlane2.value = thetaXZ;
        passesSelection.push_back(cutResultPlane2);
      }
      else {
        cutResultPlane2.result = false;
        cutResultPlane2.value = thetaXZ;
        passesSelection.push_back(cutResultPlane2);
      }
    }

    return passesSelection;

  }

  std::vector< pdsp::CutResult > SelectionCuts::IsPassesThetaYZSelection( art::Ptr<recob::Track> thisTrack){
 
    std::vector< pdsp::CutResult > passesSelection;
    pdsp::CutResult cutResultPlane0;
    pdsp::CutResult cutResultPlane1;
    pdsp::CutResult cutResultPlane2;

    TVector3 startPosition  = thisTrack->Vertex<TVector3>();
    TVector3 startDirection = thisTrack->VertexDirection<TVector3>();

    double thetaYZ = std::atan2(startDirection.Y(), startDirection.Z());
   
    // low x
    if (startPosition.X() < 0){
      
      // plane 0
      if (   thetaYZ*180./TMath::Pi() < _plane0AngularLimits.lowxThetaYZ.first
          || thetaYZ*180./TMath::Pi() > _plane0AngularLimits.lowxThetaYZ.second){
        cutResultPlane0.result = true;
        cutResultPlane0.value = thetaYZ;
        passesSelection.push_back(cutResultPlane0);
      }
      else {
        cutResultPlane0.result = false;
        cutResultPlane0.value = thetaYZ;
        passesSelection.push_back(cutResultPlane0);
      }

      // plane 1
      if (   thetaYZ*180./TMath::Pi() < _plane1AngularLimits.lowxThetaYZ.first
          || thetaYZ*180./TMath::Pi() > _plane1AngularLimits.lowxThetaYZ.second) {
        cutResultPlane1.result = true;
        cutResultPlane1.value = thetaYZ;
        passesSelection.push_back(cutResultPlane1);
      }
      else {
        cutResultPlane1.result = false;
        cutResultPlane1.value = thetaYZ;
        passesSelection.push_back(cutResultPlane1);
      }

      // plane 2
      if (   std::abs(thetaYZ*180./TMath::Pi()) < _plane2AngularLimits.lowxThetaYZ.first
          || std::abs(thetaYZ*180./TMath::Pi()) > _plane2AngularLimits.lowxThetaYZ.second){
        cutResultPlane2.result = true;
        cutResultPlane2.value = thetaYZ;
        passesSelection.push_back(cutResultPlane2);
      }
      else {
        cutResultPlane2.result = false;
        cutResultPlane2.value = thetaYZ;
        passesSelection.push_back(cutResultPlane2);
      }

    }
    // high x
    else {

      // plane 0
      if (   thetaYZ*180./TMath::Pi() < _plane0AngularLimits.highxThetaYZ.first
          || thetaYZ*180./TMath::Pi() > _plane0AngularLimits.highxThetaYZ.second){
        cutResultPlane0.result = true;
        cutResultPlane0.value = thetaYZ;
        passesSelection.push_back(cutResultPlane0);
      }
      else {
        cutResultPlane0.result = false;
        cutResultPlane0.value = thetaYZ;
        passesSelection.push_back(cutResultPlane0);
      }

      // plane 1
      if (   thetaYZ*180./TMath::Pi() < _plane1AngularLimits.highxThetaYZ.first
          || thetaYZ*180./TMath::Pi() > _plane1AngularLimits.highxThetaYZ.second){
        cutResultPlane1.result = true;
        cutResultPlane1.value = thetaYZ;
        passesSelection.push_back(cutResultPlane1);
      }
      else {
        cutResultPlane1.result = false;
        cutResultPlane1.value = thetaYZ;
        passesSelection.push_back(cutResultPlane1);
      }

      // plane 2
      if (   std::abs(thetaYZ*180./TMath::Pi()) < _plane2AngularLimits.highxThetaYZ.first
          || std::abs(thetaYZ*180./TMath::Pi()) > _plane2AngularLimits.highxThetaYZ.second){
        cutResultPlane2.result = true;
        cutResultPlane2.value = thetaYZ;
        passesSelection.push_back(cutResultPlane2);
      }
      else {
        cutResultPlane2.result = false;
        cutResultPlane2.value = thetaYZ;
        passesSelection.push_back(cutResultPlane2);
      }
    }

    return passesSelection;

  }

  pdsp::MinDistCutResult SelectionCuts::IsPassesMinimumDistanceCut( art::Ptr< recob::Track > thisTrack, std::vector< art::Ptr< recob::Track > > allTracks ){

    pdsp::MinDistCutResult cutResults;

    double starty1 = thisTrack->Start().Y();
    double startz1 = thisTrack->Start().Z();
    double endy1 = thisTrack->End().Y();
    double endz1 = thisTrack->End().Z();
  
    TVector3 startdir1 = thisTrack->StartDirection<TVector3>();
    TVector3 enddir1 = thisTrack->EndDirection<TVector3>();

    // calculate intercept 

    for (size_t i = 0; i < allTracks.size(); i++){

      art::Ptr< recob::Track > arbTrack = allTracks.at(i);
      
      // ignore thisTrack when looping allTracks
      if (arbTrack.key() == thisTrack.key()) continue;

      double starty2 = arbTrack->Start().Y();
      double startz2 = arbTrack->Start().Z();
      double endy2 = arbTrack->End().Y();
      double endz2 = arbTrack->End().Z();
  
      // calculate distance between start/ends of this track 
      // and the start/end points of every other 
      // going to want to do this only in y and z directions 
      // since we don't know the t0 of every track in the event

      double delLenstart1start2 = std::sqrt(std::pow(starty1-starty2,2) + std::pow(startz1-startz2,2));
      double delLenstart1end2   = std::sqrt(std::pow(starty1-endy2,2) + std::pow(startz1-endz2,2));
      double delLenend1start2   = std::sqrt(std::pow(endy1-starty2,2) + std::pow(endz1-startz2,2));
      double delLenend1end2     = std::sqrt(std::pow(endy1-endy2,2) + std::pow(endz1-endz2,2));

      // also calculate the angles between the startdir/enddir of
      // this track and that of every other track

      TVector3 startdir2 = arbTrack->StartDirection<TVector3>();
      TVector3 enddir2   = arbTrack->EndDirection<TVector3>();
     
      double delAngstart1start2 = std::abs(startdir1.X()*startdir2.X() + startdir1.Y()*startdir2.Y() + startdir1.Z()*startdir2.Z()); 
      double delAngstart1end2 = std::abs(startdir1.X()*enddir2.X() + startdir1.Y()*enddir2.Y() + startdir1.Z()*enddir2.Z()); 
      double delAngend1start2 = std::abs(enddir1.X()*startdir2.X() + enddir1.Y()*startdir2.Y() + enddir1.Z()*startdir2.Z()); 
      double delAngend1end2 = std::abs(enddir1.X()*enddir2.X() + enddir1.Y()*enddir2.Y() + enddir1.Z()*enddir2.Z()); 

      // we only care about the distance between track start and end points if they
      // are within some angular tolerance. There's a set of tight cuts and loose cuts
      // which are set to be fhicl configurable

      // check against tight cuts
      if ((    delLenstart1start2 < _brokenTracksTightCuts.dist
            && delAngstart1start2 > _brokenTracksTightCuts.angle)
          || ( delLenstart1end2   < _brokenTracksTightCuts.dist
            && delAngstart1end2   > _brokenTracksTightCuts.angle)
          || ( delLenend1start2   < _brokenTracksTightCuts.dist
            && delAngend1start2   > _brokenTracksTightCuts.angle)
          || ( delLenend1end2     < _brokenTracksTightCuts.dist
            && delAngend1end2     > _brokenTracksTightCuts.angle))
        cutResults.result = false;
      // check against loose cuts
      else if ((delLenstart1start2 < _brokenTracksLooseCuts.dist
            &&  delAngstart1start2 > _brokenTracksLooseCuts.angle)
          || (  delLenstart1end2   < _brokenTracksLooseCuts.dist
            &&  delAngstart1end2   > _brokenTracksLooseCuts.angle)
          || (  delLenend1start2   < _brokenTracksLooseCuts.dist
            &&  delAngend1start2   > _brokenTracksLooseCuts.angle)
          || (  delLenend1end2     < _brokenTracksLooseCuts.dist
            &&  delAngend1end2     > _brokenTracksLooseCuts.angle))
        cutResults.result = false;

      cutResults.lenValue.push_back(delLenstart1start2);
      cutResults.lenValue.push_back(delLenstart1end2);
      cutResults.lenValue.push_back(delLenend1start2);
      cutResults.lenValue.push_back(delLenend1end2);
      cutResults.angleValue.push_back(delAngstart1start2);
      cutResults.angleValue.push_back(delAngstart1end2);
      cutResults.angleValue.push_back(delAngend1start2);
      cutResults.angleValue.push_back(delAngend1end2);

    }

    return cutResults;
  }
      
  pdsp::CutResult SelectionCuts::IsPassesMinHitPeakTime( std::vector< art::Ptr< recob::Hit > > theseHits ){

    pdsp::CutResult thisCutResult;

    double minHitPeakTime = std::numeric_limits<int>::max();

    for (size_t i = 0; i < theseHits.size(); i++){

      art::Ptr< recob::Hit > thisHit = theseHits.at(i);

      if (thisHit->PeakTime() < minHitPeakTime)
        minHitPeakTime = thisHit->PeakTime();

    }


    if (minHitPeakTime < _minHitPeakTime){
      thisCutResult.result = false;
      thisCutResult.value = minHitPeakTime;
    }
    else {  
      thisCutResult.result = true;
      thisCutResult.value = minHitPeakTime;
    }

  return thisCutResult;

  }

}

#endif
