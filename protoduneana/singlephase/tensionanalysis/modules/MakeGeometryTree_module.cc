////////////////////////////////////////////////////////////////////////
// Class:       MakeGeometryTree
// Plugin Type: analyzer (Unknown Unknown)
// File:        MakeGeometryTree_module.cc
//
// Generated at Mon Jan 10 15:28:40 2022 by Adam Lister using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "cetlib/search_path.h"

// larsoft
#include "larcore/Geometry/Geometry.h"
#include "dune-raw-data/Services/ChannelMap/PdspChannelMapService.h"
#include "protoduneana/singlephase/tensionanalysis/algorithms/Structs.h"

// root
#include "TFile.h"
#include "TTree.h"

// cpp
#include <memory>

class MakeGeometryTree;


class MakeGeometryTree : public art::EDAnalyzer {
public:
  explicit MakeGeometryTree(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MakeGeometryTree(MakeGeometryTree const&) = delete;
  MakeGeometryTree(MakeGeometryTree&&) = delete;
  MakeGeometryTree& operator=(MakeGeometryTree const&) = delete;
  MakeGeometryTree& operator=(MakeGeometryTree&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;

    /// given an analysis APA number, return a struct containing the 
    /// build APA number along with the side which faces the cathode
    /// and therefore should be used for analysis
    pdsp::APAAndSide getAPAInfoFromAnaAPANumber(int anaAPANumber);

private:
    art::ServiceHandle< art::TFileService > tfs;
    art::ServiceHandle< dune::PdspChannelMapService > channelMap;
    geo::GeometryCore const* geom = lar::providerFrom<geo::Geometry>();

  // Declare member data here.
    // geometry tree information
    TTree* geometryTree;                                ///< tree which stores geometry information
    int channelNumber;                                  ///< channel number
    int channelNumberAssociatedFEMB;                    ///< channel number associated FEMB number (1-4)
    int channelNumberAssociatedWIB;                     ///< channel number associated Warm Intefcae Board
    int channelNumberAssociatedFEMBIndex;               ///< FEMB index calculation is taken from DataUtils, 
                                                        ///< from the dunetpc repo, uses above 2 variabales in 
                                                        ///< calculation
    int channelNumberAssociatedFEMBChannel;             ///< channel on the FEMB
    int channelNumberAssociatedWires;                   ///< number of wires associated with the channel
    int channelAssociatedWiresPlane;                    ///< wire plane associated with this channel
    int channelAssociatedAPABuildNumber;                ///< APA build number associated with this channel

    std::vector<double>* channelAssociatedWiresLength     = nullptr; ///< geom length of wires associated with channel
    std::vector<double>* channelAssociatedWiresStartX     = nullptr; ///< geom start x position (drift direction) of wires
    std::vector<double>* channelAssociatedWiresEndX       = nullptr; ///< geom end x position (drift direction) of wires
    std::vector<double>* channelAssociatedWiresStartY     = nullptr; ///< geom start y position (vertical) of wires
    std::vector<double>* channelAssociatedWiresEndY       = nullptr; ///< geom end y position (vertical) of wires
    std::vector<double>* channelAssociatedWiresStartZ     = nullptr; ///< geom start z position (beam direction) of wires
    std::vector<double>* channelAssociatedWiresEndZ       = nullptr; ///< geom end z position (beam direction) of wires
    std::vector<double>* channelAssociatedSegmentsStartY  = nullptr; ///< excel start y position (vertical) of wires
    std::vector<double>* channelAssociatedSegmentsEndY    = nullptr; ///< excel end y position (vertical) of wires
    std::vector<double>* channelAssociatedSegmentsStartZ  = nullptr; ///< excel start z position (beam direction) of wires
    std::vector<double>* channelAssociatedSegmentsEndZ    = nullptr; ///< excel end z position (beam direction) of wires

    // root files to read in
    TFile* tensionsFile;

    // TTrees and maps for individual planes
    TTree* treeXLayerUS001;
    TTree* treeULayerUS001;
    TTree* treeVLayerUS001;
    TTree* treeXLayerUS002;
    TTree* treeULayerUS002;
    TTree* treeVLayerUS002;
    TTree* treeXLayerUS003;
    TTree* treeULayerUS003;
    TTree* treeVLayerUS003;
    TTree* treeXLayerUS004;
    TTree* treeULayerUS004;
    TTree* treeVLayerUS004;
    TTree* treeXLayerUK001;
    TTree* treeULayerUK001;
    TTree* treeVLayerUK001;
    TTree* treeXLayerUK002;
    TTree* treeULayerUK002;
    TTree* treeVLayerUK002;



};


MakeGeometryTree::MakeGeometryTree(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void MakeGeometryTree::analyze(art::Event const& e)
{
  // Implementation of required member function here.
}

void MakeGeometryTree::beginJob()
{
  // Implementation of optional member function here.
  MF_LOG_DEBUG("TensionAnalysis")
    << "-- begin TensionAnalysis::beginJob";

  geometryTree = tfs->make<TTree>("geometry_tree"           , "geometry tree");
  geometryTree->Branch("channelNumber"                      , &channelNumber);
  geometryTree->Branch("channelNumberAssociatedFEMB"        , &channelNumberAssociatedFEMB);
  geometryTree->Branch("channelNumberAssociatedWIB"         , &channelNumberAssociatedWIB);
  geometryTree->Branch("channelNumberAssociatedFEMBIndex"   , &channelNumberAssociatedFEMBIndex);
  geometryTree->Branch("channelNumberAssociatedFEMBChannel" , &channelNumberAssociatedFEMBChannel);
  geometryTree->Branch("channelNumberAssociatedWires"       , &channelNumberAssociatedWires);
  geometryTree->Branch("channelAssociatedWiresPlane"        , &channelAssociatedWiresPlane);
  geometryTree->Branch("channelAssociatedAPABuildNumber"    , &channelAssociatedAPABuildNumber);

  geometryTree->Branch("channelAssociatedWiresLength"    , "std::vector<double>" , &channelAssociatedWiresLength);
  geometryTree->Branch("channelAssociatedWiresStartX"    , "std::vector<double>" , &channelAssociatedWiresStartX);
  geometryTree->Branch("channelAssociatedWiresEndX"      , "std::vector<double>" , &channelAssociatedWiresEndX);
  geometryTree->Branch("channelAssociatedWiresStartY"    , "std::vector<double>" , &channelAssociatedWiresStartY);
  geometryTree->Branch("channelAssociatedWiresEndY"      , "std::vector<double>" , &channelAssociatedWiresEndY);
  geometryTree->Branch("channelAssociatedWiresStartZ"    , "std::vector<double>" , &channelAssociatedWiresStartZ);
  geometryTree->Branch("channelAssociatedWiresEndZ"      , "std::vector<double>" , &channelAssociatedWiresEndZ);
  geometryTree->Branch("channelAssociatedSegmentsStartY" , "std::vector<double>" , &channelAssociatedSegmentsStartY);
  geometryTree->Branch("channelAssociatedSegmentsEndY"   , "std::vector<double>" , &channelAssociatedSegmentsEndY);
  geometryTree->Branch("channelAssociatedSegmentsStartZ" , "std::vector<double>" , &channelAssociatedSegmentsStartZ);
  geometryTree->Branch("channelAssociatedSegmentsEndZ"   , "std::vector<double>" , &channelAssociatedSegmentsEndZ);

  // read in ROOT trees containing wire information
  std::string tensionPath;
  cet::search_path sp("FW_SEARCH_PATH");
  if(!sp.find_file("tensionanalysis/data/tension_measurements_mod.root", tensionPath)){
    throw cet::exception("FileError")
      << "Cannot find tension_measurements_mod.root file "
      << " bail ungracefully\n\n"
      << __FILE__ << ":" << __LINE__;
  }
  tensionsFile = new TFile(tensionPath.c_str(), "read");

  treeXLayerUS001 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US001_XLAYER");
  treeULayerUS001 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US001_ULAYER");
  treeVLayerUS001 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US001_VLAYER");
  treeXLayerUS002 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US002_XLAYER");
  treeULayerUS002 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US002_ULAYER");
  treeVLayerUS002 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US002_VLAYER");
  treeXLayerUS003 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US003_XLAYER");
  treeULayerUS003 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US003_ULAYER");
  treeVLayerUS003 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US003_VLAYER");
  treeXLayerUS004 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US004_XLAYER");
  treeULayerUS004 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US004_ULAYER");
  treeVLayerUS004 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US004_VLAYER");
  treeXLayerUK001 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_UK001_XLAYER");
  treeULayerUK001 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_UK001_ULAYER");
  treeVLayerUK001 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_UK001_VLAYER");
  treeXLayerUK002 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_UK002_XLAYER");
  treeULayerUK002 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_UK002_ULAYER");
  treeVLayerUK002 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_UK002_VLAYER");

  MF_LOG_DEBUG("TensionAnalysis")
    << "Found tensions file and pulled out trees";

  // make a vector of TTrees for easiness
  std::vector< TTree* > treeVec = {
    treeXLayerUS001,
    treeULayerUS001,
    treeVLayerUS001,
    treeXLayerUS002,
    treeULayerUS002,
    treeVLayerUS002,
    treeXLayerUS003,
    treeULayerUS003,
    treeVLayerUS003,
    treeXLayerUS004,
    treeULayerUS004,
    treeVLayerUS004,
    treeXLayerUK001,
    treeULayerUK001,
    treeVLayerUK001,
    treeXLayerUK002,
    treeULayerUK002,
    treeVLayerUK002
  };

    MF_LOG_DEBUG("TensionAnalysis")
    << "Filling geometry tree";
  // loop all channels and get the number of associated wires, along with their lengths
  for (int i_chan = 0; i_chan < 15360; i_chan++){
    channelAssociatedWiresLength   ->resize(0);
    channelAssociatedWiresStartX   ->resize(0);
    channelAssociatedWiresStartY   ->resize(0);
    channelAssociatedWiresStartZ   ->resize(0);
    channelAssociatedWiresEndX     ->resize(0);
    channelAssociatedWiresEndY     ->resize(0);
    channelAssociatedWiresEndZ     ->resize(0);
    channelAssociatedSegmentsStartY->resize(0);
    channelAssociatedSegmentsStartZ->resize(0);
    channelAssociatedSegmentsEndY  ->resize(0);
    channelAssociatedSegmentsEndZ  ->resize(0);

    std::vector<geo::WireID> wireIDs = geom->ChannelToWire(i_chan);

    channelNumber                      = i_chan;
    channelNumberAssociatedWires       = (int)wireIDs.size();
    channelNumberAssociatedFEMB        = (int)channelMap->FEMBFromOfflineChannel(channelNumber);
    channelNumberAssociatedWIB         = (int)channelMap->WIBFromOfflineChannel(channelNumber); 
    channelNumberAssociatedFEMBIndex   = ((channelNumberAssociatedWIB*4)+(channelNumberAssociatedFEMB-1));
    channelNumberAssociatedFEMBChannel = (int)channelMap->FEMBChannelFromOfflineChannel(channelNumber);
    channelAssociatedWiresPlane        = wireIDs.at(0).Plane;

    // get TPC ID information
    pdsp::APAAndSide thisAPAAndSide = this->getAPAInfoFromAnaAPANumber(wireIDs.at(0).TPC);
    channelAssociatedAPABuildNumber = thisAPAAndSide.APABuildNumber;

    for (size_t i_wire = 0; i_wire < wireIDs.size(); i_wire++){
      geo::WireID thisWireID = wireIDs.at(i_wire);
      geo::WireGeo const& thisWireGeo = geom->Wire(thisWireID);
      channelAssociatedWiresLength-> push_back(thisWireGeo.Length());
      channelAssociatedWiresStartX-> push_back(thisWireGeo.GetStart()[0]);
      channelAssociatedWiresEndX->   push_back(thisWireGeo.GetEnd()[0]);
      channelAssociatedWiresStartY-> push_back(thisWireGeo.GetStart()[1]);
      channelAssociatedWiresEndY->   push_back(thisWireGeo.GetEnd()[1]);
      channelAssociatedWiresStartZ-> push_back(thisWireGeo.GetStart()[2]);
      channelAssociatedWiresEndZ->   push_back(thisWireGeo.GetEnd()[2]);
    }

    for (size_t itree = 0; itree < treeVec.size(); itree++){

      TTree* t = (TTree*)treeVec.at(itree);

      int   segmentSideAChannel;
      int   segmentSideBChannel;
      int   segmentNumber;
      float segmentStartY;
      float segmentEndY;
      float segmentStartZ;
      float segmentEndZ;

      // numbers from the excel sheets use a different co-ordinate system
      // the APA measurements were taken with the APA in the horizontal orientation
      // meaning that x measures the the longest edge, and y measures the shortest edge
      // 
      // when querying the geometry, z goes along the direction of the shortest edge and
      // y measures the longest edge, so rename variables here.

      t->SetBranchAddress("side_a_channel_number" , &segmentSideAChannel);
      t->SetBranchAddress("side_b_channel_number" , &segmentSideBChannel);
      t->SetBranchAddress("segment_number"        , &segmentNumber);
      t->SetBranchAddress("x_start"               , &segmentStartY);
      t->SetBranchAddress("x_end"                 , &segmentEndY);
      t->SetBranchAddress("y_start"               , &segmentStartZ);
      t->SetBranchAddress("y_end"                 , &segmentEndZ);

      // loop the tree and find all of the wires associated with this channel
      bool isCorrectPlane = false;
      for (int ientry = 0; ientry < treeVec.at(itree)->GetEntries(); ientry++){
        t->GetEntry(ientry);
        if ( segmentSideAChannel == channelNumber || segmentSideBChannel == channelNumber){
          channelAssociatedSegmentsStartY->push_back(segmentStartY);
          channelAssociatedSegmentsEndY  ->push_back(segmentEndY);
          channelAssociatedSegmentsStartZ->push_back(segmentStartZ);
          channelAssociatedSegmentsEndZ  ->push_back(segmentEndZ);
          isCorrectPlane = true;
        }
      }
      if (isCorrectPlane == true) break;
 
    }

    geometryTree->Fill();
  }


  MF_LOG_DEBUG("TensionAnalysis")
    << "-- end TensionAnalysis::beginJob";


}

pdsp::APAAndSide MakeGeometryTree::getAPAInfoFromAnaAPANumber(int anaAPANumber)
{

  pdsp::APAAndSide thisAPAAndSide;

  if (anaAPANumber == 1 || anaAPANumber == 0){
    thisAPAAndSide.APABuildNumber         =  5;
    thisAPAAndSide.APAInstallationNumber  =  3;
    thisAPAAndSide.APAAnalysisNumber      =  1;
    thisAPAAndSide.sideToUse              =  pdsp::kB;
  }
  else if (anaAPANumber == 2 || anaAPANumber == 3){
    thisAPAAndSide.APABuildNumber         =  6;
    thisAPAAndSide.APAInstallationNumber  =  5;
    thisAPAAndSide.APAAnalysisNumber      =  2;
    thisAPAAndSide.sideToUse              =  pdsp::kA;
  }
  else if (anaAPANumber == 5 || anaAPANumber == 4){
    thisAPAAndSide.APABuildNumber         =  2;
    thisAPAAndSide.APAInstallationNumber  =  2;
    thisAPAAndSide.APAAnalysisNumber      =  5;
    thisAPAAndSide.sideToUse              =  pdsp::kB;
  }
  else if (anaAPANumber == 6 || anaAPANumber == 7){
    thisAPAAndSide.APABuildNumber         =  4;
    thisAPAAndSide.APAInstallationNumber  =  6;
    thisAPAAndSide.APAAnalysisNumber      =  6;
    thisAPAAndSide.sideToUse              =  pdsp::kA;
  }
  else if (anaAPANumber == 9 || anaAPANumber == 8){
    thisAPAAndSide.APABuildNumber         =  1;
    thisAPAAndSide.APAInstallationNumber  =  1;
    thisAPAAndSide.APAAnalysisNumber      =  9;
    thisAPAAndSide.sideToUse              =  pdsp::kB;
  }
  else if (anaAPANumber == 10 || anaAPANumber == 11){
    thisAPAAndSide.APABuildNumber         =  3;
    thisAPAAndSide.APAInstallationNumber  =  4;
    thisAPAAndSide.APAInstallationNumber  =  10;
    thisAPAAndSide.sideToUse              =  pdsp::kA;
  }
  else {
    thisAPAAndSide.APABuildNumber         =  -1;
    thisAPAAndSide.APAInstallationNumber  =  -1;
    thisAPAAndSide.APAAnalysisNumber      =  -1;
    thisAPAAndSide.sideToUse              =  pdsp::kUnknown;
  }

  return thisAPAAndSide;

}


DEFINE_ART_MODULE(MakeGeometryTree)
