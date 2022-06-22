////////////////////////////////////////////////////////////////////////
// Class:       GetGeometryInformation
// Plugin Type: analyzer (art v3_01_02)
// File:        GetGeometryInformation_module.cc
//
// Generated at Thu Jul 18 14:41:29 2019 by Adam Lister using cetskelgen
// from cetlib version v3_05_01.
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

// larsoft
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/RawDigit.h"

// root
#include "TTree.h"

class GetGeometryInformation;


class GetGeometryInformation : public art::EDAnalyzer {
public:
  explicit GetGeometryInformation(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  GetGeometryInformation(GetGeometryInformation const&) = delete;
  GetGeometryInformation(GetGeometryInformation&&) = delete;
  GetGeometryInformation& operator=(GetGeometryInformation const&) = delete;
  GetGeometryInformation& operator=(GetGeometryInformation&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  // Declare member data here.
  art::ServiceHandle< art::TFileService > tfs;
  geo::GeometryCore const* geom = lar::providerFrom<geo::Geometry>();

  TTree* geomTree;
  int channelNumber;
  int channelPlane;
  int numberOfAssociatedWires;
  std::vector< int >* tpcNumber = nullptr;
  std::vector< double >* wireStartX = nullptr;
  std::vector< double >* wireEndX   = nullptr;
  std::vector< double >* wireStartY = nullptr;
  std::vector< double >* wireEndY   = nullptr;
  std::vector< double >* wireStartZ = nullptr;
  std::vector< double >* wireEndZ   = nullptr;

};


GetGeometryInformation::GetGeometryInformation(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void GetGeometryInformation::analyze(art::Event const& e)
{
  // Implementation of required member function here.

 
  // get raw::RawDigit handle
  art::Handle< std::vector< raw::RawDigit > > rdHandle;
  e.getByLabel("tpcrawdecoder:daq:DecoderandReco", rdHandle);
  std::vector< art::Ptr< raw::RawDigit > > rdPtrVector;
  art::fill_ptr_vector(rdPtrVector, rdHandle);

  // loop RawDigits in the event, and get the channel from each of them
  for (size_t i_rd = 0; i_rd < rdPtrVector.size(); i_rd++){
    // clear vectors
    tpcNumber->resize(0);
    wireStartX->resize(0);
    wireEndX->resize(0);
    wireStartY->resize(0);
    wireEndY->resize(0);
    wireStartZ->resize(0);
    wireEndZ->resize(0);
 
    art::Ptr< raw::RawDigit > thisRawDigit = rdPtrVector.at(i_rd);

    channelNumber = thisRawDigit->Channel();
    channelPlane  = geom->View(thisRawDigit->Channel());

    std::vector<geo::WireID> theseWires = geom->ChannelToWire(channelNumber);
    numberOfAssociatedWires = theseWires.size();
    
    // for every wire associated to the channel then get the wire geometry
    for (int i_wire = 0; i_wire < numberOfAssociatedWires; i_wire++){
      geo::WireID thisWire = theseWires.at(i_wire);
      geo::WireGeo const& thisWireGeo = geom->Wire(thisWire);

      double wireStart[3] = {0.,0.,0.};
      double wireEnd[3] = {0.,0.,0.};

      thisWireGeo.GetStart(wireStart);
      thisWireGeo.GetEnd(wireEnd);

      tpcNumber->push_back(thisWire.asTPCID().TPC);
      wireStartX->push_back(wireStart[0]);
      wireEndX->push_back(wireEnd[0]);
      wireStartY->push_back(wireStart[1]);
      wireEndY->push_back(wireEnd[1]);
      wireStartZ->push_back(wireStart[2]);
      wireEndZ->push_back(wireEnd[2]);

    }

    geomTree->Fill();

  }

}

void GetGeometryInformation::beginJob()
{
  // Implementation of optional member function here.
  geomTree = tfs->make<TTree>("geometry_tree", "Geometry Tree");
  geomTree->Branch("channelNumber", &channelNumber);
  geomTree->Branch("channelPlane" , &channelPlane);
  geomTree->Branch("numberOfAssociatedWires", &numberOfAssociatedWires);
  geomTree->Branch("tpcNumber" , "std::vector<int>", &tpcNumber);
  geomTree->Branch("wireStartX", "std::vector<double>", &wireStartX);
  geomTree->Branch("wireEndX"  , "std::vector<double>", &wireEndX);
  geomTree->Branch("wireStartY", "std::vector<double>", &wireStartY);
  geomTree->Branch("wireEndY"  , "std::vector<double>", &wireEndY);
  geomTree->Branch("wireStartZ", "std::vector<double>", &wireStartZ);
  geomTree->Branch("wireEndZ"  , "std::vector<double>", &wireEndZ);


}

DEFINE_ART_MODULE(GetGeometryInformation)
