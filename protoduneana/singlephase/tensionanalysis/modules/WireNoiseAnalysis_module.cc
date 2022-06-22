////////////////////////////////////////////////////////////////////////
// Class:       WireNoiseAnalysis
// Plugin Type: analyzer (art v3_02_06)
// File:        WireNoiseAnalysis_module.cc
//
// Generated at Thu Feb 27 17:58:58 2020 by Adam Lister using cetskelgen
// from cetlib version v3_07_02.
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

// art
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

// larsoft
#include "lardataobj/RawData/RawDigit.h"
#include "larcore/Geometry/Geometry.h"

// root
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TString.h"

namespace pdsp {
  class WireNoiseAnalysis;
}


class pdsp::WireNoiseAnalysis : public art::EDAnalyzer {
public:
  explicit WireNoiseAnalysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  WireNoiseAnalysis(WireNoiseAnalysis const&) = delete;
  WireNoiseAnalysis(WireNoiseAnalysis&&) = delete;
  WireNoiseAnalysis& operator=(WireNoiseAnalysis const&) = delete;
  WireNoiseAnalysis& operator=(WireNoiseAnalysis&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;

  // function to resize vectors to 0
  void resizeVectors();
private:

  // services
  geo::GeometryCore const* thisGeom = lar::providerFrom<geo::Geometry>();
  art::ServiceHandle< art::TFileService > tfs;

  // variables
  bool isUseU = false;
  bool isUseV = false;
  bool isUseX = false;
  std::map<int, float> channelTensionMap;

  // fhicl parameters
  std::string        fRawDigitLabel;
  std::vector<int>   fPlanesToUse;

  // TTree
  TTree* anaTree;
  int run;
  int subRun;
  int event;
  std::vector<int>*      rawDigitChannel  = nullptr;
  std::vector<int>*      rawDigitView     = nullptr; 
  std::vector< double >* rawDigitPedestal = nullptr;
  std::vector< double >* rawDigitNoiseRMS = nullptr;
  std::vector< double >* rawDigitWireTension = nullptr;
};


pdsp::WireNoiseAnalysis::WireNoiseAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
{

  fhicl::ParameterSet pLabels = p.get<fhicl::ParameterSet>("ProducerLabels");
  fhicl::ParameterSet pCuts   = p.get<fhicl::ParameterSet>("CutValues");

  fRawDigitLabel = pLabels.get< std::string >       ("RawDigitLabel");
  fPlanesToUse   = pCuts  .get< std::vector<int> >  ("PlanesToUse",  {2});
}

void pdsp::WireNoiseAnalysis::analyze(art::Event const& e)
{

  run    = e.run();
  subRun = e.subRun();
  event  = e.event();

  this->resizeVectors();

  art::Handle< std::vector< raw::RawDigit > > rawDigitHandle;
  e.getByLabel(fRawDigitLabel, rawDigitHandle);
  std::vector< art::Ptr< raw::RawDigit > > rawDigitPtrVector;
  art::fill_ptr_vector(rawDigitPtrVector, rawDigitHandle);

  for (auto& thisRawDigit : rawDigitPtrVector ){
  
    int thisRawDigitChannel = thisRawDigit->Channel();
    int thisRawDigitView    = thisGeom->View(thisRawDigitChannel);

    // only use raw digits from a given channel
    if (thisRawDigitView == 0 && !isUseU) continue;
    if (thisRawDigitView == 1 && !isUseV) continue;
    if (thisRawDigitView == 2 && !isUseX) continue;

    // get pedestal level
    float thisRawDigitPedestal = thisRawDigit->GetPedestal();

    MF_LOG_DEBUG("pdsp::WireNoiseAnalysis::beginJob")
      << "Pedestal for raw::RawDigit with channel "
      << rawDigitChannel
      << " is "
      << thisRawDigitPedestal;

    std::vector<short> rawDigitADCs = thisRawDigit->ADCs();
    short minADC = *std::min_element(rawDigitADCs.begin(), rawDigitADCs.end());
    short maxADC = *std::max_element(rawDigitADCs.begin(), rawDigitADCs.end());

    // fill 1 dimensional histogram with all of the ADC values
    // from the raw waveform
    TH1D* hADCVals = new TH1D("hADCVals", 
                              ";ADC Values;",
                              int(maxADC-minADC+1), 
                              minADC, 
                              maxADC);

    for (short thisADC : rawDigitADCs){
      hADCVals->Fill(thisADC);
    }

    // we can now calculate what the 1 sigma value is for this
    // distribution
    double thisRawDigitNoiseRMS = 0.0;
    double pars[3];
    double quantile = 0;
    if (hADCVals->GetEntries() > 0){
      quantile = 0.5-0.34;
      hADCVals->GetQuantiles(1, &pars[0], &quantile);

      quantile = 0.5;
      hADCVals->GetQuantiles(1, &pars[1], &quantile);

      quantile = 0.5+0.34;
      hADCVals->GetQuantiles(1, &pars[2], &quantile);

      thisRawDigitNoiseRMS = sqrt((pow(pars[1]-pars[0],2)+pow(pars[2]-pars[1],2))/2.);
    }

    rawDigitChannel ->push_back(thisRawDigitChannel);
    rawDigitView    ->push_back(thisRawDigitView);
    rawDigitPedestal->push_back(thisRawDigitPedestal);
    rawDigitNoiseRMS->push_back(thisRawDigitNoiseRMS);
    rawDigitWireTension->push_back(channelTensionMap.find(thisRawDigitChannel)->second);

  }

  anaTree->Fill();

}

void pdsp::WireNoiseAnalysis::beginJob()
{

  // which planes are we using?
  for (int plane : fPlanesToUse){
    if (plane == 0)
      isUseU = true;
    if (plane == 1)
      isUseV = true;
    if (plane == 2)
      isUseX = true;
  }

  // setup tree
  anaTree = tfs->make<TTree>("analysis_tree" , "analysis tree");

  anaTree->Branch("run"              , &run);
  anaTree->Branch("subRun"           , &subRun);
  anaTree->Branch("event"            , &event);
  anaTree->Branch("rawDigitChannel"  , "std::vector<int>"    , &rawDigitChannel);
  anaTree->Branch("rawDigitView"     , "std::vector<int>"    , &rawDigitView);
  anaTree->Branch("rawDigitPedestal" , "std::vector<double>" , &rawDigitPedestal);
  anaTree->Branch("rawDigitNoiseRMS" , "std::vector<double>" , &rawDigitNoiseRMS);
  anaTree->Branch("rawDigitWireTension", "std::vector<double>", &rawDigitWireTension);

  // loop all Trees, make map of channel number to wire tension (for X plane wires)

  // read in ROOT trees containing wire information
  std::string tensionPath;
  cet::search_path sp("FW_SEARCH_PATH");
  if(!sp.find_file("tensionanalysis/data/tension_measurements_mod.root", tensionPath)){
    throw cet::exception("FileError")
      << "Cannot find tension_measurements.root file "
      << " bail ungracefully\n\n"
      << __FILE__ << ":" << __LINE__;
  }
  TFile* tensionsFile = new TFile(tensionPath.c_str(), "read");

  TTree* treeXLayerUS001 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US001_XLAYER");
  TTree* treeULayerUS001 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US001_ULAYER");
  TTree* treeVLayerUS001 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US001_VLAYER");
  TTree* treeXLayerUS002 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US002_XLAYER");
  TTree* treeULayerUS002 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US002_ULAYER");
  TTree* treeVLayerUS002 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US002_VLAYER");
  TTree* treeXLayerUS003 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US003_XLAYER");
  TTree* treeULayerUS003 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US003_ULAYER");
  TTree* treeVLayerUS003 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US003_VLAYER");
  TTree* treeXLayerUS004 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US004_XLAYER");
  TTree* treeULayerUS004 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US004_ULAYER");
  TTree* treeVLayerUS004 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_US004_VLAYER");
  TTree* treeXLayerUK001 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_UK001_XLAYER");
  TTree* treeULayerUK001 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_UK001_ULAYER");
  TTree* treeVLayerUK001 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_UK001_VLAYER");
  TTree* treeXLayerUK002 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_UK002_XLAYER");
  TTree* treeULayerUK002 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_UK002_ULAYER");
  TTree* treeVLayerUK002 = (TTree*)tensionsFile->Get("Tension_ProtoDUNE_UK002_VLAYER");

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

  int side_a_channel_number;
  int side_b_channel_number;
  float side_a_final_tension;
  float side_b_final_tension;

  for (auto tree : treeVec){
    if (TString(tree->GetName()).Contains("ULayer") ||
        TString(tree->GetName()).Contains("VLayer")) continue;

    tree->SetBranchAddress("side_a_channel_number" , &side_a_channel_number);
    tree->SetBranchAddress("side_b_channel_number" , &side_b_channel_number);
    tree->SetBranchAddress("side_a_final_tension" , &side_a_final_tension);
    tree->SetBranchAddress("side_b_final_tension" , &side_b_final_tension);

    for (int iEnt = 0; iEnt < tree->GetEntries(); iEnt++){
      tree->GetEntry(iEnt);

      channelTensionMap.insert(std::make_pair(side_a_channel_number, side_a_final_tension));
      channelTensionMap.insert(std::make_pair(side_b_channel_number, side_b_final_tension));

    }

  }


}

void pdsp::WireNoiseAnalysis::resizeVectors(){
  rawDigitChannel ->resize(0);
  rawDigitView    ->resize(0);
  rawDigitPedestal->resize(0);
  rawDigitChannel ->resize(0);
  rawDigitWireTension->resize(0);
}

DEFINE_ART_MODULE(pdsp::WireNoiseAnalysis)
