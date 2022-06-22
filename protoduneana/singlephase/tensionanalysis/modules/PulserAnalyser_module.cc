////////////////////////////////////////////////////////////////////////
// Class:       PulserAnalyser
// Plugin Type: analyzer (art v3_01_02)
// File:        PulserAnalyser_module.cc
//
// Generated at Mon Jun 24 21:12:55 2019 by Adam Lister using cetskelgen
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
#include "lardataobj/RawData/RawDigit.h"

// root
#include "TTree.h"
#include "TH1.h"


namespace pdsp {
  class PulserAnalyser;
}


class pdsp::PulserAnalyser : public art::EDAnalyzer {
public:
  explicit PulserAnalyser(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PulserAnalyser(PulserAnalyser const&) = delete;
  PulserAnalyser(PulserAnalyser&&) = delete;
  PulserAnalyser& operator=(PulserAnalyser const&) = delete;
  PulserAnalyser& operator=(PulserAnalyser&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  // services
  art::ServiceHandle< art::TFileService > tfs;

  TH1D* waveform;

  TTree* anaTree;

  int run;
  int subRun;
  int event;
  std::vector<int>* deadChannels = nullptr;

  std::string fRawDigitLabel;
  float fADCLimit;
  bool fMakeSingleWirePlot;
  int fSelectedChannelNo;
  int eventCounter;

};


pdsp::PulserAnalyser::PulserAnalyser(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{

  fhicl::ParameterSet const pLabels = p.get< fhicl::ParameterSet >("ProducerLabels");
  fhicl::ParameterSet const pCuts   = p.get< fhicl::ParameterSet >("CutValues");

  fRawDigitLabel      = pLabels.get<std::string>("RawDigitLabel");
  fADCLimit           = pCuts.get<float>("ADCLimit");
  fMakeSingleWirePlot = pCuts.get<bool>("MakeSingleWirePlot");
  fSelectedChannelNo  = pCuts.get<int>("SelectedChannelNo");

}

void pdsp::PulserAnalyser::analyze(art::Event const& e)
{

  deadChannels->resize(0);

  run = e.run();
  subRun = e.subRun();
  event = e.event();

  art::Handle< std::vector< raw::RawDigit > >rawDigitHandle;
  e.getByLabel(fRawDigitLabel, rawDigitHandle);
  std::vector< art::Ptr< raw::RawDigit > > rawDigitPtrVector;
  art::fill_ptr_vector(rawDigitPtrVector, rawDigitHandle);

  for (size_t i_rd = 0; i_rd < rawDigitPtrVector.size(); i_rd++){

    art::Ptr< raw::RawDigit > thisRawDigit = rawDigitPtrVector.at(i_rd);

    int n_adc_gt_limit = 0;

    //loop over raw digit and look for channels
    for (size_t i_adc = 0; i_adc < thisRawDigit->NADC(); i_adc++){

      if (thisRawDigit->ADC(i_adc)-thisRawDigit->GetPedestal() > fADCLimit)
        n_adc_gt_limit++;

      if (fSelectedChannelNo == (int)thisRawDigit->Channel() && eventCounter == 1)
        waveform->SetBinContent(i_adc+1, thisRawDigit->ADC(i_adc)-thisRawDigit->GetPedestal());

    }

    if (n_adc_gt_limit == 0){
      std::cout << "-- number of dead channels: " << n_adc_gt_limit << std::endl; 
      deadChannels->push_back(thisRawDigit->Channel());
    }
  }

  anaTree->Fill();

}

void pdsp::PulserAnalyser::beginJob()
{

  anaTree = tfs->make<TTree>("analysis_tree", "analysis tree");
  anaTree->Branch("run", &run);
  anaTree->Branch("subRun", &subRun);
  anaTree->Branch("event", &event);
  anaTree->Branch("deadChannels", "std::vector<int>", &deadChannels);

  if (fMakeSingleWirePlot)
    waveform = tfs->make<TH1D>("SingleWirePlot", ";Time (ticks);Q (ADC)", 6000, 0, 6000);
  eventCounter = 1;

}

DEFINE_ART_MODULE(pdsp::PulserAnalyser)
