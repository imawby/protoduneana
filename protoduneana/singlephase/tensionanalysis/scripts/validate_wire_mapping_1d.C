/**
 * \brief Makes plots showing number of wires connected to each channel, broken up by APA
 * 
 * Usage: root -l -b <input_root_file> make_tension_var_plots_per_apa.C
 *
 * \author Adam Lister
 *
 * \date 2020-01-20
 *
 * Contact: adam.lister@wisc.edu
 **/

void style_plots(TH1D* U, TH1D* V, TH1D* X, int c){

  gStyle->SetHatchesSpacing(0.5);

  U->SetFillColor(c);
  V->SetFillColor(c);
  X->SetFillColor(c);

  U->SetLineColor(c);
  V->SetLineColor(c);
  X->SetLineColor(c);

  U->SetFillStyle(3354);
  V->SetFillStyle(3345);
  X->SetFillStyle(3395);

}

void validate_wire_mapping_1d(){

  std::vector<std::string> side = {
    "side_a",
    "side_b"
  };

  std::vector<std::string> sideTitle = {
    "Side A",
    "Side B"
  };

  TFile *fIn = new TFile("../data/tension_measurements_mod.root", "read");

  TTree* tUK001_U = (TTree*)fIn->Get("Tension_ProtoDUNE_UK001_ULAYER");
  TTree* tUK001_V = (TTree*)fIn->Get("Tension_ProtoDUNE_UK001_VLAYER");
  TTree* tUK001_X = (TTree*)fIn->Get("Tension_ProtoDUNE_UK001_XLAYER");
  TTree* tUK002_U = (TTree*)fIn->Get("Tension_ProtoDUNE_UK002_ULAYER");
  TTree* tUK002_V = (TTree*)fIn->Get("Tension_ProtoDUNE_UK002_VLAYER");
  TTree* tUK002_X = (TTree*)fIn->Get("Tension_ProtoDUNE_UK002_XLAYER");
  TTree* tUS001_U = (TTree*)fIn->Get("Tension_ProtoDUNE_US001_ULAYER");
  TTree* tUS001_V = (TTree*)fIn->Get("Tension_ProtoDUNE_US001_VLAYER");
  TTree* tUS001_X = (TTree*)fIn->Get("Tension_ProtoDUNE_US001_XLAYER");
  TTree* tUS002_U = (TTree*)fIn->Get("Tension_ProtoDUNE_US002_ULAYER");
  TTree* tUS002_V = (TTree*)fIn->Get("Tension_ProtoDUNE_US002_VLAYER");
  TTree* tUS002_X = (TTree*)fIn->Get("Tension_ProtoDUNE_US002_XLAYER");
  TTree* tUS003_U = (TTree*)fIn->Get("Tension_ProtoDUNE_US003_ULAYER");
  TTree* tUS003_V = (TTree*)fIn->Get("Tension_ProtoDUNE_US003_VLAYER");
  TTree* tUS003_X = (TTree*)fIn->Get("Tension_ProtoDUNE_US003_XLAYER");
  TTree* tUS004_U = (TTree*)fIn->Get("Tension_ProtoDUNE_US004_ULAYER");
  TTree* tUS004_V = (TTree*)fIn->Get("Tension_ProtoDUNE_US004_VLAYER");
  TTree* tUS004_X = (TTree*)fIn->Get("Tension_ProtoDUNE_US004_XLAYER");

  for (int i = 0; i < side.size(); i++){

    TH1D* hUK001_U = new TH1D("hUK001_U", (sideTitle.at(i)+std::string(";Channel Number;N Segments")).c_str(), 15360, 0, 15360);
    TH1D* hUK001_V = new TH1D("hUK001_V", (sideTitle.at(i)+std::string(";Channel Number;N Segments")).c_str(), 15360, 0, 15360);
    TH1D* hUK001_X = new TH1D("hUK001_X", (sideTitle.at(i)+std::string(";Channel Number;N Segments")).c_str(), 15360, 0, 15360);
    TH1D* hUK002_U = new TH1D("hUK002_U", (sideTitle.at(i)+std::string(";Channel Number;N Segments")).c_str(), 15360, 0, 15360);
    TH1D* hUK002_V = new TH1D("hUK002_V", (sideTitle.at(i)+std::string(";Channel Number;N Segments")).c_str(), 15360, 0, 15360);
    TH1D* hUK002_X = new TH1D("hUK002_X", (sideTitle.at(i)+std::string(";Channel Number;N Segments")).c_str(), 15360, 0, 15360);
    TH1D* hUS001_U = new TH1D("hUS001_U", (sideTitle.at(i)+std::string(";Channel Number;N Segments")).c_str(), 15360, 0, 15360);
    TH1D* hUS001_V = new TH1D("hUS001_V", (sideTitle.at(i)+std::string(";Channel Number;N Segments")).c_str(), 15360, 0, 15360);
    TH1D* hUS001_X = new TH1D("hUS001_X", (sideTitle.at(i)+std::string(";Channel Number;N Segments")).c_str(), 15360, 0, 15360);
    TH1D* hUS002_U = new TH1D("hUS002_U", (sideTitle.at(i)+std::string(";Channel Number;N Segments")).c_str(), 15360, 0, 15360);
    TH1D* hUS002_V = new TH1D("hUS002_V", (sideTitle.at(i)+std::string(";Channel Number;N Segments")).c_str(), 15360, 0, 15360);
    TH1D* hUS002_X = new TH1D("hUS002_X", (sideTitle.at(i)+std::string(";Channel Number;N Segments")).c_str(), 15360, 0, 15360);
    TH1D* hUS003_U = new TH1D("hUS003_U", (sideTitle.at(i)+std::string(";Channel Number;N Segments")).c_str(), 15360, 0, 15360);
    TH1D* hUS003_V = new TH1D("hUS003_V", (sideTitle.at(i)+std::string(";Channel Number;N Segments")).c_str(), 15360, 0, 15360);
    TH1D* hUS003_X = new TH1D("hUS003_X", (sideTitle.at(i)+std::string(";Channel Number;N Segments")).c_str(), 15360, 0, 15360);
    TH1D* hUS004_U = new TH1D("hUS004_U", (sideTitle.at(i)+std::string(";Channel Number;N Segments")).c_str(), 15360, 0, 15360);
    TH1D* hUS004_V = new TH1D("hUS004_V", (sideTitle.at(i)+std::string(";Channel Number;N Segments")).c_str(), 15360, 0, 15360);
    TH1D* hUS004_X = new TH1D("hUS004_X", (sideTitle.at(i)+std::string(";Channel Number;N Segments")).c_str(), 15360, 0, 15360);


    tUK001_U->Draw((side.at(i)+std::string("_channel_number >> hUK001_U")).c_str());
    tUK001_V->Draw((side.at(i)+std::string("_channel_number >> hUK001_V")).c_str());
    tUK001_X->Draw((side.at(i)+std::string("_channel_number >> hUK001_X")).c_str());
    tUK002_U->Draw((side.at(i)+std::string("_channel_number >> hUK002_U")).c_str());
    tUK002_V->Draw((side.at(i)+std::string("_channel_number >> hUK002_V")).c_str());
    tUK002_X->Draw((side.at(i)+std::string("_channel_number >> hUK002_X")).c_str());
    tUS001_U->Draw((side.at(i)+std::string("_channel_number >> hUS001_U")).c_str());
    tUS001_V->Draw((side.at(i)+std::string("_channel_number >> hUS001_V")).c_str());
    tUS001_X->Draw((side.at(i)+std::string("_channel_number >> hUS001_X")).c_str());
    tUS002_U->Draw((side.at(i)+std::string("_channel_number >> hUS002_U")).c_str());
    tUS002_V->Draw((side.at(i)+std::string("_channel_number >> hUS002_V")).c_str());
    tUS002_X->Draw((side.at(i)+std::string("_channel_number >> hUS002_X")).c_str());
    tUS003_U->Draw((side.at(i)+std::string("_channel_number >> hUS003_U")).c_str());
    tUS003_V->Draw((side.at(i)+std::string("_channel_number >> hUS003_V")).c_str());
    tUS003_X->Draw((side.at(i)+std::string("_channel_number >> hUS003_X")).c_str());
    tUS004_U->Draw((side.at(i)+std::string("_channel_number >> hUS004_U")).c_str());
    tUS004_V->Draw((side.at(i)+std::string("_channel_number >> hUS004_V")).c_str());
    tUS004_X->Draw((side.at(i)+std::string("_channel_number >> hUS004_X")).c_str());

    style_plots(hUK001_U, hUK001_V, hUK001_X, kBlack+0);
    style_plots(hUK002_U, hUK002_V, hUK002_X, kRed+0);
    style_plots(hUS001_U, hUS001_V, hUS001_X, kAzure+1);
    style_plots(hUS002_U, hUS002_V, hUS002_X, kOrange+1);
    style_plots(hUS003_U, hUS003_V, hUS003_X, kViolet+1);
    style_plots(hUS004_U, hUS004_V, hUS004_X, kGreen+1);

    
    //  hUK001_U->SetFillColor(kBlack);
    //  hUK001_U->SetFillStyle(3345);
    //  hUK001_V->SetFillStyle(3354);
    //  hUK002_U->SetLineColor(kRed);
    //  hUK002_V->SetLineColor(kRed);
    //  hUK002_X->SetLineColor(kRed);
    //  hUS001_U->SetLineColor(kAzure+1);
    //  hUS001_V->SetLineColor(kAzure+1);
    //  hUS001_X->SetLineColor(kAzure+1);
    //  hUS002_U->SetLineColor(kGreen+1);
    //  hUS002_V->SetLineColor(kGreen+1);
    //  hUS002_X->SetLineColor(kGreen+1);
    //  hUS003_U->SetLineColor(kOrange-1);
    //  hUS003_V->SetLineColor(kOrange-1);
    //  hUS003_X->SetLineColor(kOrange-1);
    //  hUS004_U->SetLineColor(kViolet+1);
    //  hUS004_V->SetLineColor(kViolet+1);
    //  hUS004_X->SetLineColor(kViolet+1);
    
    hUK001_U->GetYaxis()->SetRangeUser(0,3);

    TCanvas *c1 = new TCanvas("c1", "c1");
    c1->cd();

    hUK001_U->Draw();
    hUK001_V->Draw("same");
    hUK001_X->Draw("same");
    hUK002_U->Draw("same");
    hUK002_V->Draw("same");
    hUK002_X->Draw("same");
    hUS001_U->Draw("same");
    hUS001_V->Draw("same");
    hUS001_X->Draw("same");
    hUS002_U->Draw("same");
    hUS002_V->Draw("same");
    hUS002_X->Draw("same");
    hUS003_U->Draw("same");
    hUS003_V->Draw("same");
    hUS003_X->Draw("same");
    hUS004_U->Draw("same");
    hUS004_V->Draw("same");
    hUS004_X->Draw("same");

    TLegend* left = new TLegend(0.2, 0.7, 0.5, 0.85);
    left->AddEntry(hUK001_U, "UK001 U/V/Y");
    left->AddEntry(hUK002_U, "UK002 U/V/Y");
    left->AddEntry(hUS001_U, "US001 U/V/Y");
    TLegend* right = new TLegend(0.5, 0.7, 0.8, 0.85);
    right->AddEntry(hUS002_U, "US002 U/V/Y");
    right->AddEntry(hUS003_U, "US003 U/V/Y");
    right->AddEntry(hUS004_U, "US004 U/V/Y");

    left->SetBorderSize(0);
    right->SetBorderSize(0);

    left->Draw();
    right->Draw();

    c1->SaveAs((side.at(i)+std::string("_sanity_plots.png")).c_str());
    c1->SaveAs((side.at(i)+std::string("_sanity_plots.pdf")).c_str());

  }
}
