void place_fit_results(double grad){

  grad = std::round(grad*100)/100;

  std::string text = "Gradient: y = " + std::to_string(grad) + "x";

  TPaveText* pt = new TPaveText(0.2, 0.8, 0.5, 0.85, "NDC");
  pt->AddText(text.c_str());
  pt->SetFillColor(kWhite);

  pt->Draw("same");
}

void make_final_cern_tension_plots(){

  TFile* _file0 = new TFile("/dune/app/users/alister1/dunetpc_v08_35_00/srcs/protoduneana/protoduneana/singlephase/tensionanalysis/data/tension_measurements_mod.root");

  std::vector<std::string> apaLabels = 
  {
    "UK001",
    "UK002",
    "US001",
    "US002",
    "US003",
    "US004"
  };

  std::vector<std::string> planeLabels =
  {
    "ULAYER",
    "VLAYER",
    "XLAYER"
  };

  std::vector<TH2D*> histVecA;
  std::vector<TH2D*> histVecB;

  TF1* lin = new TF1("lin", "x", 0, 10);
  lin->SetLineColor(kGray);
  lin->SetLineStyle(1);
  lin->SetNpx(1000);

  TF1* linFit = new TF1("linFit", "[0]*x", 0, 10);
  linFit->SetParameter(0,1.0);
  linFit->SetLineColor(kRed);
  linFit->SetNpx(1000);

  for (int iapa = 0; iapa < apaLabels.size(); ++iapa){

    for (int iplane = 0; iplane < planeLabels.size(); ++iplane){

      std::string treeName    = "Tension_ProtoDUNE_" +
        apaLabels.at(iapa)   +
        "_"                  +
        planeLabels.at(iplane);

      TTree* tree = (TTree*)_file0->Get(treeName.c_str());

      std::string histName    = treeName+"_nomvcern";
      std::string histNameA   = histName+"_A";
      std::string histNameB   = histName+"_B";

      histVecA.push_back(new TH2D(histNameA.c_str(), 
            (histNameA+";US/UK Wire Tennsion (N);CERN Wire Tension (N)").c_str(), 
            50, 0, 10, 
            50, 0, 10));
      histVecB.push_back(new TH2D(histNameB.c_str(), 
            (histNameB+";US/UK Wire Tennsion (N);CERN Wire Tension (N)").c_str(), 
            50, 0, 10, 
            50, 0, 10));

      tree->Draw(("side_a_cern_tension:side_a_final_tension >> " + histNameA).c_str(),
          "side_a_cern_tension > 0");
      tree->Draw(("side_b_cern_tension:side_b_final_tension >> " + histNameB).c_str(),
          "side_b_cern_tension > 0");

      if (histVecA.back()->Integral() > 0 && !histVecA.back()->IsZombie()){
        TCanvas *c1a = new TCanvas("c1a", "c1a", 500, 500);
        c1a->SetRightMargin(0.15);
        c1a->cd();
        histVecA.back()->Draw("colz");
        lin->Draw("same");
        gPad->Modified();
        gPad->Update();
        TPaletteAxis *aPalette = (TPaletteAxis*)(histVecA.back()->GetListOfFunctions()->FindObject("palette"));
        aPalette->SetX1NDC(0.87);
        aPalette->SetX2NDC(0.92);
        gPad->Modified();
        histVecA.back()->Fit("linFit", "", "", 0, 10);
        place_fit_results(linFit->GetParameter(0));
        c1a->SaveAs(("side_a_"+histName+".png").c_str());
      }
      if (histVecB.back()->Integral() > 0 && !histVecB.back()->IsZombie()){
        TCanvas *c1b = new TCanvas("c1b", "c1b", 500, 500);
        c1b->SetRightMargin(0.15);
        c1b->cd();
        histVecB.back()->Draw("colz");
        lin->Draw("same");
        gPad->Modified();
        gPad->Update();
        TPaletteAxis *aPalette = (TPaletteAxis*)(histVecB.back()->GetListOfFunctions()->FindObject("palette"));
        aPalette->SetX1NDC(0.87);
        aPalette->SetX2NDC(0.92);
        gPad->Modified();
        histVecB.back()->Fit("linFit", "", "", 0, 10);
        place_fit_results(linFit->GetParameter(0));
        c1b->SaveAs(("side_b_"+histName+".png").c_str());
      }
    }
  }

  TH2D* hTotal = new TH2D("hTotal", ";US/UK Wire Tension (N);CERN Wire Tension (N)", 50, 0, 10, 50, 0, 10);
  for (int i = 0; i < histVecA.size(); i++){
    hTotal->Add(histVecA.at(i), 1);
  }
  for (int i = 0; i < histVecB.size(); i++){
    hTotal->Add(histVecB.at(i), 1);
  }

  TCanvas *c1 = new TCanvas("c1", "c1", 500, 500);
  c1->SetRightMargin(0.15);
  c1->cd();
  hTotal->Draw("colz");
  lin->Draw("same");
  gPad->Modified();
  gPad->Update();
  TPaletteAxis *palette = (TPaletteAxis*)hTotal->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.87);
  palette->SetX2NDC(0.92);
  gPad->Modified();
  hTotal->Fit("linFit", "", "", 0, 10);
  place_fit_results(linFit->GetParameter(0));
  c1->SaveAs("total.png");

}
