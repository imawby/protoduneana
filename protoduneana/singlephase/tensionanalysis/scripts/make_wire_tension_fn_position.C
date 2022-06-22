/**
 * \brief Makes plots showing tension of all wired
 * 
 * this currently only produces plots for the U adn V planes. That should be 
 * changed at some point.
 *
 * Usage: root -l -b <input_root_file> validate_wire_mapping.C
 *
 * \author Adam Lister
 *
 * \date 2020-01-20
 *
 * Contact: adam.lister@wisc.edu
 **/

void make_wire_tension_fn_position(){

  /**
   * first make gradient, this is useful for plotting later
   **/
  gStyle->SetPalette(kRainBow);
  TColor::InvertPalette();
  TH2D* hGradient = new TH2D("hGradient", ";Wire Tension (N);", 1000,3,8,1,0,1);
  for (int i = 1; i < hGradient->GetNbinsX()+1; i++){
    hGradient->SetBinContent(i,1,hGradient->GetXaxis()->GetBinCenter(i));
  }
  hGradient->SetContour(1000);

  TCanvas *cGrad = new TCanvas("cGrad", "cGrad", 1000, 200);
  cGrad->SetBottomMargin(0.3);
  cGrad->SetLeftMargin(0.05);
  cGrad->SetRightMargin(0.05);
  hGradient->GetYaxis()->SetTitleOffset(10);
  hGradient->GetYaxis()->SetLabelOffset(10);
  hGradient->GetXaxis()->SetTitleSize(0.15);
  hGradient->GetXaxis()->SetLabelSize(0.15);
  hGradient->Draw("colz");
  cGrad->Modified();
  cGrad->Update();

  TPaletteAxis *palette = (TPaletteAxis*)hGradient->GetListOfFunctions()->FindObject("palette");
  palette->SetX1(2.0);
  palette->SetX2(2.0);

  cGrad->Modified();
  cGrad->Update();

  /**
   * now just need to define some vectors and loop over them
   **/
  std::vector< std::string > APA = {
    "UK001",
    "UK002",
    "US001",
    "US002",
    "US003",
    "US004"
  };

  std::vector< std::string > layers = {
    "ULAYER",
    "VLAYER",
    "XLAYER"
  };

  std::vector< std::string > side = {
    "a",
    "b"
  };

  /**
   * pull file from the data directory
   * TODO: pull information from tension analysis geometry_tree so we can 
   * produce results for the X plane
   **/
  TFile* _file0 = new TFile("../data/tension_measurements_mod.root", "read");

  for (int iapa = 0; iapa < APA.size(); iapa++){

    for (int ilayer = 0; ilayer < layers.size(); ilayer++){

      std::string treeName = "Tension_ProtoDUNE_" + APA.at(iapa) + "_" + layers.at(ilayer);

      TTree* t = (TTree*)_file0->Get(treeName.c_str());

      for (int iside = 0; iside < side.size(); iside++){

        float tension;
        float x_start = 0;
        float x_end = 5987.6;
        float y_start = 0;
        float y_end = 0;

        std::string tensionString = "side_" + side.at(iside) + "_final_tension";

        t->SetBranchAddress(tensionString.c_str() , &tension);
        if (layers.at(ilayer) != "XLAYER"){
          t->SetBranchAddress("x_start"             , &x_start);
          t->SetBranchAddress("x_end"               , &x_end);
          t->SetBranchAddress("y_start"             , &y_start);
          t->SetBranchAddress("y_end"               , &y_end);
        }

        TCanvas* c1 = new TCanvas("c1", "c1", 500, 1000);

        std::string titleString =  APA.at(iapa) + " " + layers.at(ilayer) + ", side " + side.at(iside) + ";z (mm);y (mm)";

        TH2D* h_background = new TH2D("h_background", titleString.c_str(), 100, -500, 2800, 100, -500, 6487);
        h_background->GetYaxis()->SetTitleOffset(1.6);
        h_background->Draw();

        TLine* bottom = new TLine(0    , 0    , 2300 , 0);
        TLine* top    = new TLine(0    , 5988 , 2300 , 5988);
        TLine* left   = new TLine(0    , 0    , 0    , 5988);
        TLine* right  = new TLine(2300 , 0    , 2300 , 5988);

        for (int i = 0; i < t->GetEntries(); i++){
          t->GetEntry(i);
          if (tension <= 0) continue;

          if (layers.at(ilayer) == "XLAYER"){
            y_start += 4.75;
            y_end   += 4.75;
          }

      
          /// pull colour from the gradient histogram
          Int_t ci = palette->GetBinColor(hGradient->GetXaxis()->FindBin(tension),1);
          TLine *line = new TLine(y_start, x_start, y_end, x_end);
          line->SetLineColor(ci);
          line->Draw("same");

        }
        bottom->Draw();
        top->Draw();
        left->Draw();
        right->Draw();

        std::string printString =  "low_tension_wires_"+APA.at(iapa) + "_" + layers.at(ilayer) + "_side_" + side.at(iside);
        c1->SaveAs((printString+".png").c_str());
        c1->SaveAs((printString+".pdf").c_str());

        h_background->Delete();

      }
    }
  }

  cGrad->cd();
  hGradient->Draw("col");
  hGradient->GetYaxis()->SetNdivisions(505);
  hGradient->GetXaxis()->CenterTitle();
  cGrad->SaveAs("gradient.png");
  cGrad->SaveAs("gradient.pdf");

}
