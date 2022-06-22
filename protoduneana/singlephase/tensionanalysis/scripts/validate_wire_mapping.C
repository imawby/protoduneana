/**
 * \brief Validates wire to channel mapping
 * 
 * This code uses the output geometry_tree from the tensionanalysis module 
 * to compare the wire information stored in the geometry and that which
 * was initially stored in excel sheets (and ported to root trees).
 *
 * In particular here, the number of wires connected to each channel is 
 * compared, in addition to the difference in wire lengths.
 * 
 * Note that this only compares the U and V plane wires since there is a
 * 1-to-1 correspondence for the X planes.
 *
 * Usage: root -l -b <input_root_file> validate_wire_mapping.C
 *
 * \author Adam Lister
 *
 * \date 2020-01-20
 *
 * Contact: adam.lister@wisc.edu
 **/

/**
 * struct used to store information about the wire
 * mostly useful for sorting
 **/
struct WireAnchors {
  double starty = -1;
  double endy   = -1;
  double startz = -1;
  double endz   = -1;
  double length = -1;
};

/**
 * used to calculate the length of the wire based on the start and 
 * end points
 **/
double calculate_length(double x1, double y1, double x2, double y2){
  double length = std::sqrt(std::pow(x2-x1,2)+std::pow(y2-y1,2));
  return length;
}

/**
 * main function
 **/
void validate_wire_mapping(){

  /// pass the file as described at top of file
  TTree* geom_tree = (TTree*)_file0->Get("tensionanalysis/geometry_tree");

  TH1D* hLengthDiffUTop = new TH1D("hLengthDiffUTop",
      ";Length Difference (mm);Number of Wires",
      100, -10, 20);

  TH1D* hLengthDiffUMid = new TH1D("hLengthDiffUMid",
      ";Length Difference (mm);Number of Wires",
      100, -10, 20);

  TH1D* hLengthDiffUBot = new TH1D("hLengthDiffUBot",
      ";Length Difference (mm);Number of Wires",
      100, -10, 20);

  TH1D* hLengthDiffVTop = new TH1D("hLengthDiffVTop",
      ";Length Difference (mm);Number of Wires",
      100, -10, 20);

  TH1D* hLengthDiffVMid = new TH1D("hLengthDiffVMid",
      ";Length Difference (mm);Number of Wires",
      100, -10, 20);

  TH1D* hLengthDiffVBot = new TH1D("hLengthDiffVBot",
      ";Length Difference (mm);Number of Wires",
      100, -10, 20);

  TH2D* hNumSegments    = new TH2D("hNumSegments",
      ";Number of Segments/Channel, Geometry;Number of Segments/Channel, Measured",
      4, 0, 4,
      4, 0, 4);

  int channelNumber;
  int channelAssociatedWiresPlane;
  std::vector< double >* channelAssociatedWiresStartY    = nullptr;
  std::vector< double >* channelAssociatedWiresEndY      = nullptr;
  std::vector< double >* channelAssociatedWiresStartZ    = nullptr;
  std::vector< double >* channelAssociatedWiresEndZ      = nullptr;
  std::vector< double >* channelAssociatedSegmentsStartY = nullptr;
  std::vector< double >* channelAssociatedSegmentsEndY   = nullptr;
  std::vector< double >* channelAssociatedSegmentsStartZ = nullptr;
  std::vector< double >* channelAssociatedSegmentsEndZ   = nullptr;

  geom_tree->SetBranchAddress("channelNumber"                   , &channelNumber                   );
  geom_tree->SetBranchAddress("channelAssociatedWiresPlane"     , &channelAssociatedWiresPlane     );
  geom_tree->SetBranchAddress("channelAssociatedWiresStartY"    , &channelAssociatedWiresStartY    );
  geom_tree->SetBranchAddress("channelAssociatedWiresEndY"      , &channelAssociatedWiresEndY      );
  geom_tree->SetBranchAddress("channelAssociatedWiresStartZ"    , &channelAssociatedWiresStartZ    );
  geom_tree->SetBranchAddress("channelAssociatedWiresEndZ"      , &channelAssociatedWiresEndZ      );
  geom_tree->SetBranchAddress("channelAssociatedSegmentsStartY" , &channelAssociatedSegmentsStartY );
  geom_tree->SetBranchAddress("channelAssociatedSegmentsEndY"   , &channelAssociatedSegmentsEndY   );
  geom_tree->SetBranchAddress("channelAssociatedSegmentsStartZ" , &channelAssociatedSegmentsStartZ );
  geom_tree->SetBranchAddress("channelAssociatedSegmentsEndZ"   , &channelAssociatedSegmentsEndZ   );

  for (int i = 0; i < geom_tree->GetEntries(); i++){
    geom_tree->GetEntry(i);

    int wiresSize    = channelAssociatedWiresStartY   ->size();
    int segmentsSize = channelAssociatedSegmentsStartY->size();

    hNumSegments->Fill(segmentsSize, wiresSize);

    if (wiresSize != segmentsSize){
      std::cout 
        << "channel " << channelNumber << " has " 
        << wiresSize << "wires simulated, but "
        << segmentsSize << "measured wires".
      throw std::logic_error("uh oh");
    }

    std::vector< WireAnchors > wireAnchorsMeasured;
    std::vector< WireAnchors > wireAnchorsSimulated;

    double maxExt = 0;
    if (channelAssociatedWiresPlane == 0)
      maxExt = 5991;
    if (channelAssociatedWiresPlane == 1)
      maxExt = 5987.6; 

    /**
     * fill vector of WireAnchors information for wires attached to 
     * this channel
     *
     * note that the direction of the measured channels must be flipped to
     * coincide with the simulated wires, and that the simulated wires
     * are stored in cm, not mm, so correct for that
     **/
    for (int iw = 0; iw < wiresSize; iw++){
      WireAnchors thisWireAnchorMeasured;
      thisWireAnchorMeasured.starty = maxExt-channelAssociatedSegmentsStartY->at(iw);
      thisWireAnchorMeasured.startz = maxExt-channelAssociatedSegmentsStartZ->at(iw);
      thisWireAnchorMeasured.endy   = maxExt-channelAssociatedSegmentsEndY  ->at(iw);
      thisWireAnchorMeasured.endz   = maxExt-channelAssociatedSegmentsEndZ  ->at(iw);
      thisWireAnchorMeasured.length = calculate_length(thisWireAnchorMeasured.starty,
          thisWireAnchorMeasured.startz,
          thisWireAnchorMeasured.endy,
          thisWireAnchorMeasured.endz);
      WireAnchors thisWireAnchorSimulated;
      thisWireAnchorSimulated.starty = channelAssociatedWiresStartY->at(iw)*10;
      thisWireAnchorSimulated.startz = channelAssociatedWiresStartZ->at(iw)*10;
      thisWireAnchorSimulated.endy   = channelAssociatedWiresEndY  ->at(iw)*10;
      thisWireAnchorSimulated.endz   = channelAssociatedWiresEndZ  ->at(iw)*10;
      thisWireAnchorSimulated.length = calculate_length(thisWireAnchorSimulated.starty,
          thisWireAnchorSimulated.startz,
          thisWireAnchorSimulated.endy,
          thisWireAnchorSimulated.endz);

      wireAnchorsMeasured .push_back(thisWireAnchorMeasured);
      wireAnchorsSimulated.push_back(thisWireAnchorSimulated);
    }

    /**
     * sort them based on their proximity to the top of the TPC
     **/
    std::sort(wireAnchorsMeasured.begin(), wireAnchorsMeasured.end(), 
        [](auto const &a, auto const &b) { return a.starty < b.starty; });

    std::sort(wireAnchorsSimulated.begin(), wireAnchorsSimulated.end(), 
        [](auto const &a, auto const &b) { return a.starty < b.starty; });

    /**
     * ... and plot the information
     **/
    for (int iw = 0; iw < wireAnchorsSimulated.size(); iw++){
      if (channelAssociatedWiresPlane == 0){
        if (iw == 0)
          hLengthDiffUTop->Fill(wireAnchorsSimulated.at(iw).length - wireAnchorsMeasured.at(iw).length);
        if (iw == 1)
          hLengthDiffUMid->Fill(wireAnchorsSimulated.at(iw).length - wireAnchorsMeasured.at(iw).length);
        if (iw == 2)
          hLengthDiffUBot->Fill(wireAnchorsSimulated.at(iw).length - wireAnchorsMeasured.at(iw).length);
      }
      if (channelAssociatedWiresPlane == 1){
        if (iw == 0)
          hLengthDiffVTop->Fill(wireAnchorsSimulated.at(iw).length - wireAnchorsMeasured.at(iw).length);
        if (iw == 1)
          hLengthDiffVMid->Fill(wireAnchorsSimulated.at(iw).length - wireAnchorsMeasured.at(iw).length);
        if (iw == 2)
          hLengthDiffVBot->Fill(wireAnchorsSimulated.at(iw).length - wireAnchorsMeasured.at(iw).length);
      }

    }
  }

  TCanvas *c1 = new TCanvas("c1", "c1", 600, 400);

  hLengthDiffUTop->SetLineColor(kRed+1);
  hLengthDiffUMid->SetLineColor(kOrange+1);
  hLengthDiffUBot->SetLineColor(kYellow-3);
  hLengthDiffUTop->SetMarkerColor(kRed+1);
  hLengthDiffUMid->SetMarkerColor(kOrange+1);
  hLengthDiffUBot->SetMarkerColor(kYellow-3);
  hLengthDiffUTop->SetFillColor(kRed+1);
  hLengthDiffUMid->SetFillColor(kOrange+1);
  hLengthDiffUBot->SetFillColor(kYellow-3);
  hLengthDiffUTop->SetFillStyle(3345);
  hLengthDiffUMid->SetFillStyle(3345);
  hLengthDiffUBot->SetFillStyle(3345);

  hLengthDiffVTop->SetLineColor(kGreen+1);
  hLengthDiffVMid->SetLineColor(kAzure+1);
  hLengthDiffVBot->SetLineColor(kViolet+1);
  hLengthDiffVTop->SetMarkerColor(kGreen+1);
  hLengthDiffVMid->SetMarkerColor(kAzure+1);
  hLengthDiffVBot->SetMarkerColor(kViolet+1);
  hLengthDiffVTop->SetFillColor(kGreen+1);
  hLengthDiffVMid->SetFillColor(kAzure+1);
  hLengthDiffVBot->SetFillColor(kViolet+1);
  hLengthDiffVTop->SetFillStyle(3354);
  hLengthDiffVMid->SetFillStyle(3354);
  hLengthDiffVBot->SetFillStyle(3354);

  hLengthDiffUTop->Draw();
  hLengthDiffUMid->Draw("same");
  hLengthDiffUBot->Draw("same");
  hLengthDiffVTop->Draw("same");
  hLengthDiffVMid->Draw("same");
  hLengthDiffVBot->Draw("same");

  TLegend* leg = new TLegend(0.65, 0.55, 0.85, 0.85);
  leg->AddEntry(hLengthDiffUTop, "U Plane, Top Seg.");
  leg->AddEntry(hLengthDiffUMid, "U Plane, Mid Seg.");
  leg->AddEntry(hLengthDiffUBot, "U Plane, Bot Seg.");
  leg->AddEntry(hLengthDiffVTop, "V Plane, Top Seg.");
  leg->AddEntry(hLengthDiffVMid, "V Plane, Mid Seg.");
  leg->AddEntry(hLengthDiffVBot, "V Plane, Bot Seg.");
  leg->SetFillStyle(0);
  leg->SetLineWidth(0);
  leg->Draw("same");

  c1->SaveAs("LengthDifferences.png");
  c1->SaveAs("LengthDifferences.pdf");

  TCanvas *c2 = new TCanvas("c2", "c2", 600, 500);
  c2->SetRightMargin(0.18);
  hNumSegments->Draw("colz");
  gPad->Modified();
  c2->SaveAs("NumberOfSegments.png");
  c2->SaveAs("NumberOfSegments.pdf");

}
