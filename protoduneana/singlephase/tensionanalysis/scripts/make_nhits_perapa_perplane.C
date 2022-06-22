void make_nhits_perapa_perplane(){

  TTree* tree = (TTree*)_file0->Get("tensionanalysis/analysis_tree");

  // define cuts

  std::vector< std::string > cutsTPC = {
    "trackHitAPABuildNumber>-1", // all tpcs 
    "trackHitAPABuildNumber==1", // tpc 1
    "trackHitAPABuildNumber==2", // tpc 2
    "trackHitAPABuildNumber==3", // tpc 3
    "trackHitAPABuildNumber==4", // tpc 4
    "trackHitAPABuildNumber==5", // tpc 5
    "trackHitAPABuildNumber==6"  // tpc 6
    };

  for (int itpc = 0; itpc < cutsTPC.size(); itpc++){

      std::string nameString  = cutsTPC.at(itpc);
      std::string uNameString = "uPlane_" + nameString;
      std::string vNameString = "vPlane_" + nameString;
      std::string xNameString = "xPlane_" + nameString;

      TH1D* uPlane = new TH1D(uNameString.c_str(), ";Wire Tension (N);Number of Hits", 11, 2.5, 8);
      TH1D* vPlane = new TH1D(vNameString.c_str(), ";Wire Tension (N);Number of Hits", 11, 2.5, 8);
      TH1D* xPlane = new TH1D(xNameString.c_str(), ";Wire Tension (N);Number of Hits", 11, 2.5, 8);

      tree->Draw(("trackHitWireTension >> "+uNameString).c_str(), ("trackHitView == 0 && " + cutsTPC.at(itpc)).c_str());
      tree->Draw(("trackHitWireTension >> "+vNameString).c_str(), ("trackHitView == 1 && " + cutsTPC.at(itpc)).c_str());
      tree->Draw(("trackHitWireTension >> "+xNameString).c_str(), ("trackHitView == 2 && " + cutsTPC.at(itpc)).c_str());

      uPlane->SetLineColor(kGreen+1);
      vPlane->SetLineColor(kAzure+1);
      xPlane->SetLineColor(kOrange+1);
      uPlane->SetLineWidth(2);
      vPlane->SetLineWidth(2);
      xPlane->SetLineWidth(2);
      xPlane->GetYaxis()->SetRangeUser(1, std::max({uPlane->GetMaximum(), vPlane->GetMaximum(), xPlane->GetMaximum()})*100);

      TCanvas *c1 = new TCanvas("c1", "c1", 600, 500);
      xPlane->Draw();
      uPlane->Draw("same");
      vPlane->Draw("same");
     
      TLegend* leg = new TLegend(0.3, 0.7, 0.5, 0.85);
      leg->AddEntry(uPlane, "U Plane");
      leg->AddEntry(vPlane, "V Plane");
      leg->AddEntry(xPlane, "x Plane");
      leg->SetLineWidth(0);
      leg->SetFillStyle(0);
      leg->Draw("same");

      c1->SetLogy();
      c1->SaveAs(("nhits"+nameString+".png").c_str());

  }

}
