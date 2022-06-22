void make_geom_v_excel_plots(){

  std::vector< std::string > names =
  {
    "trackHitWireGeomYStart-trackHitWireSegmentYStart",
    "trackHitWireGeomZStart-trackHitWireSegmentZStart",
    "trackHitWireGeomLength*10-trackHitWireSegmentLength",
  };

  std::vector< std::string > labels =
  {
    ";Geom - Excel Y Start (mm);",
    ";Geom - Excel Z Start (mm);",
    ";Geom - Excel Length  (mm);",
  };

  std::vector< std::vector< int > > bins = 
  {
    {1000, -50, 50},
    {1000, -50, 50},
    {1000, -50, 50},
  };

  TTree* tree = (TTree*)_file0->Get("tensionanalysis/analysis_tree");

  for (int i  = 0; i < names.size(); i++){

    TH1D* XPlane = new TH1D("XPlane", labels.at(i).c_str(), bins.at(i).at(0), bins.at(i).at(1), bins.at(i).at(2));
    TH1D* UPlane = new TH1D("UPlane", labels.at(i).c_str(), bins.at(i).at(0), bins.at(i).at(1), bins.at(i).at(2));
    TH1D* VPlane = new TH1D("VPlane", labels.at(i).c_str(), bins.at(i).at(0), bins.at(i).at(1), bins.at(i).at(2));

    tree->Draw((names.at(i)+std::string(">> XPlane")).c_str(), "trackHitView==2");
    tree->Draw((names.at(i)+std::string(">> UPlane")).c_str(), "trackHitView==0");
    tree->Draw((names.at(i)+std::string(">> VPlane")).c_str(), "trackHitView==1");

    TCanvas *c1 = new TCanvas();
    XPlane->SetLineColor(kBlack);
    UPlane->SetLineColor(kGreen+1);
    VPlane->SetLineColor(kAzure+1);
    XPlane->DrawNormalized();
    UPlane->DrawNormalized("same");
    VPlane->DrawNormalized("same");
    c1->SaveAs((names.at(i)+".png").c_str());

    XPlane->Delete();
    UPlane->Delete();
    VPlane->Delete();
  }

}
