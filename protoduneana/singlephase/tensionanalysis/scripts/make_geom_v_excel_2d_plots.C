void make_geom_v_excel_2d_plots(){

  std::vector< std::string > names =
  {
    //"trackHitWireGeomYStart-trackHitWireSegmentYStart:trackHitWireGeomLength*10-trackHitWireSegmentLength",
    //"trackHitWireGeomYStart:trackHitWireGeomLength*10-trackHitWireSegmentLength",
    //"trackHitWireGeomYStart:trackHitWireGeomYStart-trackHitWireSegmentYStart",
    "trackHitWireNo:trackHitWireGeomLength*10-trackHitWireSegmentLength",
    "trackHitWireNo:trackHitWireGeomYStart-trackHitWireSegmentYStart"
  };

  std::vector< std::string > labels =
  {
    //";Geom - Excel Length (mm);Geom - Excel Y Start (mm)",
    //";Geom - Excel Length (mm); Geom Y Start",
    //";Geom - Excel Y Start (mm); Geom Y Start",
    ";Geom - Excel Length (mm); Excel Wire Number",
    ";Geom - Excel Y Start (mm); Excel Wire Number"
  };

  std::vector< std::vector< int > > bins = 
  {
    //{500, -20, 20, 1000, -50, 50},
    //{500, -20, 20, 1000, 0, 6000},
    //{500, -20, 20, 1000, 0, 6000},
    {500, -20, 20, 1151, 0, 1151},
    {500, -20, 20, 1151, 0, 1151}
  };

  TTree* tree = (TTree*)_file0->Get("tensionanalysis/analysis_tree");

  for (int i  = 0; i < names.size(); i++){

    TH2D* XPlane = new TH2D("XPlane", labels.at(i).c_str(), bins.at(i).at(0), bins.at(i).at(1), bins.at(i).at(2), bins.at(i).at(3), bins.at(i).at(4), bins.at(i).at(5));
    TH2D* UPlane = new TH2D("UPlane", labels.at(i).c_str(), bins.at(i).at(0), bins.at(i).at(1), bins.at(i).at(2), bins.at(i).at(3), bins.at(i).at(4), bins.at(i).at(5));
    TH2D* VPlane = new TH2D("VPlane", labels.at(i).c_str(), bins.at(i).at(0), bins.at(i).at(1), bins.at(i).at(2), bins.at(i).at(3), bins.at(i).at(4), bins.at(i).at(5));

    tree->Draw((names.at(i)+std::string(">> XPlane")).c_str(), "trackHitView==2 && trackHitWireSegmentYStart != -1");
    tree->Draw((names.at(i)+std::string(">> UPlane")).c_str(), "trackHitView==0 && trackHitWireSegmentYStart != -1");
    tree->Draw((names.at(i)+std::string(">> VPlane")).c_str(), "trackHitView==1 && trackHitWireSegmentYStart != -1");

    TCanvas *c1 = new TCanvas();
    XPlane->Draw("colz");
    c1->SaveAs((names.at(i)+"XPlane.png").c_str());

    UPlane->Draw("colz");
    c1->SaveAs((names.at(i)+"UPlane.png").c_str());

    VPlane->Draw("colz");
    c1->SaveAs((names.at(i)+"VPlane.png").c_str());

    XPlane->Delete();
    UPlane->Delete();
    VPlane->Delete();
  }

}
