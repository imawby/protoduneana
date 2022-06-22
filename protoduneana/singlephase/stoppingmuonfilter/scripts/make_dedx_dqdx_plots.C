void plot(TH2D* hist, std::string name){
  TCanvas* c1 = new TCanvas("c1", "c1", 700, 500);
  c1->cd();
  c1->SetRightMargin(0.15);
  hist->Draw("colz");
  gPad->Modified();
  gPad->Update();
  TPaletteAxis *aPalette = (TPaletteAxis*)(hist->GetListOfFunctions()->FindObject("palette"));
  gPad->Modified(); gPad->Update();
  aPalette->SetX1NDC(0.85);
  aPalette->SetX2NDC(0.90);
  gPad->Modified();
  c1->SaveAs((name+".png").c_str());
}

void make_dedx_dqdx_plots(){

  //TFile* _file0 = new TFile("/dune/data/users/alister1/2020-Jan-CM/stoppingmuonfilter.root", "read");

  TTree* tree = (TTree*)_file0->Get("stopmufilter/analysis_tree");

  TH2D* uplanedqdxresrg = new TH2D("uplanedqdxresrg", "U Plane;Residual Range (cm); dQ/dx", 300, 0, 300, 100, 0, 1500);
  TH2D* vplanedqdxresrg = new TH2D("vplanedqdxresrg", "V Plane;Residual Range (cm); dQ/dx", 300, 0, 300, 100, 0, 1500);
  TH2D* xplanedqdxresrg = new TH2D("xplanedqdxresrg", "X Plane;Residual Range (cm); dQ/dx", 300, 0, 300, 100, 0, 1500);
  TH2D* uplanededxresrg = new TH2D("uplanededxresrg", "U Plane;Residual Range (cm); dE/dx", 300, 0, 300, 100, 0, 10);
  TH2D* vplanededxresrg = new TH2D("vplanededxresrg", "V Plane;Residual Range (cm); dE/dx", 300, 0, 300, 100, 0, 10);
  TH2D* xplanededxresrg = new TH2D("xplanededxresrg", "X Plane;Residual Range (cm); dE/dx", 300, 0, 300, 100, 0, 10);
  TH2D* uplanededxthetaxz = new TH2D("uplanededxthetaxz", "U Plane;#theta_{XZ}; dE/dx", 100, -3.15, 3.15, 100, 0, 10);
  TH2D* vplanededxthetaxz = new TH2D("vplanededxthetaxz", "V Plane;#theta_{XZ}; dE/dx", 100, -3.15 ,3.15, 100, 0, 10);
  TH2D* xplanededxthetaxz = new TH2D("xplanededxthetaxz", "X Plane;#theta_{XZ}; dE/dx", 100, -3.15, 3.15, 100, 0, 10);
  TH2D* uplanededxthetayz = new TH2D("uplanededxthetayz", "U Plane;#theta_{YZ}; dE/dx", 100, -3.15, 3.15, 100, 0, 10);
  TH2D* vplanededxthetayz = new TH2D("vplanededxthetayz", "V Plane;#theta_{YZ}; dE/dx", 100, -3.15 ,3.15, 100, 0, 10);
  TH2D* xplanededxthetayz = new TH2D("xplanededxthetayz", "X Plane;#theta_{YZ}; dE/dx", 100, -3.15, 3.15, 100, 0, 10);

  tree->Draw("trackdQdxByHitPlane0:trackResRangeByHitPlane0 >> uplanedqdxresrg");
  tree->Draw("trackdQdxByHitPlane1:trackResRangeByHitPlane1 >> vplanedqdxresrg");
  tree->Draw("trackdQdxByHitPlane2:trackResRangeByHitPlane2 >> xplanedqdxresrg");
  tree->Draw("trackdEdxByHitPlane0:trackResRangeByHitPlane0 >> uplanededxresrg");
  tree->Draw("trackdEdxByHitPlane1:trackResRangeByHitPlane1 >> vplanededxresrg");
  tree->Draw("trackdEdxByHitPlane2:trackResRangeByHitPlane2 >> xplanededxresrg");
  tree->Draw("trackdEdxByHitPlane0:trackThetaXZ >> uplanededxthetaxz");
  tree->Draw("trackdEdxByHitPlane1:trackThetaXZ >> vplanededxthetaxz");
  tree->Draw("trackdEdxByHitPlane2:trackThetaXZ >> xplanededxthetaxz");
  tree->Draw("trackdEdxByHitPlane0:trackThetaYZ >> uplanededxthetayz");
  tree->Draw("trackdEdxByHitPlane1:trackThetaYZ >> vplanededxthetayz");
  tree->Draw("trackdEdxByHitPlane2:trackThetaYZ >> xplanededxthetayz");

  plot(uplanedqdxresrg, "uplanedqdxresrg");
  plot(vplanedqdxresrg, "vplanedqdxresrg");
  plot(xplanedqdxresrg, "xplanedqdxresrg");
  plot(uplanededxresrg, "uplanededxresrg");
  plot(vplanededxresrg, "vplanededxresrg");
  plot(xplanededxresrg, "xplanededxresrg");
  plot(uplanededxthetaxz, "uplanededxthetaxz");
  plot(vplanededxthetaxz, "vplanededxthetaxz");
  plot(xplanededxthetaxz, "xplanededxthetaxz");
  plot(uplanededxthetayz, "uplanededxthetayz");
  plot(vplanededxthetayz, "vplanededxthetayz");
  plot(xplanededxthetayz, "xplanededxthetayz");

}
