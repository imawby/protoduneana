void compare_v_plane(){

  TFile* _file_min2V = new TFile("2dplotsandgraphs_-3U_-2V.root", "read");
  TFile* _file_min3V = new TFile("2dplotsandgraphs_-3U_-3V.root", "read");
  TFile* _file_min4V = new TFile("2dplotsandgraphs_-3U_-4V.root", "read");

  TGraphErrors* gr_min2V = (TGraphErrors*)_file_min2V->Get("trackHitCaloChargeDep_trackHitView==1_trackHitAPABuildNumber>-1_mean"); 
  TGraphErrors* gr_min3V = (TGraphErrors*)_file_min3V->Get("trackHitCaloChargeDep_trackHitView==1_trackHitAPABuildNumber>-1_mean");
  TGraphErrors* gr_min4V = (TGraphErrors*)_file_min4V->Get("trackHitCaloChargeDep_trackHitView==1_trackHitAPABuildNumber>-1_mean");

  gr_min3V->SetMarkerStyle(20);
  gr_min2V->SetMarkerColor(kAzure+1);
  gr_min2V->SetLineColor(kAzure+1);
  gr_min2V->SetMarkerStyle(23);
  gr_min4V->SetMarkerColor(kGreen+1);
  gr_min4V->SetLineColor(kGreen+1);
  gr_min4V->SetMarkerStyle(22);

  gr_min3V->SetTitle(";Wire Tension (N);dQ/dx (ADC*/cm)");

  TCanvas *c1 = new TCanvas("c1", "c1");

  gr_min3V->Draw("");
  gr_min2V->Draw("lp same");
  gr_min4V->Draw("lp same");

  TLegend *leg = new TLegend(0.18, 0.70, 0.35, 0.85);
  leg->AddEntry(gr_min3V, "Nominal Map");
  leg->AddEntry(gr_min2V, "-1 Wire");
  leg->AddEntry(gr_min4V, "+1 Wire");
  leg->Draw("same");

  c1->SaveAs("vPlaneComparison.png");
  
}
