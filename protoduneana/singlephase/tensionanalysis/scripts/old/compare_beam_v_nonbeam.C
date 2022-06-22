void compare_beam_v_nonbeam(){

  TFile* _file_min2V = new TFile("2dplotsandgraphs.root", "read");

  TGraphErrors* gr_min2V_1 = (TGraphErrors*)_file_min2V->Get("trackHitCaloChargeDep_trackHitView>-1_trackHitAPABuildNumber==1_mean"); 
  TGraphErrors* gr_min2V_2 = (TGraphErrors*)_file_min2V->Get("trackHitCaloChargeDep_trackHitView>-1_trackHitAPABuildNumber==2_mean"); 
  TGraphErrors* gr_min2V_3 = (TGraphErrors*)_file_min2V->Get("trackHitCaloChargeDep_trackHitView>-1_trackHitAPABuildNumber==3_mean"); 
  TGraphErrors* gr_min2V_4 = (TGraphErrors*)_file_min2V->Get("trackHitCaloChargeDep_trackHitView>-1_trackHitAPABuildNumber==4_mean"); 
  TGraphErrors* gr_min2V_5 = (TGraphErrors*)_file_min2V->Get("trackHitCaloChargeDep_trackHitView>-1_trackHitAPABuildNumber==5_mean"); 
  TGraphErrors* gr_min2V_6 = (TGraphErrors*)_file_min2V->Get("trackHitCaloChargeDep_trackHitView>-1_trackHitAPABuildNumber==6_mean"); 

  gr_min2V_1->SetMarkerColor(kAzure+1);
  gr_min2V_1->SetLineColor(kAzure+1);
  gr_min2V_1->SetMarkerStyle(23);
  gr_min2V_2->SetMarkerColor(kAzure+1);
  gr_min2V_2->SetLineColor(kAzure+1);
  gr_min2V_2->SetMarkerStyle(23);
  gr_min2V_3->SetMarkerColor(kAzure+1);
  gr_min2V_3->SetLineColor(kAzure+1);
  gr_min2V_3->SetMarkerStyle(23);
  gr_min2V_4->SetMarkerColor(kBlack);
  gr_min2V_4->SetLineColor(kBlack);
  gr_min2V_4->SetMarkerStyle(23);
  gr_min2V_5->SetMarkerColor(kBlack);
  gr_min2V_5->SetLineColor(kBlack);
  gr_min2V_5->SetMarkerStyle(23);
  gr_min2V_6->SetMarkerColor(kBlack);
  gr_min2V_6->SetLineColor(kBlack);
  gr_min2V_6->SetMarkerStyle(23);

  gr_min2V_1->SetTitle(";Wire Tension (N);dQ/dx (ADC*/cm)");
  gr_min2V_1->GetYaxis()->SetRangeUser(0,500);

  TCanvas *c1 = new TCanvas("c1", "c1");

  gr_min2V_1->Draw("");
  gr_min2V_2->Draw("lp same");
  gr_min2V_3->Draw("lp same");
  gr_min2V_4->Draw("lp same");
  gr_min2V_5->Draw("lp same");
  gr_min2V_6->Draw("lp same");

  TLegend *leg = new TLegend(0.18, 0.70, 0.35, 0.85);
  leg->AddEntry(gr_min2V_1, "Beam Side APAs");
  leg->AddEntry(gr_min2V_4, "Non-Beam Side APAs");
  leg->Draw("same");

  c1->SaveAs("BeamNonBeamSideComparison.png");
  
}
