struct GraphInfo{

  int View = -1;
  int APA  = -1;

};

GraphInfo parse_name(std::string thisName){

  GraphInfo thisGraphInfo;

  std::string hv = "trackHitView";
  size_t viewPos = thisName.find(hv);
  std::string viewStr;
  if (thisName.substr(viewPos+hv.length()+1,1) == "-")
    viewStr = thisName.substr(viewPos+hv.length()+1, 2);
  else 
    viewStr = thisName.substr(viewPos+hv.length()+2, 1);
  thisGraphInfo.View = std::stoi(viewStr, nullptr);

  std::string hapa = "trackHitAPABuildNumber";
  size_t apaPos = thisName.find(hapa);
  std::string apaStr;
  if (thisName.substr(apaPos+hapa.length()+1,1) == "-")
    apaStr = thisName.substr(apaPos+hapa.length()+1, 2);
  else 
    apaStr = thisName.substr(apaPos+hapa.length()+2, 1);

  thisGraphInfo.APA = std::stoi(apaStr, nullptr);

  return thisGraphInfo;

}

void plot_from_file(){
  std::cout << __LINE__ << std::endl;

  std::vector< std::string > type = {
    "median",
    "fitmpv",
    "mean",
    "rms"
  };
  std::cout << __LINE__ << std::endl;

  std::vector< std::string >  vars = {
    //"trackHitCaloEnergyDep",
    //"trackHitCaloChargeDep",
    //"trackThetaXZ",
    //"trackThetaYZ"
    //"trackHitRMS",
    //"trackHitPeakTime",
    //"trackHitPeakAmplitude",
    "trackHitIntegral"
  };
  std::cout << __LINE__ << std::endl;

  std::vector<std::vector<int>> yrange ={
    //{1, 3},
    //{0, 150},
    //{1,3},
    //{1,3}
    //{0, 20},
    //{0, 10000},
    //{0, 100},
    {0, 100}
  };
  std::cout << __LINE__ << std::endl;

  std::vector< std::string > labels = {
    //"Hit dE/dx (MeV/cm)",
    //"Hit dQ/dx (ADC*/cm)",
    //"Track #theta_{XZ}",
    //"Track #theta_{YZ}"
    //"Hit RMS",
    //"Hit Peak Time (#mu s)",
    //"Hit Peak Amplitude (ADC*)",
    "Hit Integral"
  };
  std::cout << __LINE__ << std::endl;

  std::vector< int > apas = {
    -1
    //1,
    //2,
    //3,
    //4,
    //5,
    //6
  };


  std::cout << __LINE__ << std::endl;

  TFile* thisFile = (TFile*)_file0;
  std::cout << __LINE__ << std::endl;

  // get list of things in the file
  TList* thisList = thisFile->GetListOfKeys();
  std::cout << __LINE__ << std::endl;

  std::vector<TGraphErrors*> thisVec(4);
  TObject* thisObj;
  std::cout << __LINE__ << std::endl;

  //// loop vars
  for (int iv = 0; iv < vars.size(); iv++){

    std::cout << vars.at(iv) << std::endl;

    // loop types
    for (int it = 0; it < type.size(); it++){

      std::cout << vars.at(iv) << " " << type.at(it) << std::endl;

      // loop APAs
      for (int ia = 0; ia < apas.size(); ia++){

        std::cout << ia << std::endl;

        // loop over the stuff in the file and if the object doesn't
        // inherit from a TGraphErrors, ignore it
        for (int i = 0; i < thisList->GetSize(); i++){

          TKey* thisKey = (TKey*)thisList->At(i);

          thisObj = thisKey->ReadObj();


          if (!thisObj->InheritsFrom("TGraphErrors")) continue;


          TGraphErrors* thisTGE = (TGraphErrors*)thisFile->Get(thisObj->GetName());

          std::string thisName = std::string(thisTGE->GetName());

          // pick out specific type of plot
          if (thisName.find(type.at(it)) != std::string::npos) {

            if (thisName.find(vars.at(iv)) != std::string::npos){


              GraphInfo thisGraphInfo = parse_name(thisName);


              if (thisGraphInfo.APA == apas.at(ia) && thisGraphInfo.View == -1)
                thisVec.at(3) = thisTGE;
              if (thisGraphInfo.APA == apas.at(ia) && thisGraphInfo.View == 0)
                thisVec.at(0) = thisTGE;
              if (thisGraphInfo.APA == apas.at(ia) && thisGraphInfo.View == 1)
                thisVec.at(1) = thisTGE;
              if (thisGraphInfo.APA == apas.at(ia) && thisGraphInfo.View == 2)
                thisVec.at(2) = thisTGE;

            }
          }
        }

        thisVec.at(0)->SetLineColor(kGreen+1);
        thisVec.at(0)->SetMarkerColor(kGreen+1);
        thisVec.at(0)->SetMarkerStyle(22);
        thisVec.at(1)->SetLineColor(kAzure+1);
        thisVec.at(1)->SetMarkerColor(kAzure+1);
        thisVec.at(1)->SetMarkerStyle(23);
        thisVec.at(2)->SetLineColor(kOrange+1);
        thisVec.at(2)->SetMarkerColor(kOrange+1);
        thisVec.at(2)->SetMarkerStyle(21);
        thisVec.at(3)->SetLineWidth(2);
        thisVec.at(3)->SetMarkerStyle(20);

        // plot!
        TCanvas* c1 = new TCanvas("c1", "c1");
        c1->SetLeftMargin(0.15);
        c1->SetBottomMargin(0.15);

        double minimum = std::min({
            *std::min_element(thisVec.at(0)->GetY(), thisVec.at(0)->GetY()+9), 
            *std::min_element(thisVec.at(1)->GetY(), thisVec.at(1)->GetY()+9), 
            *std::min_element(thisVec.at(2)->GetY(), thisVec.at(2)->GetY()+9), 
            *std::min_element(thisVec.at(3)->GetY(), thisVec.at(3)->GetY()+9)});

        double maximum = std::max({
            *std::max_element(thisVec.at(0)->GetY(), thisVec.at(0)->GetY()+9), 
            *std::max_element(thisVec.at(1)->GetY(), thisVec.at(1)->GetY()+9), 
            *std::max_element(thisVec.at(2)->GetY(), thisVec.at(2)->GetY()+9), 
            *std::max_element(thisVec.at(3)->GetY(), thisVec.at(3)->GetY()+9)});

        if (minimum < 10e-10) minimum = 0;

        thisVec.at(0)->GetYaxis()->SetRangeUser(yrange.at(iv).at(0),yrange.at(iv).at(1));
        //thisVec.at(0)->GetYaxis()->SetRangeUser(200,400);
        thisVec.at(0)->GetYaxis()->SetNdivisions(505);

        std::string axis;
        if (apas.at(ia) == -1){
          axis = 
            //std::string("All APAs") +
            std::string(";Wire Tension (N);")+
            type.at(it) +
            std::string(" ") +
            labels.at(iv);
        }
        else{
          axis = 
            std::string("APA ") +
            std::to_string(apas.at(ia)) +
            std::string(";Wire Tension (N);")+
            type.at(it) +
            std::string(" ") +
            labels.at(iv);
        }

        thisVec.at(0)->SetTitle(axis.c_str());
        thisVec.at(0)->GetYaxis()->SetTitleFont(43); 
        thisVec.at(0)->GetYaxis()->SetLabelFont(43); 
        thisVec.at(0)->GetXaxis()->SetTitleFont(43); 
        thisVec.at(0)->GetXaxis()->SetLabelFont(43); 
        thisVec.at(0)->GetYaxis()->SetTitleSize(28); 
        thisVec.at(0)->GetYaxis()->SetLabelSize(28); 
        thisVec.at(0)->GetXaxis()->SetTitleSize(28); 
        thisVec.at(0)->GetXaxis()->SetLabelSize(28); 
        thisVec.at(0)->GetXaxis()->SetTitleOffset(0.9); 
        thisVec.at(0)->GetYaxis()->SetTitleOffset(0.95);
        thisVec.at(0)->GetXaxis()->CenterTitle(); 
        thisVec.at(0)->GetYaxis()->CenterTitle(); 

        thisVec.at(0)->Draw("ap");
        for ( int i = 1; i < thisVec.size(); i++){
          thisVec.at(i)->Draw("same p");
        }

        TLegend* leg = new TLegend(0.18, 0.68, 0.5, 0.88);
        leg->AddEntry(thisVec.at(0), "U Plane");
        leg->AddEntry(thisVec.at(1), "V Plane");
        leg->AddEntry(thisVec.at(2), "Y Plane");
        leg->AddEntry(thisVec.at(3), "All Planes");
        leg->Draw("same");

        std::string plotName = 
          vars.at(iv) + 
          std::string("_") + 
          type.at(it) + 
          std::string("_apa_") + 
          std::to_string(apas.at(ia));

        c1->SaveAs((plotName+".png").c_str());
        c1->SaveAs((plotName+".pdf").c_str());

      }
    }
  }

  // now loop over TH2s in the file and change the range on the Y axis so
  // that it makes sense
  for (int i = 0; i < thisList->GetSize(); i++){

    TKey* thisKey = (TKey*)thisList->At(i);
    thisObj = thisKey->ReadObj();

    if (!thisObj->InheritsFrom("TH2D")) continue;

    TH2D* thisTH2 = (TH2D*)thisFile->Get(thisObj->GetName());
    std::string thisName = std::string(thisTH2->GetName());

    for (int j = 0; j < vars.size(); j++){
      if (thisName.find(vars.at(j)) != std::string::npos) {
        thisTH2->GetYaxis()->SetRangeUser(yrange.at(j).at(0), yrange.at(j).at(1));
      }
    }

    TCanvas* c2 = new TCanvas("c2", "c2");
    c2->SetLeftMargin(0.15);
    c2->SetBottomMargin(0.15);
    c2->SetRightMargin(0.12);
    c2->cd();

    thisTH2->Draw("colz");
    //TPaletteAxis *palette = (TPaletteAxis*)thisTH2->GetListOfFunctions()->FindObject("palette");
    //palette->SetX1NDC(0.89);
    //palette->SetX2NDC(0.93);
    //palette->SetY1NDC(0.12);
    //palette->SetY2NDC(0.94);
    gPad->Modified();
    gPad->Update();

    c2->SaveAs((thisName+std::string(".png")).c_str());

  }
}
