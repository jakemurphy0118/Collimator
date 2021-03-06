{
  TFile f("1.root");
  gROOT->SetBatch(kTRUE);
  TCanvas *c1 = new TCanvas("c1","",800,600);

  // Vertex Histograms
  hk1Spectrum->Draw("hist");
  //  c1->SaveAs("hk1Spectrum.png");
  c1->SaveAs("1.png");
  hk2Spectrum->Draw("hist");
  c1->SaveAs("2.png");
  //  c1->SaveAs("hk2Spectrum.png","png");
  hk3Spectrum->Draw("hist");
  c1->SaveAs("3.png");
  //  c1->SaveAs("hk3Spectrum.png","png");
  hSumSpectrum->Draw("hist");
  c1->SaveAs("4.png");
  //  c1->SaveAs("hSumSpectrum.png","png");
  hCosTheta->Draw("hist");
  //  c1->SaveAs("hCosTheta.png","png");
  c1->SaveAs("5.png");
  hTheta->Draw("hist");
  c1->SaveAs("6.png");
  //  c1->SaveAs("hTheta.png","png");
  hPsi->Draw("hist");
  c1->SaveAs("7.png");
  //c1->SaveAs("hPsi.png","png");
  hPhi->Draw("hist");
  c1->SaveAs("8.png");
  //  c1->SaveAs("hPhi.png","png");
  hSin2Theta->Draw("hist");
  c1->SaveAs("9.png");
  //  c1->SaveAs("hSin2Theta.png","png");
  hCosPhi->Draw("hist");
  c1->SaveAs("10.png");
  //  c1->SaveAs("hCosPhi.png","png");
  hSinPsi->Draw("hist");
  c1->SaveAs("11.png");
  //  c1->SaveAs("hSinPsi.png","png");
  hAP1->Draw("hist");
  c1->SaveAs("12.png");
  
  // Outside of Source Histograms
  hk1Spectrum_wS->Draw("hist");
  c1->SaveAs("13.png");
  //  c1->SaveAs("hk1Spectrum_wS.png","png");  
  hk2Spectrum_wS->Draw("hist");
  c1->SaveAs("14.png");
  //  c1->SaveAs("hk2Spectrum_wS.png","png");
  hk3Spectrum_wS->Draw("hist");
  c1->SaveAs("15.png");
  //  c1->SaveAs("hk3Spectrum_wS.png","png");
  hSumSpectrum_wS->Draw("hist");
  c1->SaveAs("16.png");
  //  c1->SaveAs("hSumSpectrum_wS.png","png");
  hCosTheta_wS->Draw("hist");
  c1->SaveAs("17.png");
  //  c1->SaveAs("hCosTheta_wS.png","png");
  hTheta_wS->Draw("hist");
  c1->SaveAs("18.png");
  //  c1->SaveAs("hTheta_wS.png","png");
  hPsi_wS->Draw("hist");
  c1->SaveAs("19.png");
  //  c1->SaveAs("hPsi_wS.png","png");
  hPhi_wS->Draw("hist");
  c1->SaveAs("20.png");
  //  c1->SaveAs("hPhi_wS.png","png");
  hSin2Theta_wS->Draw("hist");
  c1->SaveAs("21.png");
  //  c1->SaveAs("hSin2Theta_wS.png","png");
  hCosPhi_wS->Draw("hist");
  c1->SaveAs("22.png");
  //  c1->SaveAs("hCosPhi_wS.png","png");
  hSinPsi_wS->Draw("hist");
  c1->SaveAs("23.png");
  //  c1->SaveAs("hSinPsi_wS.png","png");
  hAP1_wS->Draw("hist");
  c1->SaveAs("24.png");
  //  c1->SaveAs("hAP1_wS.png","png");

  // APEX Reconstruction Histograms

  hk1Spectrum_APEX->Draw("hist");
  c1->SaveAs("25.png");
  //  c1->SaveAs("hk1Spectrum_APEX.png","png");
  hk2Spectrum_APEX->Draw("hist");
  c1->SaveAs("26.png");
  //  c1->SaveAs("hk2Spectrum_APEX.png","png");
  hk3Spectrum_APEX->Draw("hist");
  c1->SaveAs("27.png");
  //  c1->SaveAs("hk3Spectrum_APEX.png","png");
  hSumSpectrum_APEX->Draw("hist");
  c1->SaveAs("28.png");
  //  c1->SaveAs("hSumSpectrum_APEX.png","png");
  hCosTheta_APEX->Draw("hist");
  c1->SaveAs("29.png");
  //c1->SaveAs("hCosTheta_APEX.png","png");
  hTheta_APEX->Draw("hist");
  c1->SaveAs("30.png");
  //  c1->SaveAs("hTheta_APEX.png","png");
  hPsi_APEX->Draw("hist");
  c1->SaveAs("31.png");
  //  c1->SaveAs("hPsi_APEX.png","png");
  hPhi_APEX->Draw("hist");
  c1->SaveAs("32.png");
  //  c1->SaveAs("hPhi_APEX.png","png");
  hSin2Theta_APEX->Draw("hist");
  c1->SaveAs("32.png");
  //  c1->SaveAs("hSin2Theta_APEX.png","png");
  hCosPhi_APEX->Draw("hist");
  c1->SaveAs("33.png");
  //  c1->SaveAs("hCosPhi_APEX.png","png");
  hSinPsi_APEX->Draw("hist");
  c1->SaveAs("34.png");
  //  c1->SaveAs("hSinPsi_APEX.png","png");
  hAP1_APEX->Draw("hist");
  c1->SaveAs("34.png");
  //  c1->SaveAs("hAP1_APEX.png","png");
  // Set to log scale
  totalEnergy_APEX->Draw("hist");
  c1->SetLogy();
  c1->SaveAs("35.png");
  //  c1->SaveAs("totalEnergy_APEX.png","png");
  numBarsHit->Draw("hist");
  c1->SetLogy(0);
  c1->SaveAs("36.png");
  //  c1->SaveAs("numBarsHit.png","png");
  hitbars->Draw("hist");
  c1->SaveAs("37.png");
  //  c1->SaveAs("hitbars.png","png");
  hitTime->Draw("hist");
  c1->SaveAs("38.png");
  //c1->SaveAs("hitTime.png","png");
  totalEnergy_1Bar->Draw("hist");
  c1->SaveAs("39.png");
  //c1->SaveAs("totalEnergy_1Bar.png","png");
  totalEnergy_2Bars->Draw("hist");
  c1->SaveAs("40.png");
  //  c1->SaveAs("totalEnergy_2Bars.png","png");
  totalEnergy_3Bars->Draw("hist");
  c1->SaveAs("41.png");
  //  c1->SaveAs("totalEnergy_3Bars.png","png");
  totalEnergy_4Bars->Draw("hist");
  c1->SaveAs("42.png");
  //  c1->SaveAs("totalEnergy_4Bars.png","png");
  totalEnergy_5Bars->Draw("hist");
  c1->SaveAs("43.png");
  //  c1->SaveAs("totalEnergy_5Bars.png","png");
  totalEnergy_6Bars->Draw("hist");
  c1->SaveAs("44.png");
  //  c1->SaveAs("totalEnergy_6Bars.png","png");
  totalEnergy_7Bars->Draw("hist");
  c1->SaveAs("45.png");
  //  c1->SaveAs("totalEnergy_7Bars.png","png");
  totalEnergy_8Bars->Draw("hist");
  c1->SaveAs("46.png");
  //  c1->SaveAs("totalEnergy_8Bars.png","png");

  hS2TCP->Draw("colz");
  c1->SaveAs("47.png");
  //  c1->SaveAs("hS2TCP.png","png");
  hS2TCP_wS->Draw("colz");
  c1->SaveAs("48.png");
  //c1->SaveAs("hS2TCP_wS.png","png");
  hS2TCP_APEX->Draw("colz");
  c1->SaveAs("49.png");
  //  c1->SaveAs("hS2TCP_APEX.png","png");
  k1Energy_vs_k2Energy->Draw("colz");
  c1->SaveAs("50.png");
  //  c1->SaveAs("k1Energy_vs_k2Energy.png","png");
  theta_v_phi->Draw("colz");
  c1->SaveAs("51.png");
  //  c1->SaveAs("theta_v_phi.png","png");
  theta_v_psi->Draw("colz");
  c1->SaveAs("52.png");
  //  c1->SaveAs("theta_v_psi.png","png");
  phi_v_psi->Draw("colz");
  c1->SaveAs("53.png");
  //  c1->SaveAs("phi_v_psi.png","png");

  asymmetryParameter = -(hAP1->Integral(1,100) - hAP1->Integral(101,200))/hAP1->Integral();
  asymmetryParameter_wS = -(hAP1_wS->Integral(1,100) - hAP1_wS->Integral(101,200))/hAP1_wS->Integral();
  asymmetryParameter_APEX = -(hAP1_APEX->Integral(1,100) - hAP1_APEX->Integral(101,200))/hAP1_APEX->Integral();
  std::cout << "Asymmetry Parameter at Vertex: " << asymmetryParameter << " +- " << 1./TMath::Sqrt(hAP1->Integral()) << std::endl;
  std::cout << "Asymmetry Parameter at Outside Source Holder: " << asymmetryParameter_wS << " +- " << 1./TMath::Sqrt(hAP1_wS->Integral()) << std::endl;
  std::cout << "Asymmetry Parameter with APEX: " << asymmetryParameter_APEX << " +- " << 1./TMath::Sqrt(hAP1_APEX->Integral()) << std::endl;
}
