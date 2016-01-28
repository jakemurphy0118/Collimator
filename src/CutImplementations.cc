//
// CutImplementations
//
#include "CutImplementations.hh"
#include <TMath.h>
#include "AnalysisManager.hh"
#include <TVector3.h>

CutImplementations::CutImplementations()
{
}

CutImplementations::~CutImplementations()
{
}

void CutImplementations::DoAnalysisCuts(G4double totalEnergy, G4double time, G4double magVectorSum, G4double k1Energy, G4double k2Energy, G4double k3Energy, std::vector<G4double> k1GammaRayData, std::vector<G4double> k2GammaRayData, std::vector<G4double> k3GammaRayData, std::vector<G4int> threeBars, G4int hitBarsSetSize)
{
  G4String identifier;

  // No cuts
  identifier="nocuts";
  RyanReconstruction(k1GammaRayData,k2GammaRayData,k3GammaRayData,threeBars,hitBarsSetSize,identifier);
  YamazakiReconstruction(k1GammaRayData,k2GammaRayData,k3GammaRayData,threeBars,hitBarsSetSize,identifier);

  // Energy Cut
  G4double totalEnergyCutMax=2500;
  G4double totalEnergyCutMin=2100;
  if ( (totalEnergy<totalEnergyCutMax)&&(totalEnergy>totalEnergyCutMin) ) {
    identifier="cut1";
    RyanReconstruction(k1GammaRayData,k2GammaRayData,k3GammaRayData,threeBars,hitBarsSetSize,identifier);
    YamazakiReconstruction(k1GammaRayData,k2GammaRayData,k3GammaRayData,threeBars,hitBarsSetSize,identifier);
  }

  G4double numBarCut=0;
  // Timing Cut
  if ( time==numBarCut ) {
    identifier="cut2";
    RyanReconstruction(k1GammaRayData,k2GammaRayData,k3GammaRayData,threeBars,hitBarsSetSize,identifier);
    YamazakiReconstruction(k1GammaRayData,k2GammaRayData,k3GammaRayData,threeBars,hitBarsSetSize,identifier);    
  }
  
}

void CutImplementations::RyanReconstruction(std::vector<G4double> k1Data, std::vector<G4double> k2Data, std::vector<G4double> k3Data, std::vector<G4int> Bars, G4double numberOfBars, G4String identifier)
{
  // Get Analysis Manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  // Initialize variables
  double xi1, xi2, xi3;
  double delta;
  double k1_mag, k2_mag, k3_mag;
  double sinxi1, cosxi1, sinxi2, cosxi2, sinxi3, cosxi3, sindelta;
  
  G4double k1_energy=k1Data[0];
  G4double k2_energy=k2Data[0];
  G4double k3_energy=k3Data[0];

  TVector3 k1(k1Data[1],k1Data[2],k1Data[3]);
  TVector3 k2(k2Data[1],k2Data[2],k2Data[3]);
  TVector3 k3(k3Data[1],k3Data[2],k3Data[3]);

  IDNum = analysisManager->GetH2Id("k1zPosition_v_k2zPosition",1);
  analysisManager->FillH2(IDNum,k1Data[3]*275,k2Data[3]*275);
  
  //===========================
  // Ryan Reconstruction
  //===========================

  k1_mag = k1.Mag(); k2_mag = k2.Mag(); k3_mag = k3.Mag();
  sinxi1 = k1.Z()/k1_mag;
  xi1 = TMath::ASin(sinxi1);
  cosxi1 = TMath::Sqrt(1-k1.Z()*k1.Z()/k1_mag/k1_mag);
  sinxi2 = k2.Z()/k2_mag;
  cosxi2 = TMath::Sqrt(1-k2.Z()*k2.Z()/k2_mag/k2_mag);
  xi2 = TMath::ASin(sinxi2);
  sinxi3 = k3.Z()/k3_mag;
  cosxi3 = TMath::Sqrt(1-k3.Z()*k3.Z()/k3_mag/k3_mag);
  xi3 = TMath::ASin(sinxi3);
  sindelta = (k1.X()*k2.Y()-k1.Y()*k2.X())/cosxi1/cosxi2/k1_mag/k2_mag;
  delta = TMath::ASin(sindelta);
    
  //Copy Numbers run from 0 to 23
  G4int k1Bar=Bars[0];
  G4int k2Bar=Bars[1];
  G4int oppositeBar=(k1Bar+12)%24;

  //Check if k2Bar is
  if (oppositeBar==(k2Bar+1)%24) {
    IDNum = analysisManager->GetH1Id("cosXi_rBP1",1);   
    analysisManager->FillH1(IDNum,TMath::Cos(xi1));
  }
  if (oppositeBar==(k2Bar+2)%24) {
    IDNum = analysisManager->GetH1Id("cosXi_rBP2",1);   
    analysisManager->FillH1(IDNum,TMath::Cos(xi1));
  }
  if (oppositeBar==(k2Bar+3)%24) {
    IDNum = analysisManager->GetH1Id("cosXi_rBP3",1);   
    analysisManager->FillH1(IDNum,TMath::Cos(xi1));
  }
  if (oppositeBar==(k2Bar+4)%24) {
    IDNum = analysisManager->GetH1Id("cosXi_rBP4",1);   
    analysisManager->FillH1(IDNum,TMath::Cos(xi1));
  }
  if (oppositeBar==(k2Bar-1)%24) {
    IDNum = analysisManager->GetH1Id("cosXi_rBM1",1);   
    analysisManager->FillH1(IDNum,TMath::Cos(xi1));
  }
  if (oppositeBar==(k2Bar-2)%24) {
    IDNum = analysisManager->GetH1Id("cosXi_rBM2",1);   
    analysisManager->FillH1(IDNum,TMath::Cos(xi1));
  }
  if (oppositeBar==(k2Bar-3)%24) {
    IDNum = analysisManager->GetH1Id("cosXi_rBM3",1);   
    analysisManager->FillH1(IDNum,TMath::Cos(xi1));
  }
  if (oppositeBar==(k2Bar-4)%24) {
    IDNum = analysisManager->GetH1Id("cosXi_rBM4",1);   
    analysisManager->FillH1(IDNum,TMath::Cos(xi1));
  }
  
  G4double Xmomentum=std::abs(k1.X()+k2.X()+k3.X());
  G4double Ymomentum=std::abs(k1.Y()+k2.Y()+k3.Y());
  G4double Zmomentum=std::abs(k1.Z()+k2.Z()+k3.Z());

  if (numberOfBars==4) {

    if ( identifier=="nocuts" ) {
      
      IDNum = analysisManager->GetH1Id("hk1Spectrum_APEX",1);
      analysisManager->FillH1(IDNum,k1_energy);
      IDNum = analysisManager->GetH1Id("hk2Spectrum_APEX",1);
      analysisManager->FillH1(IDNum,k2_energy);
      IDNum = analysisManager->GetH1Id("hk3Spectrum_APEX",1);
      analysisManager->FillH1(IDNum,k3_energy);
      IDNum = analysisManager->GetH1Id("hSumSpectrum_APEX",1);
      analysisManager->FillH1(IDNum,k1_energy+k2_energy+k3_energy);
      IDNum = analysisManager->GetH1Id("hXi1_APEX",1);
      analysisManager->FillH1(IDNum,xi1);
      IDNum = analysisManager->GetH1Id("hXi2_APEX",1);
      analysisManager->FillH1(IDNum,xi2);
      IDNum = analysisManager->GetH1Id("hXi3_APEX",1);
      analysisManager->FillH1(IDNum,xi3);
      //      IDNum = analysisManager->GetH1Id("sinXi_1_APEX",1);
      //      analysisManager->FillH1(IDNum,TMath::Sin(xi1));
      //      IDNum = analysisManager->GetH1Id("cosXi_2_APEX",1);
      //      analysisManager->FillH1(IDNum,TMath::Cos(xi2));
      //      IDNum = analysisManager->GetH1Id("sin_delta_APEX",1);
      //      analysisManager->FillH1(IDNum,sindelta);
      //      IDNum = analysisManager->GetH1Id("sinxi_1xsindelta_APEX",1);
      //      analysisManager->FillH1(IDNum,TMath::Sin(xi1)*TMath::Sin(delta));    
      IDNum = analysisManager->GetH1Id("rh_asymmetry_APEX",1);
      analysisManager->FillH1(IDNum,TMath::Sin(xi1)*TMath::Cos(xi1)*TMath::Cos(xi1)*TMath::Sin(delta));
      IDNum = analysisManager->GetH1Id("delta_APEX",1);
      analysisManager->FillH1(IDNum,delta);
      
    }

    if ( identifier=="cut1" ) {
      IDNum = analysisManager->GetH1Id("hk1Spectrum_APEX_Cut1",1);
      analysisManager->FillH1(IDNum,k1_energy);
      IDNum = analysisManager->GetH1Id("hk2Spectrum_APEX_Cut1",1);
      analysisManager->FillH1(IDNum,k2_energy);
      IDNum = analysisManager->GetH1Id("hk3Spectrum_APEX_Cut1",1);
      analysisManager->FillH1(IDNum,k3_energy);
      IDNum = analysisManager->GetH1Id("hSumSpectrum_APEX_Cut1",1);
      analysisManager->FillH1(IDNum,k1_energy+k2_energy+k3_energy);
      IDNum = analysisManager->GetH1Id("hXi1_APEX_Cut1",1);
      analysisManager->FillH1(IDNum,xi1);
      IDNum = analysisManager->GetH1Id("hXi2_APEX_Cut1",1);
      analysisManager->FillH1(IDNum,xi2);
      IDNum = analysisManager->GetH1Id("hXi3_APEX_Cut1",1);
      analysisManager->FillH1(IDNum,xi3);
      //      IDNum = analysisManager->GetH1Id("sinXi_1_APEX",1);
      //      analysisManager->FillH1(IDNum,TMath::Sin(xi1));
      //      IDNum = analysisManager->GetH1Id("cosXi_2_APEX",1);
      //      analysisManager->FillH1(IDNum,TMath::Cos(xi2));
      //      IDNum = analysisManager->GetH1Id("sin_delta_APEX",1);
      //      analysisManager->FillH1(IDNum,sindelta);
      //      IDNum = analysisManager->GetH1Id("sinxi_1xsindelta_APEX",1);
      //      analysisManager->FillH1(IDNum,TMath::Sin(xi1)*TMath::Sin(delta));    
      IDNum = analysisManager->GetH1Id("rh_asymmetry_APEX",1);
      analysisManager->FillH1(IDNum,TMath::Sin(xi1)*TMath::Cos(xi1)*TMath::Cos(xi1)*TMath::Sin(delta));
      IDNum = analysisManager->GetH1Id("delta_APEX",1);
      analysisManager->FillH1(IDNum,delta);     
    }
    
  }

}

void CutImplementations::YamazakiReconstruction(std::vector<G4double> k1Data, std::vector<G4double> k2Data, std::vector<G4double> k3Data, std::vector<G4int> Bars, G4double numberOfBars, G4String identifier)
{
  // Get Analysis Manager 
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  // Initialize variables
  TVector3 normal;
  double theta;
  double psi;
  double phi;
  TVector3 sprojection;
  double x1;
  double y1;
  double z1;
  double x2;
  double y2;
  double z2;
  double pos_z;
  double pos_y;
  double pos_x;

  // Get energy and momentum of the gamma rays
  G4double k1_energy=k1Data[0];
  G4double k2_energy=k2Data[0];
  G4double k3_energy=k3Data[0];

  TVector3 k1(k1Data[1],k1Data[2],k1Data[3]);
  TVector3 k2(k2Data[1],k2Data[2],k2Data[3]);
  TVector3 k3(k3Data[1],k3Data[2],k3Data[3]);

  IDNum = analysisManager->GetH2Id("k1zPosition_v_k2zPosition",1);
  analysisManager->FillH2(IDNum,k1Data[3]*275,k2Data[3]*275);
  
  //===========================
  // Yamazaki Reconstruction
  //===========================
  normal = k1.Cross(k2);
  theta = TMath::ACos(normal.CosTheta());
  psi = k1.Angle(k2);
  sprojection = normal;
  sprojection.Rotate(TMath::Pi()*0.5, normal.Cross(TVector3(0,0,1))); // projection of z-axis (S) onto decay plane.
  phi = k1.Angle(sprojection);

  //====================================
  // Fill Histograms for 4 bars hit only
  //====================================
  
  if (numberOfBars==4) {

    if (identifier=="nocuts") {
    
      IDNum = analysisManager->GetH1Id("hk1Spectrum_APEX",1);
      analysisManager->FillH1(IDNum,k1_energy);
      IDNum = analysisManager->GetH1Id("hk2Spectrum_APEX",1);
      analysisManager->FillH1(IDNum,k2_energy);
      IDNum = analysisManager->GetH1Id("hk3Spectrum_APEX",1);
      analysisManager->FillH1(IDNum,k3_energy);
      IDNum = analysisManager->GetH1Id("hSumSpectrum_APEX",1);
      analysisManager->FillH1(IDNum,k1_energy+k2_energy+k3_energy);
      //      IDNum = analysisManager->GetH1Id("hCosTheta_APEX",1);
      //      analysisManager->FillH1(IDNum,k1.CosTheta());
      IDNum = analysisManager->GetH1Id("hTheta_APEX",1);
      analysisManager->FillH1(IDNum,theta);
      IDNum = analysisManager->GetH1Id("hPsi_APEX",1);
      analysisManager->FillH1(IDNum,psi);    
      IDNum = analysisManager->GetH1Id("hPhi_APEX",1);
      analysisManager->FillH1(IDNum,phi);
      //      IDNum = analysisManager->GetH1Id("hSin2Theta_APEX",1);
      //      analysisManager->FillH1(IDNum,TMath::Sin(2.*theta));
      //      IDNum = analysisManager->GetH1Id("hCosPhi_APEX",1);
      //      analysisManager->FillH1(IDNum,TMath::Cos(phi));
      //      IDNum = analysisManager->GetH1Id("hSinPsi_APEX",1);
      //      analysisManager->FillH1(IDNum,TMath::Sin(psi));
      IDNum = analysisManager->GetH1Id("hAP1_APEX",1);
      analysisManager->FillH1(IDNum,TMath::Cos(phi) * TMath::Sin(2.*theta));
      IDNum = analysisManager->GetH2Id("theta_v_phi_APEX",1);
      analysisManager->FillH2(IDNum,theta,phi);
      IDNum = analysisManager->GetH2Id("theta_v_psi_APEX",1);
      analysisManager->FillH2(IDNum,theta,psi);
      IDNum = analysisManager->GetH2Id("phi_v_psi_APEX",1);
      analysisManager->FillH2(IDNum,phi,psi);
    }

    if ( identifier == "cut1" ) {
      // Fill other histograms
    }
  }
}


