//--------------------------------------------------------------------
// RunAction.cc
//
// Description: Defines the run itself.
// Changes: 7/18/05 None yet
//--------------------------------------------------------------------
#include <iostream>
#include <TF2.h>
#include <TF1.h>
#include <TH2.h>
#include <TH1D.h>
#include <TFile.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TVector3.h>
#include <TROOT.h>
#include "G4UImanager.hh"
#include "G4Polyline.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4VisExecutive.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "Randomize.hh"
#include <time.h>
#include "RunActionMessenger.hh"
#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "EventAction.hh"
#include "G4GeneralParticleSource.hh"
#include "G4UIcommand.hh"
#include "TimeDistribution.hh"
#include "AnalysisManager.hh"
#include "PhysicsList.hh"

RunAction::RunAction(PrimaryGeneratorAction* gen,PhysicsList* physList)
: genAction(gen)
{
  runActMessenger = new RunActionMessenger(this);
  
  BField=physList->GetBFieldValue();

  // Set up timing histograms
  timedistr = new TimeDistribution();
  //Set B Field
  timingHistoPS = timedistr->buildHist("oPSDecay",5000);
  timingHistpPS = timedistr->buildHist("pPSDecay",5000);
  timingHistoPS_2Gamma = timedistr->buildHist("oPSDecay_2Gamma",5000);
  timingHistpPS_2Gamma = timedistr->buildHist("pPSDecay_2Gamma",5000);
  timingHistoPS_3Gamma = timedistr->buildHist("oPSDecay_3Gamma",5000);
  timingHistpPS_3Gamma = timedistr->buildHist("pPSDecay_3Gamma",5000);

  BRatio=timedistr->GetBranchingRatio("oPSDecay",5000);
  BRatio=timedistr->GetBranchingRatio("pPSDecay",5000);

  // Create analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetFirstHistoId(1);

  // Creating histograms
  analysisManager->CreateH1("hk1Spectrum","k1 Spectrum (Vertex)", 220, 0., 1.1); 
  analysisManager->CreateH1("hk2Spectrum","k2 Spectrum (Vertex)", 220, 0., 1.1); 
  analysisManager->CreateH1("hk3Spectrum","k3 Spectrum (Vertex)", 220, 0., 1.1); 
  analysisManager->CreateH1("hSumSpectrum","Sum Spectrum (Vertex)", 200, 0., 2.5); 
  //  analysisManager->CreateH1("hCosTheta","cos(theta) (Vertex)", 200, -1, 1); 
  analysisManager->CreateH1("hTheta","theta (Vertex)", 200, 0., TMath::Pi()); 
  analysisManager->CreateH1("hPsi","psi (Vertex)", 100, 0., 2.*TMath::Pi()); 
  analysisManager->CreateH1("hPhi","phi (Vertex)", 100, 0., TMath::Pi()); 
  //  analysisManager->CreateH1("hSin2Theta","sin(2theta) (Vertex)", 200, -1,1); 
  //  analysisManager->CreateH1("hCosPhi","cos(phi) (Vertex)", 200, -1,1); 
  //  analysisManager->CreateH1("hSinPsi","sin(psi) (Vertex)", 200, -1,1); 
  analysisManager->CreateH1("hXi1","xi 1 (Vertex)",200,-TMath::Pi()/2,TMath::Pi()/2); 
  analysisManager->CreateH1("hXi2","xi 2 (Vertex)",200,-TMath::Pi()/2,TMath::Pi()/2); 
  analysisManager->CreateH1("hXi3","xi 3 (Vertex)",200,-TMath::Pi()/2,TMath::Pi()/2); 
  //  analysisManager->CreateH1("sinXi_1","sin(xi1) (Vertex)",200,-1,1); 
  //  analysisManager->CreateH1("cosXi_2","cos(xi2) (Vertex)",200,-1,1);
  //  analysisManager->CreateH1("sin_delta","sin(delta) (Vertex)",200,-1,1); 
  //  analysisManager->CreateH1("sinxi_1xsindelta","sin(xi1)*sin(delta) (Vertex)",200,-1,1); 
  analysisManager->CreateH1("delta","delta (Vertex)",200,-TMath::Pi()/2,TMath::Pi()/2); 
  analysisManager->CreateH1("hAP1","sin(2 theta) * cos(phi) (Vertex)", 200, -1,1);
  analysisManager->CreateH1("PsSubType","PsSubType", 6, 0,6);
  
  //Create histograms (with source)
  analysisManager->CreateH1("hk1Spectrum_wS","k1 Spectrum (with source holder)", 220, 0., 1.1); 
  analysisManager->CreateH1("hk2Spectrum_wS","k2 Spectrum (with source holder)", 220, 0., 1.1); 
  analysisManager->CreateH1("hk3Spectrum_wS","k3 Spectrum (with source holder)", 220, 0., 1.1); 
  analysisManager->CreateH1("hSumSpectrum_wS","Sum Spectrum (with source holder)", 200, 0., 2.5); 
  //  analysisManager->CreateH1("hCosTheta_wS","cos(theta) (with source holder)", 200, -1, 1); 
  analysisManager->CreateH1("hTheta_wS","theta (with source holder)", 200, 0., TMath::Pi()); 
  analysisManager->CreateH1("hPsi_wS","psi (with source holder)", 100, 0., 2.*TMath::Pi()); 
  analysisManager->CreateH1("hPhi_wS","phi (with source holder)", 100, 0., TMath::Pi()); 
  //  analysisManager->CreateH1("hSin2Theta_wS","sin(2theta) (with source holder)", 200, -1,1); 
  //  analysisManager->CreateH1("hCosPhi_wS","cos(phi) (with source holder)", 200, -1,1); 
  //  analysisManager->CreateH1("hSinPsi_wS","sin(psi) (with source holder)", 200, -1,1); 
  analysisManager->CreateH1("hXi1_wS","xi 1 (with source holder)",200,-TMath::Pi()/2,TMath::Pi()/2); 
  analysisManager->CreateH1("hXi2_wS","xi 2 (with source holder)",200,-TMath::Pi()/2,TMath::Pi()/2); 
  analysisManager->CreateH1("hXi3_wS","xi 3 (with source holder)",200,-TMath::Pi()/2,TMath::Pi()/2); 
  //  analysisManager->CreateH1("sinXi_1_wS","sin(xi1) (with source holder)",200,-1,1); 
  //  analysisManager->CreateH1("cosXi_2_wS","cos(xi2) (with source holder)",200,-1,1); 
  //  analysisManager->CreateH1("sin_delta_wS","sin(delta) (with source holder)",200,-1,1); 
  //  analysisManager->CreateH1("sinxi_1xsindelta_wS","sin(xi1)*sin(delta) (with source holder)",200,-1,1); 
  analysisManager->CreateH1("delta_wS","delta (with source holder)",200,-TMath::Pi()/2,TMath::Pi()/2); 
  analysisManager->CreateH1("hAP1_wS","sin(2 theta) * cos(phi) (with source holder)", 200, -1,1); 
  //  analysisManager->CreateH2("position_xy","position_xy", 100, -100, 100, 100, -100, 100);
  //  analysisManager->CreateH2("position_zy","position_zy", 100, -100, 100, 100, -100, 100);
  //  analysisManager->CreateH2("position_zx","position_zx", 100, -100, 100, 100, -100, 100);

  //  analysisManager->CreateH2("phi_vs_energy","phi_vs_energy", 200, 0, TMath::Pi(), 200, 0, 1.5); 
  //  analysisManager->CreateH2("phi_vs_zVert","phi_vs_zVert", 200, 0, TMath::Pi(), 200, -2.5, 2.5); 
  //  analysisManager->CreateH2("phi_vs_zVert_atVertex","phi_vs_zVert_atVertex", 200, 0, TMath::Pi(), 200, -2.5, 2.5); 
  //  analysisManager->CreateH2("phi_vs_phiVert","phi_vs_phiVert", 200, 0, TMath::Pi(), 200, 0, TMath::Pi()); 


  //=================================================

  // (with APEX)
  analysisManager->CreateH1("hk1Spectrum_APEX","k1 Spectrum (APEX)", 220, 0., 1.1); 
  analysisManager->CreateH1("hk2Spectrum_APEX","k2 Spectrum (APEX)", 220, 0., 1.1); 
  analysisManager->CreateH1("hk3Spectrum_APEX","k3 Spectrum (APEX)", 220, 0., 1.1); 
  analysisManager->CreateH1("hSumSpectrum_APEX","Sum Spectrum (APEX)", 200, 0., 2.5); 
  //  analysisManager->CreateH1("hCosTheta_APEX","cos(theta) (APEX)", 200, -1, 1); 
  analysisManager->CreateH1("hTheta_APEX","theta (APEX)", 200, 0., TMath::Pi()); 
  analysisManager->CreateH1("hPsi_APEX","psi (APEX)", 100, 0., 2.*TMath::Pi()); 
  analysisManager->CreateH1("hPhi_APEX","phi (APEX)", 100, 0., TMath::Pi()); 
  //  analysisManager->CreateH1("hSin2Theta_APEX","sin(2theta) (APEX)", 200, -1,1); 
  //  analysisManager->CreateH1("hCosPhi_APEX","cos(phi) (APEX)", 200, -1,1); 
  //  analysisManager->CreateH1("hSinPsi_APEX","sin(psi) (APEX)", 200, -1,1); 
  analysisManager->CreateH1("hAP1_APEX","sin(2 theta) * cos(phi) (APEX)", 200, -1,1);
  
  // Histograms relevant to APEX reconstruction only
  analysisManager->CreateH1("totalEnergy_APEX","Total Energy (APEX)",5000,0,3000); 
  analysisManager->CreateH1("numBarsHit","Number of Bars Hit (APEX)",24,0,24); 
  analysisManager->CreateH1("hitbars", "APEX Bars Hit (APEX)",24,0,24); 
  analysisManager->CreateH1("X_mom", "X momentum (APEX)",1000,0,1); 
  analysisManager->CreateH1("Y_mom", "Y momentum (APEX)",1000,0,1); 
  analysisManager->CreateH1("Z_mom", "Z momentum (APEX)",1000,0,1); 
  analysisManager->CreateH1("hitTime","Hit Time (APEX)",1000,0,600); 
  analysisManager->CreateH1("totalEnergy_1Bar","Total Energy with 1 Bar Hit (APEX)",5000,0,3000); 
  analysisManager->CreateH1("totalEnergy_2Bars","Total Energy with 2 Bars Hit (APEX)",5000,0,3000);
  analysisManager->CreateH1("totalEnergy_3Bars","Total Energy with 3 Bars Hit (APEX)",5000,0,3000);
  analysisManager->CreateH1("totalEnergy_4Bars","Total Energy with 4 Bars Hit (APEX)",5000,0,3000);
  analysisManager->CreateH1("totalEnergy_5Bars","Total Energy with 5 Bars Hit (APEX)",5000,0,3000);
  analysisManager->CreateH1("totalEnergy_6Bars","Total Energy with 6 Bars Hit (APEX)",5000,0,3000);
  analysisManager->CreateH1("totalEnergy_7Bars","Total Energy with 7 Bars Hit (APEX)",5000,0,3000);
  analysisManager->CreateH1("totalEnergy_8Bars","Total Energy with 8 Bars Hit (APEX)",5000,0,3000);
  
  //Reyco wanted these (all bars)
  analysisManager->CreateH1("zPos_allHits","Z Position for All Hits (APEX)",1000,-275,275); 
  analysisManager->CreateH1("singles_eSpec","Single Hits Energy Spectrum (APEX)",1000,0,3000);
  
  //Reyco wanted these for 3 bar hits only
  analysisManager->CreateH1("singles_eSpec_3Bars","Single Hits Energy Spectrum (3 Bars Only) (APEX)",1000,0,3000); 
  analysisManager->CreateH1("zPos_allHits_3BarsOnly","Z Position for All Hits (3 Bars Only) (APEX)",1000,-275,275); 
  analysisManager->CreateH1("zPos_k1_3BarsOnly","Z Position of K1 (3 Bars Only) (APEX)",1000,-275,275); 
  analysisManager->CreateH1("zPos_k2_3BarsOnly","Z Position of K2 (3 Bars Only) (APEX)",1000,-275,275); 
  analysisManager->CreateH1("zPos_k3_3BarsOnly","Z Position of K3 (3 Bars Only) (APEX)",1000,-275,275);
  
  analysisManager->CreateH1("mag_vector_sum","Magnitude of Vector Sum (APEX)",200,-10.,10.); 
  analysisManager->CreateH1("hXi1_APEX","xi 1 (APEX)",200,-TMath::Pi()/2,TMath::Pi()/2); 
  analysisManager->CreateH1("hXi2_APEX","xi 2 (APEX)",200,-TMath::Pi()/2,TMath::Pi()/2); 
  analysisManager->CreateH1("hXi3_APEX","xi 3 (APEX)",200,-TMath::Pi()/2,TMath::Pi()/2);
  //  analysisManager->CreateH1("sinXi_1_APEX","sin(xi1) (APEX)",200,-1,1); 
  //  analysisManager->CreateH1("cosXi_2_APEX","cos(xi2) (APEX)",200,-1,1); 
  //  analysisManager->CreateH1("sin_delta_APEX","sin(delta) (APEX)",200,-1,1);
  //  analysisManager->CreateH1("sinxi_1xsindelta_APEX","sin(xi1)*sin(delta) (APEX)",200,-1,1); 
  analysisManager->CreateH1("delta_APEX","delta (APEX)",200,-TMath::Pi()/2,TMath::Pi()/2); 
  
  analysisManager->CreateH1("cosXi_rBP1","cos(xi)--relativeBarNumber: +1 (APEX)",200,-1,1); 
  analysisManager->CreateH1("cosXi_rBP2","cos(xi)--relativeBarNumber: +2 (APEX)",200,-1,1); 
  analysisManager->CreateH1("cosXi_rBP3","cos(xi)--relativeBarNumber: +3 (APEX)",200,-1,1); 
  analysisManager->CreateH1("cosXi_rBP4","cos(xi)--relativeBarNumber: +4 (APEX)",200,-1,1); 
  analysisManager->CreateH1("cosXi_rBM1","cos(xi)--relativeBarNumber: -1 (APEX)",200,-1,1); 
  analysisManager->CreateH1("cosXi_rBM2","cos(xi)--relativeBarNumber: -2 (APEX)",200,-1,1); 
  analysisManager->CreateH1("cosXi_rBM3","cos(xi)--relativeBarNumber: -3 (APEX)",200,-1,1); 
  analysisManager->CreateH1("cosXi_rBM4","cos(xi)--relativeBarNumber: -4 (APEX)",200,-1,1); 

  analysisManager->CreateH1("rh_asymmetry_Vertex","RH Asymmetry (Vertex)",200,-1,1); 
  analysisManager->CreateH1("rh_asymmetry_wS","RH Asymmetry (with source holder)",200,-1,1);
  analysisManager->CreateH1("rh_asymmetry_APEX","RH Asymmetry (APEX)",200,-1,1); 
  

  analysisManager->CreateH1("hk1Spectrum_APEX_Cut1","k1 Energy Spectrum (APEX Cut1)",200,0,1.1);
  analysisManager->CreateH1("hk2Spectrum_APEX_Cut1","k2 Energy Spectrum (APEX Cut1)",200,0,1.1);
  analysisManager->CreateH1("hk3Spectrum_APEX_Cut1","k3 Energy Spectrum (APEX Cut1)",200,0,1.1);
  analysisManager->CreateH1("hSumSpectrum_APEX_Cut1","Sum Energy Spectrum (APEX Cut1)",200,0,1.1);
  analysisManager->CreateH1("hXi1_APEX_Cut1","xi 1 (APEX Cut1)",200,-TMath::Pi()/2,TMath::Pi()/2);
  analysisManager->CreateH1("hXi2_APEX_Cut1","xi 2 (APEX Cut1)",200,-TMath::Pi()/2,TMath::Pi()/2);
  analysisManager->CreateH1("hXi3_APEX_Cut1","xi 3 (APEX Cut1)",200,-TMath::Pi()/2,TMath::Pi()/2); 



  //histos for collimator investigation 
  analysisManager->CreateH1("zPos_someHits","zpos from single bar within 511 range", 1000, -275, 275);
  //analysisManager->CreateH1("zPos_someHits2","zpos from total event energy dep sum",1000, -275, 275); 
  analysisManager->CreateH1("zPos_someHits3","zpos if two bars within 511 range and opposite(in APEX)",1000,-275,275);
  //analysisManager->CreateH1("all_bar_edep","all total Energy in APEX bars", 5000, 0., 3000);
  //analysisManager->CreateH1("bar_edep","total Energy in APEX bars within range", 5000, 0., 3000);
  //analysisManager->CreateH1("bar_edep2","Enery in APEX bars within range", 5000, 0., 3000); 

  
  //2D
  analysisManager->CreateH2("k1Energy_vs_k2Energy","k1Energy_vs_k2Energy", 500, 0, 1200, 500, 0, 1200);
  analysisManager->SetH2XAxisTitle(analysisManager->GetH2Id("k1Energy_vs_k2Energy",1),"k1 Energy");
  analysisManager->SetH2YAxisTitle(analysisManager->GetH2Id("k1Energy_vs_k2Energy",1),"k2 Energy");

  analysisManager->CreateH2("theta_v_phi_Vertex","Theta vs Phi", 200, 0, TMath::Pi(), 200, 0, TMath::Pi()); 
  analysisManager->SetH2XAxisTitle(analysisManager->GetH2Id("theta_v_phi_Vertex",1),"Theta");
  analysisManager->SetH2YAxisTitle(analysisManager->GetH2Id("theta_v_phi_Vertex",1),"Phi");
  analysisManager->CreateH2("theta_v_psi_Vertex","Theta vs Psi", 200, 0, TMath::Pi(), 200, 0, TMath::Pi()); 
  analysisManager->SetH2XAxisTitle(analysisManager->GetH2Id("theta_v_psi_Vertex",1),"Theta");
  analysisManager->SetH2YAxisTitle(analysisManager->GetH2Id("theta_v_psi_Vertex",1),"Psi");
  //Psi is restricted due to kinematics
  analysisManager->CreateH2("phi_v_psi_Vertex","Phi vs Psi", 200, 0, TMath::Pi(), 200, 0, TMath::Pi()); 
  analysisManager->SetH2XAxisTitle(analysisManager->GetH2Id("phi_v_psi_Vertex",1),"Phi");
  analysisManager->SetH2YAxisTitle(analysisManager->GetH2Id("phi_v_psi_Vertex",1),"Psi");

  analysisManager->CreateH2("theta_v_phi_wS","Theta vs Phi", 200, 0, TMath::Pi(), 200, 0, TMath::Pi()); 
  analysisManager->SetH2XAxisTitle(analysisManager->GetH2Id("theta_v_phi_wS",1),"Theta");
  analysisManager->SetH2YAxisTitle(analysisManager->GetH2Id("theta_v_phi_wS",1),"Phi");
  analysisManager->CreateH2("theta_v_psi_wS","Theta vs Psi", 200, 0, TMath::Pi(), 200, 0, TMath::Pi()); 
  analysisManager->SetH2XAxisTitle(analysisManager->GetH2Id("theta_v_psi_wS",1),"Theta");
  analysisManager->SetH2YAxisTitle(analysisManager->GetH2Id("theta_v_psi_wS",1),"Psi");
  //Psi is restricted due to kinematics
  analysisManager->CreateH2("phi_v_psi_wS","Phi vs Psi", 200, 0, TMath::Pi(), 200, 0, TMath::Pi()); 
  analysisManager->SetH2XAxisTitle(analysisManager->GetH2Id("phi_v_psi_wS",1),"Phi");
  analysisManager->SetH2YAxisTitle(analysisManager->GetH2Id("phi_v_psi_wS",1),"Psi");

  analysisManager->CreateH2("theta_v_phi_APEX","Theta vs Phi", 200, 0, TMath::Pi(), 200, 0, TMath::Pi()); 
  analysisManager->SetH2XAxisTitle(analysisManager->GetH2Id("theta_v_phi_APEX",1),"Theta");
  analysisManager->SetH2YAxisTitle(analysisManager->GetH2Id("theta_v_phi_APEX",1),"Phi");
  analysisManager->CreateH2("theta_v_psi_APEX","Theta vs Psi", 200, 0, TMath::Pi(), 200, 0, TMath::Pi()); 
  analysisManager->SetH2XAxisTitle(analysisManager->GetH2Id("theta_v_psi_APEX",1),"Theta");
  analysisManager->SetH2YAxisTitle(analysisManager->GetH2Id("theta_v_psi_APEX",1),"Psi");
  //Psi is restricted due to kinematics
  analysisManager->CreateH2("phi_v_psi_APEX","Phi vs Psi", 200, 0, TMath::Pi(), 200, 0, TMath::Pi()); 
  analysisManager->SetH2XAxisTitle(analysisManager->GetH2Id("phi_v_psi_APEX",1),"Phi");
  analysisManager->SetH2YAxisTitle(analysisManager->GetH2Id("phi_v_psi_APEX",1),"Psi");
  
  analysisManager->CreateH2("energy_v_numbars","Energy v Number of Bars", 10, 0, 10, 1000, 0, 2000); 
  analysisManager->SetH2XAxisTitle(analysisManager->GetH2Id("energy_v_numbars",1),"Energy");
  analysisManager->SetH2YAxisTitle(analysisManager->GetH2Id("energy_v_numbars",1),"Number of Bars");
  analysisManager->CreateH2("k1zPosition_v_k2zPosition","k1zPosition_vs_k2zPosition",1000,-275,275,1000,-275,275); 
  analysisManager->SetH2XAxisTitle(analysisManager->GetH2Id("k1zPosition_v_k2zPosition",1),"k1 Z");
  analysisManager->SetH2YAxisTitle(analysisManager->GetH2Id("k1zPosition_v_k2zPosition",1),"k2 Z");


  
}


RunAction::~RunAction()
{
  delete G4AnalysisManager::Instance();
}

void RunAction::UpdateRunAction()
{
  physListUpdate=G4RunManager::GetRunManager()->GetUserPhysicsList();
  myPhysListUpdate=dynamic_cast<PhysicsList*>(const_cast<G4VUserPhysicsList*>(physListUpdate));
  BField=myPhysListUpdate->GetBFieldValue();

  // Set up timing histograms
  //Set B Field
  //  timingHist = timedistr->buildHist("oPSDecay",5000);
  timingHistoPS = timedistr->buildHist("oPSDecay",5000);
  timingHistpPS = timedistr->buildHist("pPSDecay",5000);
  BRatio=timedistr->GetBranchingRatio("oPSDecay",5000);
  BRatio=timedistr->GetBranchingRatio("pPSDecay",5000);
}
void RunAction::SetOutputFile(G4String filePath)
{
  myPath=filePath;

}

void RunAction::BeginOfRunAction(const G4Run* run)
{
  //  G4long seed=time(0);
  //CLHEP::HepRandom::setTheSeed(seed);

  
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);

  // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
  //  std::string
  G4String fFileName = AnalysisManager::getInstance()->GetFileName();

  std::string pathString = myPath;
  //  std::string filename = pathString.append(fFileName);
  G4String filename = myPath+fFileName;//myPath.append(fFileName);
  // Open an output file
  analysisManager->OpenFile(filename.c_str());

  //  rootanalysistruth->BeginOfRunAction();
  //  rootanalysis->BeginOfRunAction();
}

void RunAction::EndOfRunAction(const G4Run*)
{ 
  /*G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  G4Polyline x_axis;
  G4Polyline y_axis;
  G4Polyline z_axis;
  //Set red line color                                                                                    
  G4Colour red(1.0,0.0,0.0);
  G4VisAttributes att(red);
  x_axis.push_back( G4Point3D (0.,0.,0.));
  x_axis.push_back( G4Point3D (0.*cm,500.*cm,0.));
  x_axis.SetVisAttributes(att);

  y_axis.push_back( G4Point3D (0.,0.,0.));
  y_axis.push_back( G4Point3D (500.*cm,0.*cm,0.));
  y_axis.SetVisAttributes(att);

  z_axis.push_back( G4Point3D (0.,0.,0.));
  z_axis.push_back( G4Point3D (0.*cm,0.*cm,500.*cm));
  z_axis.SetVisAttributes(att);

  //Get the pointer to the User Interface manager                                                         
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  //Let's make it interactive mode for now                                                                
  //  G4UIExecutive* ui = new G4UIExecutive(argc, argv);

  UImanager->ApplyCommand("/vis~/clear/view");
  UImanager->ApplyCommand("/vis~/draw/current");
  pVVisManager->Draw(x_axis);
  pVVisManager->Draw(y_axis);
  pVVisManager->Draw(z_axis);
  UImanager->ApplyCommand(".vis~/show/view");*/



  // Save histograms
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  analysisManager->Write();
  analysisManager->CloseFile();
  //  rootanalysistruth->EndOfRunAction();
  //  rootanalysis->EndOfRunAction();
}

TH1F* RunAction::GetoPSTimingHistogramPtr(){
  return timingHistoPS;
}
TH1F* RunAction::GetpPSTimingHistogramPtr(){
  return timingHistpPS;
}
TH1F* RunAction::GetoPS_2GammaTimingHistogramPtr(){
  return timingHistoPS_2Gamma;
}
TH1F* RunAction::GetpPS_2GammaTimingHistogramPtr(){
  return timingHistpPS_2Gamma;
}
TH1F* RunAction::GetoPS_3GammaTimingHistogramPtr(){
  return timingHistoPS_3Gamma;
}
TH1F* RunAction::GetpPS_3GammaTimingHistogramPtr(){
  return timingHistpPS_3Gamma;
}



Double_t RunAction::GetBranchingRatio(){
  return BRatio;
}
