///--------------------------------------------------------------------
// EventAction.cc
//
// Description: Events!
// Changes: 7/18/05 None yet
//--------------------------------------------------------------------
#include <boost/bind.hpp>
#include <iostream>
#include "CutImplementations.hh"
#include "UserEventInformation.hh"
#include "G4EventManager.hh"
#include "G4RunManager.hh"
#include <TF2.h>
#include <TF1.h>
#include <TH2.h>
#include <TH1D.h>
#include <TFile.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TVector3.h>
#include <TROOT.h>
#include "G4Event.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4GeneralParticleSource.hh"
#include <vector>
#include <iostream>
#include <iterator>
#include <algorithm>
#include "Randomize.hh"
#include "AnalysisManager.hh"
#include "EventAction.hh"
#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include <cmath>
#include "TimeDistribution.hh"

EventAction::EventAction(RunAction* run,PrimaryGeneratorAction* gen)
  :runaction(run),genaction(gen)
{
  // Initialize o-Ps and p-Ps timing histograms and branching ratios
  timingHistoPS = run->GetoPSTimingHistogramPtr();
  timingHistpPS = run->GetpPSTimingHistogramPtr();
  timingHistoPS_2Gamma = run->GetoPS_2GammaTimingHistogramPtr();
  timingHistpPS_2Gamma = run->GetpPS_2GammaTimingHistogramPtr();
  timingHistoPS_3Gamma = run->GetoPS_3GammaTimingHistogramPtr();
  timingHistpPS_3Gamma = run->GetpPS_3GammaTimingHistogramPtr();
  BRatio = run->GetBranchingRatio();
  //  N = 0;  
  std::string fFileName = AnalysisManager::getInstance()->GetFileName();
  std::string filename = "/nas02/home/c/b/cbartram/Caliope_Work/CALIOPE-BUILD/events";
  eventFile.open(filename.c_str(),std::fstream::app);

}

EventAction::~EventAction()
{ 
  gammaRayStorage_withSource.clear();
  gammaRayStorage_noSource.clear();
  gammaRayCoordinates.clear();
}

void EventAction::BeginOfEventAction(const G4Event* evt)
{
  pSSubtype=10;
  
  int eventID=G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
  if (eventID%100==0) {
    G4cout << eventID << G4endl;
  }

  G4EventManager::GetEventManager()->SetUserInformation(new UserEventInformation);
  UserEventInformation *eventinfo=(UserEventInformation*)evt->GetUserInformation();
  eventinfo->SetCopyNo(-2);

  //Initialize an array to contain information pertaining to the highest and second highest energy gamma rays

  for (int i = 0; i < 3; i++) {
    std::vector<double> row_withSource; // Create an empty row
    for (int j = 0; j < 5; j++) {
      row_withSource.push_back(-100); // Add an element (column) to the row
    }
    gammaRayStorage_withSource.push_back(row_withSource); // Add the row to the main vector
  }
  for (int i = 0; i < 3; i++) {
    std::vector<double> row_noSource; // Create an empty row
    for (int j = 0; j < 5; j++) {
      row_noSource.push_back(-100); // Add an element (column) to the row
    }
    gammaRayStorage_noSource.push_back(row_noSource); // Add the row to the main vector
  }
  for (int i = 0; i < 5; i++) {
    std::vector<double> row_Coordinates; // Create an empty row
    for (int j = 0; j < 5; j++) {
      row_Coordinates.push_back(-100); // Add an element (column) to the row
    }
    gammaRayCoordinates.push_back(row_Coordinates); // Add the row to the main vector
  }

  event_id = evt->GetEventID();
}

void EventAction::EndOfEventAction(const G4Event* evt)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  std::vector< std::vector<double> >* gammastorage_Coordinates = GetGammaRayCoordinates();
  std::vector< std::vector<double> >* gammastorage_noSource = GetGammaRayVector_NoSource();
  //  G4double* PSSubType;
  //  *PSSubType=-1;
  G4double* PSSubType = GetPsSubTypePtr();
  
  //  if(PSSubType!=NULL){
  IDNum = analysisManager->GetH1Id("PsSubType",1);
  if(*PSSubType!=10){
    analysisManager->FillH1(IDNum,*PSSubType);
  }
      //}
      //}

  
  AlternativeReconstruction(gammastorage_noSource, gammastorage_Coordinates, "noSource");
  std::vector< std::vector<double> >* gammastorage_withSource = GetGammaRayVector_WithSource();
  AlternativeReconstruction(gammastorage_withSource, gammastorage_Coordinates, "withSource");
  APEXReconstruction();  
  SetGammaRayVector_NoSource(gammastorage_noSource);
  SetGammaRayVector_WithSource(gammastorage_withSource);
  //Clear Gamma Ray Vectors
  gammaRayStorage_withSource.clear();
  gammaRayStorage_noSource.clear();
  gammaRayCoordinates.clear();
}

std::vector< std::vector<double> >* EventAction::GetGammaRayVector_NoSource()
{
  return &gammaRayStorage_noSource;
}

std::vector< std::vector<double> >* EventAction::GetGammaRayVector_WithSource()
{
  return &gammaRayStorage_withSource;
}

std::vector< std::vector<double> >* EventAction::GetGammaRayCoordinates()
{
  return &gammaRayCoordinates;
}

G4double* EventAction::GetPsSubTypePtr()
{
  return &pSSubtype;
}
void EventAction::SetPsSubTypePtr(G4double* subtype)
{
  //Do nothing for now
  //  pSSubtype=subtype;
}

void EventAction::SetGammaRayVector(int val)
{
  num=val;
}

void EventAction::SetGammaRayVector_NoSource(std::vector< std::vector<double> >* gammaRayData)
{
  (*gammaRayData)[0][0]=-100;
  (*gammaRayData)[0][1]=-100;
  (*gammaRayData)[0][2]=-100;
  (*gammaRayData)[0][3]=-100;

  (*gammaRayData)[1][0]=-100;
  (*gammaRayData)[1][1]=-100;
  (*gammaRayData)[1][2]=-100;
  (*gammaRayData)[1][3]=-100;
}

void EventAction::SetGammaRayVector_WithSource(std::vector< std::vector<double> >* gammaRayData)
{
  (*gammaRayData)[0][0]=-100;
  (*gammaRayData)[0][1]=-100;
  (*gammaRayData)[0][2]=-100;
  (*gammaRayData)[0][3]=-100;

  (*gammaRayData)[1][0]=-100;
  (*gammaRayData)[1][1]=-100;
  (*gammaRayData)[1][2]=-100;
  (*gammaRayData)[1][3]=-100;
}

void EventAction::AlternativeReconstruction(std::vector< std::vector<double> >* gammaRayData, std::vector< std::vector<double> >* positions, std::string reconType)//,int Num)
{
  TVector3 k1(-100,-100,-100);
  TVector3 k2(-100,-100,-100);
  TVector3 k3(-100,-100,-100);

  //=================================
  // Check the size of the vector
  //=================================

  int vecSize=(*gammaRayData).size();
  
  std::vector<double> energyVector;
  std::vector<double> energyVectorDuplicate;
  for(int i=0;i<(*gammaRayData).size();i++)
    {
      energyVector.push_back((*gammaRayData)[i][0]);
    }
  energyVectorDuplicate=energyVector;
  std::sort(energyVector.begin(),energyVector.end());
  std::reverse(energyVector.begin(),energyVector.end());

  double highest=energyVector[0];
  double secondhighest=energyVector[1];
  double thirdhighest=energyVector[2];
  std::vector<double>::iterator it;

  it=std::find(energyVectorDuplicate.begin(),energyVectorDuplicate.end(),highest);
  int index=std::distance(energyVectorDuplicate.begin(), it);

  double HighestKin=(*gammaRayData)[index][0];

  double kin_energy_1=(*gammaRayData)[index][0];
  dir_x_1=(*gammaRayData)[index][1];
  dir_y_1=(*gammaRayData)[index][2];
  dir_z_1=(*gammaRayData)[index][3];

  it=std::find(energyVectorDuplicate.begin(),energyVectorDuplicate.end(),secondhighest);
  index=std::distance(energyVectorDuplicate.begin(), it);

  double SecondHighestKin=(*gammaRayData)[index][0];

  double kin_energy_2=(*gammaRayData)[index][0];
  dir_x_2=(*gammaRayData)[index][1];
  dir_y_2=(*gammaRayData)[index][2];
  dir_z_2=(*gammaRayData)[index][3];

  it=std::find(energyVectorDuplicate.begin(),energyVectorDuplicate.end(),thirdhighest);
  index=std::distance(energyVectorDuplicate.begin(), it);

  double ThirdHighestKin=(*gammaRayData)[index][0];

  double kin_energy_3=(*gammaRayData)[index][0];
  dir_x_3=(*gammaRayData)[index][1];
  dir_y_3=(*gammaRayData)[index][2];
  dir_z_3=(*gammaRayData)[index][3];

  energyVector.clear();
  energyVectorDuplicate.clear();

  double pos_x;
  double pos_y;
  double pos_z;
 
  double k1_energy=kin_energy_1;
  double k2_energy=kin_energy_2;
  double k3_energy=kin_energy_3;

  if (((dir_x_1)!=-100)&&((dir_y_1)!=-100)&&((dir_z_1)!=-100)&&((dir_x_2)!=-100)&&((dir_y_2)!=-100)&&((dir_z_2)!=-100)) {
    k1.SetXYZ(dir_x_1,dir_y_1,dir_z_1);
    k2.SetXYZ(dir_x_2,dir_y_2,dir_z_2);
    k3.SetXYZ(dir_x_3,dir_y_3,dir_z_3);
    k1.Unit();
    k2.Unit();
    k3.Unit();
    //FIX THIS
    if ((k1.X()==k2.X())&&(k1.Y()==k2.Y())&&(k1.Z()==k2.Z()))
      {}
    else{
      //=========================
      // Yamazaki Coordinates
      //=========================
      
    TVector3 normal;
    double theta;
    double psi;
    double ksi;
    TVector3 sprojection;
    double phi=-5000;	
    double x1;
    double y1;
    double z1;
    double x2;
    double y2;
    double z2;

    //=========================
    // Ryan Coordinates
    //=========================
    double xi1, xi2, xi3;
    double delta;
    double k1_mag, k2_mag, k3_mag;
    double sinxi1, cosxi1, sinxi2, cosxi2, sinxi3, cosxi3, sindelta;
    
    k1_mag = k1.Mag(); k2_mag = k2.Mag(); k3_mag = k3.Mag();
    sinxi1 = k1.Z()/k1_mag;
    xi1 = TMath::ASin(sinxi1);
    cosxi1 = TMath::Sqrt(1-k1.Z()*k1.Z()/k1_mag/k1_mag);
    sinxi2 = k2.Z()/k2_mag;
    cosxi2 = TMath::Sqrt(1-k2.Z()*k2.Z()/k2_mag/k2_mag);
    xi2 = TMath::ASin(sinxi2);
    cosxi3 = TMath::Sqrt(1-k3.Z()*k3.Z()/k3_mag/k3_mag);
    sinxi3 = k3.Z()/k3_mag;
    xi3 = TMath::ASin(sinxi3);
    sindelta = (k1.X()*k2.Y()-k1.Y()*k2.X())/cosxi1/cosxi2/k1_mag/k2_mag;
    delta = TMath::ASin(sindelta);

    if ((*positions).size()!=0) {
      pos_x=(*positions)[0][0];
      pos_y=(*positions)[0][1];
      pos_z=(*positions)[0][2];
    }    
    
    //===========================
    // Reconstruction
    //===========================
    
    normal = k1.Cross(k2);
    theta = TMath::ACos(normal.CosTheta());
    psi = k1.Angle(k2);
    sprojection = normal;
    sprojection.Rotate(TMath::Pi()*0.5, normal.Cross(TVector3(0,0,1))); // projection of z-axis (S) onto decay plane.
    phi = k1.Angle(sprojection);
    
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    if (reconType.compare("noSource")==0) {
      IDNum = analysisManager->GetH1Id("hk1Spectrum",1);
      analysisManager->FillH1(IDNum,k1_energy);
      IDNum = analysisManager->GetH1Id("hk2Spectrum",1);
      analysisManager->FillH1(IDNum,k2_energy);
      IDNum = analysisManager->GetH1Id("hk3Spectrum",1);
      analysisManager->FillH1(IDNum,k3_energy);
      IDNum = analysisManager->GetH1Id("hSumSpectrum",1);
      analysisManager->FillH1(IDNum,k1_energy+k2_energy+k3_energy);
      //      IDNum = analysisManager->GetH1Id("hCosTheta",1);
      //      analysisManager->FillH1(IDNum,k1.CosTheta());
      IDNum = analysisManager->GetH1Id("hTheta",1);
      analysisManager->FillH1(IDNum,theta);
      IDNum = analysisManager->GetH1Id("hPsi",1);
      analysisManager->FillH1(IDNum,psi);
      IDNum = analysisManager->GetH1Id("hPhi",1);
      analysisManager->FillH1(IDNum,phi);
      //      IDNum = analysisManager->GetH1Id("hSin2Theta",1);
      //      analysisManager->FillH1(IDNum,TMath::Sin(2.*theta));
      //      IDNum = analysisManager->GetH1Id("hCosPhi",1);
      //      analysisManager->FillH1(IDNum,TMath::Cos(phi));
      //      IDNum = analysisManager->GetH1Id("hSinPsi",1);
      //      analysisManager->FillH1(IDNum,TMath::Sin(psi));
      IDNum = analysisManager->GetH1Id("hXi1",1);
      analysisManager->FillH1(IDNum,xi1);
      IDNum = analysisManager->GetH1Id("hXi2",1);
      analysisManager->FillH1(IDNum,xi2);
      IDNum = analysisManager->GetH1Id("hXi3",1);
      analysisManager->FillH1(IDNum,xi3);
      //      IDNum = analysisManager->GetH1Id("sin_delta",1);
      //      analysisManager->FillH1(IDNum,sindelta);
      IDNum = analysisManager->GetH1Id("delta",1);
      analysisManager->FillH1(IDNum,delta);
      //      IDNum = analysisManager->GetH1Id("sinXi_1",1);
      //      analysisManager->FillH1(IDNum,TMath::Sin(xi1));
      //      IDNum = analysisManager->GetH1Id("cosXi_2",1);
      //      analysisManager->FillH1(IDNum,TMath::Cos(xi2));
      //      IDNum = analysisManager->GetH1Id("sinxi_1xsindelta",1);
      //      analysisManager->FillH1(IDNum,TMath::Sin(xi1)*TMath::Sin(delta));
      IDNum = analysisManager->GetH1Id("rh_asymmetry_Vertex",1);
      analysisManager->FillH1(IDNum,TMath::Sin(xi1)*TMath::Cos(xi1)*TMath::Cos(xi1)*TMath::Sin(delta));
      IDNum = analysisManager->GetH2Id("theta_v_phi_Vertex",1);
      analysisManager->FillH2(IDNum,theta,phi);
      IDNum = analysisManager->GetH2Id("theta_v_psi_Vertex",1);
      analysisManager->FillH2(IDNum,theta,psi);
      IDNum = analysisManager->GetH2Id("phi_v_psi_Vertex",1);
      analysisManager->FillH2(IDNum,phi,psi);
      IDNum = analysisManager->GetH1Id("hAP1",1);
      analysisManager->FillH1(IDNum,TMath::Cos(phi) * TMath::Sin(2.*theta));
    }
        
    if (reconType.compare("withSource")==0) {
      IDNum = analysisManager->GetH1Id("hk1Spectrum_wS",1);
      analysisManager->FillH1(IDNum,k1_energy);
      IDNum = analysisManager->GetH1Id("hk2Spectrum_wS",1);
      analysisManager->FillH1(IDNum,k2_energy);
      IDNum = analysisManager->GetH1Id("hk3Spectrum_wS",1);
      analysisManager->FillH1(IDNum,k3_energy);
      IDNum = analysisManager->GetH1Id("hSumSpectrum_wS",1);
      analysisManager->FillH1(IDNum,k1_energy+k2_energy+k3_energy);
      //      IDNum = analysisManager->GetH1Id("hCosTheta_wS",1);
      //      analysisManager->FillH1(IDNum,k1.CosTheta());
      IDNum = analysisManager->GetH1Id("hTheta_wS",1);
      analysisManager->FillH1(IDNum,theta);
      IDNum = analysisManager->GetH1Id("hPsi_wS",1);
      analysisManager->FillH1(IDNum,psi);
      IDNum = analysisManager->GetH1Id("hPhi_wS",1);
      analysisManager->FillH1(IDNum,phi);
      //      IDNum = analysisManager->GetH1Id("hSin2Theta_wS",1);
      //      analysisManager->FillH1(IDNum,TMath::Sin(2.*theta));
      //      IDNum = analysisManager->GetH1Id("hCosPhi_wS",1);
      //      analysisManager->FillH1(IDNum,TMath::Cos(phi));
      //      IDNum = analysisManager->GetH1Id("hSinPsi_wS",1);
      //      analysisManager->FillH1(IDNum,TMath::Sin(psi));
      IDNum = analysisManager->GetH1Id("hXi1_wS",1);
      analysisManager->FillH1(IDNum,xi1);
      IDNum = analysisManager->GetH1Id("hXi2_wS",1);
      analysisManager->FillH1(IDNum,xi2);
      IDNum = analysisManager->GetH1Id("hXi3_wS",1);
      analysisManager->FillH1(IDNum,xi3);
      //      IDNum = analysisManager->GetH1Id("sin_delta_wS",1);
      //      analysisManager->FillH1(IDNum,sindelta);
      IDNum = analysisManager->GetH1Id("delta_wS",1);
      analysisManager->FillH1(IDNum,delta);
      //      IDNum = analysisManager->GetH1Id("sinXi_1_wS",1);
      //      analysisManager->FillH1(IDNum,TMath::Sin(xi1));
      //      IDNum = analysisManager->GetH1Id("cosXi_2_wS",1);
      //      analysisManager->FillH1(IDNum,TMath::Cos(xi2));
      //      IDNum = analysisManager->GetH1Id("sinxi_1xsindelta_wS",1);
      //      analysisManager->FillH1(IDNum,TMath::Sin(xi1)*TMath::Sin(delta));
      IDNum = analysisManager->GetH1Id("rh_asymmetry_wS",1);
      analysisManager->FillH1(IDNum,TMath::Sin(xi1)*TMath::Cos(xi1)*TMath::Cos(xi1)*TMath::Sin(delta));
      IDNum = analysisManager->GetH2Id("theta_v_phi_wS",1);
      analysisManager->FillH2(IDNum,theta,phi);
      IDNum = analysisManager->GetH2Id("theta_v_psi_wS",1);
      analysisManager->FillH2(IDNum,theta,psi);
      IDNum = analysisManager->GetH2Id("phi_v_psi_wS",1);
      analysisManager->FillH2(IDNum,phi,psi);
      IDNum = analysisManager->GetH1Id("hAP1_wS",1);
      analysisManager->FillH1(IDNum,TMath::Cos(phi) * TMath::Sin(2.*theta));
    }
      }
  }
}

void EventAction::APEXReconstruction()
{
  CutImplementations reconstruction;
  
  // Get the Analysis Manager for writing to histograms
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  // Get UserEventInfo to access data
  UserEventInformation *eventinfo=(UserEventInformation*)G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetUserInformation();

  //Initialize variables
  G4int copyNum=-1;
  std::set<G4int> hitBarsSet=eventinfo->GetSetOfHitBars();
  TRandom3 gRandom(0);
  G4double randEdep=0;  
  std::vector<std::pair<G4int,G4double> > gammaRayVector;
  std::vector<G4int> threeBars;
  std::vector<G4double> k1GammaRayData;
  std::vector<G4double> k2GammaRayData;
  std::vector<G4double> k3GammaRayData;
  G4double energy=0;
  G4double angle=0;
  G4double k1zRecon, k1xRecon, k1yRecon;
  k1zRecon=k1xRecon=k1yRecon=0;
  G4double k2zRecon, k2xRecon, k2yRecon;
  k2zRecon=k2xRecon=k2yRecon=0;
  G4double k3zRecon, k3xRecon, k3yRecon;
  k3zRecon=k3xRecon=k3yRecon=0;
  G4double k1energy, k2energy, k3energy;
  k1energy=k2energy=k3energy=0;
  G4double kzRecon=0;
  G4double k1hitTime=0;
  G4double k2hitTime=0;
  G4double k3hitTime=0;
  int barCount=0;
  int bar1=-1;
  int bar2=-1;
  int bar3=-1;
  G4double totalEnergy=0;
  G4double hitTime=0;
  G4int numStepsInBar=0;
  std::set<G4int>::iterator it;
  G4double k1TotZ, k2TotZ, k3TotZ;
  k1TotZ=k2TotZ=k3TotZ=0;
  G4double kTotZ=0;
  G4double  k1NumSteps, k2NumSteps, k3NumSteps;
  k1NumSteps=k2NumSteps=k3NumSteps=0;
  G4double kNumSteps=0;
  G4double pos_z=0;
  TF1* cauchyFunction;
  G4double phi1, theta1, phi2, theta2, phi3, theta3;
  G4double k1XMomentum, k1YMomentum, k1ZMomentum, k1Momentum;
  G4double k2XMomentum, k2YMomentum, k2ZMomentum, k2Momentum;
  G4double k3XMomentum, k3YMomentum, k3ZMomentum, k3Momentum;
  G4double kXMomentum, kYMomentum, kZMomentum;
  G4double kXMomentumAbs, kYMomentumAbs, kZMomentumAbs;
  
  if((hitBarsSet.size()<24)&&(hitBarsSet.size()>0)){
    for (it=hitBarsSet.begin(); it!=hitBarsSet.end(); it++){
      copyNum=*it;
      energy=eventinfo->GetEnergyDeposition(copyNum);
      randEdep = gRandom.Gaus(energy, 0.14*energy/2.355);
      totalEnergy+=randEdep;
	
	//COLLIMATOR--511 SINGLE BAR EDEP HISTO--START
	//if single bar edep is within 511 range, it fills the hit bar's position
	if(((randEdep*1000)<611)&&((randEdep*1000)>411)){
		kTotZ=(eventinfo->GetTotZ(copyNum));
      		kNumSteps=eventinfo->GetNumStepsInBar(copyNum);
      		pos_z=kTotZ/kNumSteps;
      		cauchyFunction = new TF1("cauchyFn","TMath::CauchyDist(x,[0],[1])",-275,275);
      		cauchyFunction->SetParameters(pos_z,35);
      		kzRecon = cauchyFunction->GetRandom();
      		delete cauchyFunction;
      		IDNum = analysisManager->GetH1Id("zPos_someHits",1);
      		analysisManager->FillH1(IDNum,kzRecon);	
	}
	//END

      gammaRayVector.push_back(std::make_pair(*it,randEdep));
      IDNum = analysisManager->GetH1Id("hitbars",1);
      analysisManager->FillH1(IDNum,copyNum);      
      /*      hitTime=eventinfo->GetHitTime(copyNum);
      analysisManager->FillH1(43,hitTime);
      numStepsInBar=eventinfo->GetNumStepsInBar(copyNum);
      hitTime=hitTime/numStepsInBar;     */
      kTotZ=(eventinfo->GetTotZ(copyNum));
      kNumSteps=eventinfo->GetNumStepsInBar(copyNum);
      pos_z=kTotZ/kNumSteps;
      cauchyFunction = new TF1("cauchyFn","TMath::CauchyDist(x,[0],[1])",-275,275);
      cauchyFunction->SetParameters(pos_z,35);
      kzRecon = cauchyFunction->GetRandom();
      //      kzRecon = kzRecon;
      delete cauchyFunction;
      IDNum = analysisManager->GetH1Id("zPos_allHits",1);
      analysisManager->FillH1(IDNum,kzRecon);
      IDNum = analysisManager->GetH1Id("singles_eSpec",1);
      analysisManager->FillH1(IDNum,randEdep*1000);
    }

    if(hitBarsSet.size()==1){
      IDNum = analysisManager->GetH1Id("totalEnergy_1Bar",1);
      analysisManager->FillH1(IDNum,totalEnergy*1000);
    }
    if(hitBarsSet.size()==2){
      IDNum = analysisManager->GetH1Id("totalEnergy_2Bars",1);
      analysisManager->FillH1(IDNum,totalEnergy*1000);

	//COLLIMATOR--FILL IF OPPOSITE BARS W/IN 511 RANGE
	//two bars were hit, if both are within 511 range and the abs. value of 
	//the difference of their copynums is 12, it fills both their positions
	//not very efficient due to 2 hit bar restriction
		G4int barnum[2];
		barnum[0]=25;
		barnum[1]=25;
		for(it=hitBarsSet.begin(); it!=hitBarsSet.end(); it++){
		copyNum=*it;	
		energy=eventinfo->GetEnergyDeposition(copyNum);
      		randEdep = gRandom.Gaus(energy, 0.14*energy/2.355);
			if((randEdep*1000 > 411) && (randEdep*1000 < 611)){
				if(barnum[0]==25) barnum[0]=copyNum;
				else barnum[1]=copyNum;
			}
		}
		if(abs(barnum[0]-barnum[1])==12){
		for(it=hitBarsSet.begin(); it!=hitBarsSet.end(); it++){
		copyNum=*it;
		kTotZ=(eventinfo->GetTotZ(copyNum));
      		kNumSteps=eventinfo->GetNumStepsInBar(copyNum);
      		pos_z=kTotZ/kNumSteps;
		cauchyFunction = new TF1("cauchyFn","TMath::CauchyDist(x,[0],[1])",-275,275);
      		cauchyFunction->SetParameters(pos_z,35);
      		kzRecon = cauchyFunction->GetRandom();
      		delete cauchyFunction;
      		IDNum = analysisManager->GetH1Id("zPos_someHits3",1);
		analysisManager->FillH1(IDNum,kzRecon);
		}}
		
	//END

    }
    if(hitBarsSet.size()==3){
      IDNum = analysisManager->GetH1Id("totalEnergy_3Bars",1);
      analysisManager->FillH1(IDNum,totalEnergy*1000);
    }
    if(hitBarsSet.size()==4){
      IDNum = analysisManager->GetH1Id("totalEnergy_4Bars",1);
      analysisManager->FillH1(IDNum,totalEnergy*1000);
    }
    if(hitBarsSet.size()==5){
      IDNum = analysisManager->GetH1Id("totalEnergy_5Bars",1);
      analysisManager->FillH1(IDNum,totalEnergy*1000);
    }
    if(hitBarsSet.size()==6){
      IDNum = analysisManager->GetH1Id("totalEnergy_6Bars",1);
      analysisManager->FillH1(IDNum,totalEnergy*1000);
    }
    if(hitBarsSet.size()==7){
      IDNum = analysisManager->GetH1Id("totalEnergy_7Bars",1);
      analysisManager->FillH1(IDNum,totalEnergy*1000);
    }
    if(hitBarsSet.size()==8){
      IDNum = analysisManager->GetH1Id("totalEnergy_8Bars",1);
      analysisManager->FillH1(IDNum,totalEnergy*1000);
    }
    
    IDNum = analysisManager->GetH2Id("energy_v_numbars",1);
    analysisManager->FillH2(IDNum,hitBarsSet.size(),totalEnergy*1000);
    //Sorts the bars from highest energy deposition to lowest
    std::sort(gammaRayVector.begin(), gammaRayVector.end(),
	      boost::bind(&std::pair<G4int, G4double>::second, _1) >
	      boost::bind(&std::pair<G4int, G4double>::second, _2));
    totalEnergy=totalEnergy*1000;
    IDNum = analysisManager->GetH1Id("totalEnergy_APEX",1);
    analysisManager->FillH1(IDNum,totalEnergy);

    for(int num=0;num<gammaRayVector.size();num++){

      // Get Important Information
      copyNum=gammaRayVector[num].first;
      energy=gammaRayVector[num].second;
      hitTime=eventinfo->GetHitTime(copyNum);
      kTotZ=eventinfo->GetTotZ(copyNum);
      kNumSteps=eventinfo->GetNumStepsInBar(copyNum);
      pos_z=kTotZ/kNumSteps;
      cauchyFunction = new TF1("cauchyFn","TMath::CauchyDist(x,[0],[1])",-275,275);
      cauchyFunction->SetParameters(pos_z,35);
      kzRecon = cauchyFunction->GetRandom();
      kzRecon = kzRecon/275.;
      delete cauchyFunction;
      
      if (barCount==1) {
	k1hitTime=hitTime;
	bar1=copyNum;
	k1energy=energy;
	k1zRecon=kzRecon*275.;
	IDNum = analysisManager->GetH1Id("zPos_k1_3BarsOnly",1);
	analysisManager->FillH1(IDNum,k1zRecon);
	k1zRecon=kzRecon;
	angle=15*(copyNum)+7.5;
	angle=angle*3.14/180.;
	k1xRecon=cos(angle);
	k1yRecon=sin(angle);
	k1GammaRayData.push_back(k1energy);
	k1GammaRayData.push_back(k1xRecon);
	k1GammaRayData.push_back(k1yRecon);
	k1GammaRayData.push_back(k1zRecon);
	theta1=TMath::ACos(k1zRecon/TMath::Sqrt((k1xRecon*k1xRecon)+(k1yRecon*k1yRecon)+(k1zRecon*k1zRecon)));
	phi1=TMath::ATan(k1yRecon/k1xRecon);
	if ((k1yRecon>0)&&(k1xRecon>0)){
	  phi1=phi1;
	}
	if ((k1yRecon>0)&&(k1xRecon<0)) {
	  phi1=TMath::Pi()-phi1;
	}
	if ((k1yRecon<0)&&(k1xRecon<0)) {
	  phi1=TMath::Pi()+phi1;
	}
	if ((k1yRecon<0)&&(k1xRecon>0)) {
	  phi1=2*TMath::Pi()-phi1;
	}
	k1XMomentum=k1energy*1000*TMath::Cos(theta1)*TMath::Sin(phi1);
	k1YMomentum=k1energy*1000*TMath::Sin(theta1)*TMath::Sin(phi1);
	k1ZMomentum=k1energy*1000*TMath::Cos(phi1);
      }
      
      if (barCount==2) {
	
	k2hitTime=hitTime;
	bar2=copyNum;
	k2energy=energy;
	k2zRecon=kzRecon*275.;
	IDNum = analysisManager->GetH1Id("zPos_k2_3BarsOnly",1);
	analysisManager->FillH1(IDNum,k2zRecon);
	k2zRecon=kzRecon;
	IDNum = analysisManager->GetH2Id("k1Energy_vs_k2Energy",1);
	analysisManager->FillH2(IDNum,k1energy*1000,k2energy*1000);
	angle=15*(copyNum)+7.5;
	angle=angle*3.14/180.;
	k2xRecon=cos(angle);
	k2yRecon=sin(angle);
	k2GammaRayData.push_back(k2energy);
	k2GammaRayData.push_back(k2xRecon);
	k2GammaRayData.push_back(k2yRecon);
	k2GammaRayData.push_back(k2zRecon);
	theta2=TMath::ACos(k2zRecon/TMath::Sqrt((k2xRecon*k2xRecon)+(k2yRecon*k2yRecon)+(k2zRecon*k2zRecon)));
	phi2=TMath::ATan(k2yRecon/k2xRecon);
	if((k2yRecon>0)&&(k2xRecon>0)){
	  phi2=phi2;
	}
	if((k2yRecon>0)&&(k2xRecon<0)){
	  phi2=TMath::Pi()-phi2;
	}
	if((k2yRecon<0)&&(k2xRecon<0)){
	  phi2=TMath::Pi()+phi2;
	}
	if((k2yRecon<0)&&(k2xRecon>0)){
	  phi2=2*TMath::Pi()-phi2;
	}
	k2XMomentum=k2energy*1000*TMath::Cos(theta2)*TMath::Sin(phi2);
	k2YMomentum=k2energy*1000*TMath::Sin(theta2)*TMath::Sin(phi2);
	k2ZMomentum=k2energy*1000*TMath::Cos(phi2);
      }
      if (barCount==3) {
	
	IDNum = analysisManager->GetH1Id("singles_eSpec_3Bars",1);
	analysisManager->FillH1(IDNum,energy*1000);
	k3hitTime=hitTime;
	IDNum = analysisManager->GetH1Id("hitTime",1);
	analysisManager->FillH1(IDNum,k1hitTime);
	analysisManager->FillH1(IDNum,k2hitTime);
	analysisManager->FillH1(IDNum,k3hitTime);
	bar3=copyNum;
	threeBars.push_back(bar1);
	threeBars.push_back(bar2);
	threeBars.push_back(bar3);
	k3energy=energy;
	k3zRecon=kzRecon*275.;
	IDNum = analysisManager->GetH1Id("zPos_k3_3BarsOnly",1);
	analysisManager->FillH1(IDNum,k3zRecon);
	IDNum = analysisManager->GetH1Id("zPos_allHits_3BarsOnly",1);
	analysisManager->FillH1(IDNum,k3zRecon);
	k3zRecon=kzRecon;
	angle=15*(copyNum)+7.5;
	angle=angle*3.14/180.;
	k3xRecon=cos(angle);
	k3yRecon=sin(angle);
	k3GammaRayData.push_back(k3energy);
	k3GammaRayData.push_back(k3xRecon);
	k3GammaRayData.push_back(k3yRecon);
	k3GammaRayData.push_back(k3zRecon);
	theta3=TMath::ACos(k3zRecon/TMath::Sqrt((k3xRecon*k3xRecon)+(k3yRecon*k3yRecon)+(k3zRecon*k3zRecon)));
	phi3=TMath::ATan(k3yRecon/k3xRecon);
	if ((k3yRecon>0)&&(k3xRecon>0)) {
	  phi3=phi3;
	}
	if ((k3yRecon>0)&&(k3xRecon<0)) {
	  phi3=TMath::Pi()-phi3;
	}
	if ((k3yRecon<0)&&(k3xRecon<0)) {
	  phi3=TMath::Pi()+phi3;
	}
	if ((k3yRecon<0)&&(k3xRecon>0)) {
	  phi3=2*TMath::Pi()-phi3;
	}
	k3XMomentum=k3energy*1000*TMath::Cos(theta3)*TMath::Sin(phi3);
	k3YMomentum=k3energy*1000*TMath::Sin(theta3)*TMath::Sin(phi3);
	k3ZMomentum=k3energy*1000*TMath::Cos(phi3);
	kXMomentumAbs=std::abs(k1XMomentum+k2XMomentum+k3XMomentum);
	kXMomentum=k1XMomentum+k2XMomentum+k3XMomentum;	
	kYMomentumAbs=std::abs(k1YMomentum+k2YMomentum+k3YMomentum);
	kYMomentum=k1YMomentum+k2YMomentum+k3YMomentum;
	kZMomentumAbs=std::abs(k1ZMomentum+k2ZMomentum+k3ZMomentum);
	kZMomentum=k1ZMomentum+k2ZMomentum+k3ZMomentum;
	IDNum = analysisManager->GetH1Id("X_mom",1);
	analysisManager->FillH1(IDNum,kXMomentum/1000.);
	IDNum = analysisManager->GetH1Id("Y_mom",1);
	analysisManager->FillH1(IDNum,kYMomentum/1000.);
	IDNum = analysisManager->GetH1Id("Z_mom",1);
	analysisManager->FillH1(IDNum,kZMomentum/1000.);
	IDNum = analysisManager->GetH1Id("mag_vector_sum",1);
	analysisManager->FillH1(IDNum,(kXMomentumAbs+kYMomentum+kZMomentum)/1000.);
      }
      
      if (barCount>2) {
	reconstruction.DoAnalysisCuts(totalEnergy,hitTime,(kXMomentumAbs+kYMomentum+kZMomentum)/1000.,k1energy,k2energy,k3energy,k1GammaRayData,k2GammaRayData,k3GammaRayData,threeBars,hitBarsSet.size());
	//	RyanReconstruction(k1GammaRayData,k2GammaRayData,k3GammaRayData,threeBars,hitBarsSet.size());
	//YamazakiReconstruction(k1GammaRayData,k2GammaRayData,k3GammaRayData,threeBars,hitBarsSet.size());
      }
      barCount++;
    }
    IDNum = analysisManager->GetH1Id("numBarsHit",1);
    analysisManager->FillH1(IDNum,gammaRayVector.size());    
  }

}
TH1F* EventAction::GetoPSTimingHistogramPtr(){
  return timingHistoPS;
}

TH1F* EventAction::GetpPSTimingHistogramPtr(){
  return timingHistpPS;
}

TH1F* EventAction::GetoPS_2GammaTimingHistogramPtr(){
  return timingHistoPS_2Gamma;
}

TH1F* EventAction::GetpPS_2GammaTimingHistogramPtr(){
  return timingHistpPS_2Gamma;
}

TH1F* EventAction::GetoPS_3GammaTimingHistogramPtr(){
  return timingHistoPS_3Gamma;
}

TH1F* EventAction::GetpPS_3GammaTimingHistogramPtr(){
  return timingHistpPS_3Gamma;
}

