//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id$
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4eplusAnnihilation
//
// Author:        Vladimir Ivanchenko on base of Michel Maire code
//
// Creation date: 02.08.2004
//
// Modifications:
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivanchenko)
// 08-04-05 Major optimisation of internal interfaces (V.Ivanchenko)
// 03-05-05 suppress Integral option (mma)
// 04-05-05 Make class to be default (V.Ivanchenko)
// 25-01-06 remove cut dependance in AtRestDoIt (mma)
// 09-08-06 add SetModel(G4VEmModel*) (mma)
// 12-09-06, move SetModel(G4VEmModel*) in G4VEmProcess (mma)
// 30-05-12 propagate parent weight to secondaries (D. Sawkey)
//

//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include <time.h>
#include "PsGammaDecayGenerator.hh"
#include "G4PSGeneration.hh"
#include "G4PhysicalConstants.hh"
//#include "G4Positronium.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4eeToTwoGammaModel.hh"
#include "G4ThreeVector.hh"
#include <iostream>
#include <fstream>
#include <TROOT.h>
#include <TRandom3.h>
#include <math.h>
#include <sys/time.h>
#include <TVector3.h>
#include <TRandom.h>
#include <TRotation.h>
#include <TH3F.h>
#include <TF3.h>
//#include "TimeDistribution.h"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#define PI 3.14159
using namespace CLHEP;
using namespace std;

G4PSGeneration::G4PSGeneration(const G4String& name,double BField)
  : G4VEmProcess(name), isInitialised(false)
{
  //Set the magnetic field
  H=BField;

  //Calculate the branching ratio
  double mu = 9.27e-21;
  double g = 2;
  double deltaE=1.28e-15;
  double y = (mu*g*H)/(deltaE/2);
  double gamma_t=(1/(142e-9));
  double gamma_s=(1/(.125e-9));
  BR_pPS=1;
  BR_oPS=1;

  //Set the branching ratios for the perturbed triplet and perturbed singlet states
  BR_oPS=(((2+y*y)-2*(sqrt(1+y*y)))/(y*y))*(gamma_s/gamma_t);
  BR_pPS=(((2+y*y)+2*(sqrt(1+y*y)))/(y*y))*(gamma_s/gamma_t);

  outfile.open("~/Caliope_Work/CALIOPE-BUILD/file.dat");
  theGamma = G4Gamma::Gamma();
  //  thePositronium = G4Positronium::Positronium();
  SetIntegral(true);
  SetBuildTableFlag(false);
  SetStartFromNullFlag(false);
  SetSecondaryParticle(theGamma);
  SetProcessSubType(fAnnihilation);
  enableAtRestDoIt = true;
  //  thetadistr = new ThetaDistribution();
  //spinHist=thetadistr->buildHist("spin0");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PSGeneration::~G4PSGeneration()
{
  //  h1->Write();
  //  hfile->Write();
  //  delete h1;
  //  delete spinHist;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4PSGeneration::IsApplicable(const G4ParticleDefinition& p)
{
  return (&p == G4Positron::Positron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PSGeneration::AtRestGetPhysicalInteractionLength(
                              const G4Track&, G4ForceCondition* condition)
{

  *condition = NotForced;
  return 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PSGeneration::InitialiseProcess(const G4ParticleDefinition*)
{
  const G4long* table_entry;
  G4long theSeed;

  theSeed = CLHEP::HepRandom::getTheSeed();
  const long* seeds;
  seeds = CLHEP::HepRandom::getTheSeeds();
  gen = new PsGammaDecayGenerator();
  gen->Initialize(); // This is important!             
  if(!isInitialised) {
    isInitialised = true;
    if(!EmModel(1)) { SetEmModel(new G4eeToTwoGammaModel(),1); }
    EmModel(1)->SetLowEnergyLimit(MinKinEnergy());
    EmModel(1)->SetHighEnergyLimit(MaxKinEnergy());
    AddEmModel(1, EmModel(1));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PSGeneration::PrintInfo()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4PSGeneration::AtRestDoIt(const G4Track& aTrack,
                                                     const G4Step& stepData){
	G4VParticleChange* fParticleChange;

      fParticleChange=Do2GammaDecay(aTrack, stepData, 0, 0, "none");
    
              
  return fParticleChange;
  
}

G4VParticleChange* G4PSGeneration::Do2GammaDecay(const G4Track& aTrack,const G4Step& stepData, G4int m, G4bool BField, G4String perturbation){

  //p-PS
  if(BField!=0) {
    // Perturbed p-Ps: 2 gamma decay
    if(perturbation=="psi_00")
      SetProcessSubType(fppsDecay_2Gamma);
    // Perturbed o-Ps: 2 gamma decay
    if(perturbation=="psi_10")
      SetProcessSubType(fopsDecay_2Gamma);
  }
  else{
    SetProcessSubType(fppsDecay);
  }

  
  //Not sure why we have to initialize for post step
  fParticleChange.InitializeForPostStep(aTrack);
  //Not sure what this does either
  fParticleChange.ProposeWeight(aTrack.GetWeight());
  
  G4double xdir=2.*G4UniformRand()-1.;
  G4double ydir=2.*G4UniformRand()-1.;
  G4double zdir=2.*G4UniformRand()-1.;
  G4double norm=sqrt(xdir*xdir+ydir*ydir+zdir*zdir);
  xdir=xdir/norm;
  ydir=ydir/norm;
  zdir=zdir/norm;
  
  
  G4double cosTeta = 2.*G4UniformRand()-1.;
  G4double sinTeta = sqrt((1.-cosTeta)*(1.0 + cosTeta));
  G4double phi     = twopi * G4UniformRand();
  G4ThreeVector dir(sinTeta*cos(phi), sinTeta*sin(phi), cosTeta);      
  
  //      G4ThreeVector kappa1(xdir,ydir,zdir);
  
  fParticleChange.SetNumberOfSecondaries(2);
  G4DynamicParticle* dp = 
    new G4DynamicParticle(theGamma, dir, 511*keV);
  
  fParticleChange.AddSecondary(dp);
  dp = new G4DynamicParticle(theGamma,-dir, 511*keV);
  
  fParticleChange.AddSecondary(dp);
  //      dp = new G4DynamicParticle(theGamma,kappa3, 10*keV);
  
  //      fParticleChange.AddSecondary(dp);
  
  //kill the incident positron
  fParticleChange.ProposeTrackStatus(fStopAndKill);
  
  return &fParticleChange;
}

G4VParticleChange* G4PSGeneration::Do3GammaDecay(const G4Track& aTrack,const G4Step& stepData,G4int m,G4bool BField,G4String perturbation){
  
  //o-PS
  if(BField!=0) {
    if(m==0){
      // Perturbed p-Ps: 3 gamma decay
      if(perturbation=="psi_00")
	SetProcessSubType(fppsDecay_3Gamma);
      // Perturbed o-Ps: 3 gamma decay
      if(perturbation=="psi_10")
	SetProcessSubType(fopsDecay_3Gamma);
    }
    else{
      SetProcessSubType(fopsDecay);
    }
  }
  else{
    SetProcessSubType(fopsDecay);
  }
  
  //Not sure why we have to initialize for post step
  fParticleChange.InitializeForPostStep(aTrack);
  //Not sure what this does either
  fParticleChange.ProposeWeight(aTrack.GetWeight());
  
  TVector3 k1, k2, k3, normal, sprojection;
  Double_t theta, psi, phi, asymmetryParameter;
  //==========
  //Get the momenta
  
  // The next line generates the vertices. Eveything else is diagnostic.                                         
  //      gen->GenerateoPs_CP(&k1, &k2, &k3); // CP-violating case                                             
  
  //      clockid_t clockID;
  //struct timespec res;
  //struct timespec tp;
  //int clockRes=clock_getres(CLOCK_MONOTONIC,&res);
  //      if(!clockRes){
  //      G4cout << "Clock Resolution: " << res.tv_nsec << G4endl;
  //clock_gettime(CLOCK_MONOTONIC,&tp);
  //unsigned before=tp.tv_nsec;
  //      G4cout << "Before: " << before << G4endl;
  
  //      if(BField==0)
  //      gen->GenerateoPs(&k1, &k2, &k3, 0); // Bernreuther m=0
  //      if(BField==1)
  //	gen->GenerateoPs(&k1, &k2, &k3, 1); // Bernreuther m=1
  //      if(BField==-1)
  //	gen->GenerateoPs(&k1, &k2, &k3, 1); // Bernreuther m=-1
  
  
  //clock_gettime(CLOCK_MONOTONIC,&tp);
  //unsigned after=tp.tv_nsec;
  //      unsigned after=res.tv_nsec;
  //      G4cout << "After: " << after << G4endl;
  //float duration = ((float)(before-after))/((float)(res.tv_nsec));
  ////      G4cout << "Duration: " << duration << G4endl;
  //}
  gen->GenerateoPs(&k1, &k2, &k3, 0); // Bernreuther m=+-1    
  
  //=========================
  //======================
  // Test angle recon
  // By setting psi to 90 degrees
  //======================
  
  double k1_Energy=k1.Mag()*511*keV;
  double k2_Energy=k2.Mag()*511*keV;
  double k3_Energy=k3.Mag()*511*keV;
  
  G4ThreeVector kappa1(k1.X(),k1.Y(),k1.Z());
  G4ThreeVector kappa2(k2.X(),k2.Y(),k2.Z());
  //G4ThreeVector kappa1(1,0,0);
  //G4ThreeVector kappa2(-1,0,0);
  G4ThreeVector kappa3(k3.X(),k3.Y(),k3.Z());
  kappa1=kappa1.unit();
  kappa2=kappa2.unit();
  kappa3=kappa3.unit();
  
  if(isnan(k1.Mag())==0){
    fParticleChange.SetNumberOfSecondaries(3);
    G4DynamicParticle* dp = 
      new G4DynamicParticle(theGamma, kappa1, k1_Energy);

    fParticleChange.AddSecondary(dp);
    G4DynamicParticle* dp1 = 
      new G4DynamicParticle(theGamma, kappa2, k2_Energy);
    
    fParticleChange.AddSecondary(dp1);
    G4DynamicParticle* dp2 = 
      new G4DynamicParticle(theGamma, kappa3, k3_Energy);
    
    fParticleChange.AddSecondary(dp2);    
  }
  
  fParticleChange.ProposeTrackStatus(fStopAndKill);
  
  return &fParticleChange;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
