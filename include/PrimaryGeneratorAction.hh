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
//
// $Id: PrimaryGeneratorAction.hh,v 1.3 2010/07/16 07:37:48 maire Exp $
// GEANT4 tag $Name: geant4-09-04 $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include <string>
#include "G4ThreeVector.hh"
#include <fstream>
//#include "ThetaDistribution.hh"
class PrimaryGeneratorMessenger;
class PrimaryGeneratorActionTwoGamma;
class PrimaryGeneratorActionThreeGamma;
class PrimaryGeneratorActionOPS;
class PrimaryGeneratorActionGammaSource;
class PrimaryGeneratorActionGammaDist;
class PrimaryGeneratorAction208Tl;
class PrimaryGeneratorAction228Ac;
class PrimaryGeneratorAction214Bi;
class PrimaryGeneratorActionFlat;
class PrimaryGeneratorAction40K;
class PrimaryGeneratorActionMuon_plus;
class PrimaryGeneratorActionMuon_minus;
class PrimaryGeneratorActionCRY;
class PrimaryGeneratorActionMeteorite;
class PrimaryGeneratorAction14N;
class PrimaryGeneratorActionGasJet;
class PrimaryGeneratorActionGen;

class G4ParticleGun;
class G4GeneralParticleSource;
class G4Event;
class DetectorConstruction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction();    
  PrimaryGeneratorAction(DetectorConstruction*);    
  ~PrimaryGeneratorAction();
  
public:
  void GeneratePrimaries(G4Event*);
  
public:
  //G4ParticleGun* GetParticleGun() { return particleGun; };
  
  void SelectAction(G4int i) { selectedAction = i; };    
  G4int GetSelectedAction()  { return selectedAction; }; 
  //void SelectEnergy(G4double ene) { userEnergy = ene; };    
  //G4double GetUserEnergy()  { return userEnergy; };    
  PrimaryGeneratorActionGammaSource*  GetActionGammaSource() { return actionGammaSource; };
  PrimaryGeneratorActionGammaDist*  GetActionGammaDist() { return actionGammaDist; };
  PrimaryGeneratorAction208Tl*  GetAction208Tl() { return action208Tl; };
  PrimaryGeneratorAction228Ac*  GetAction228Ac() { return action228Ac; };     
  PrimaryGeneratorAction214Bi*  GetAction214Bi() { return action214Bi; };     
  PrimaryGeneratorActionFlat*  GetActionFlat() { return actionFlat; };     
  PrimaryGeneratorAction40K*  GetAction40K() { return action40K; };     
  PrimaryGeneratorActionMuon_plus*  GetActionMuon_plus() { return actionMuon_plus; };     
  PrimaryGeneratorActionMuon_minus*  GetActionMuon_minus() { return actionMuon_minus; };     
  PrimaryGeneratorActionCRY*  GetActionCRY() { return actionCRY; };     
  PrimaryGeneratorActionMeteorite*  GetActionMeteorite() { return actionMeteorite; }; 
  PrimaryGeneratorAction14N*  GetAction14N() { return action14N; };         
  PrimaryGeneratorActionGasJet*  GetActionGasJet() { return actionGasJet; };     
  PrimaryGeneratorActionGen*  GetActionGen() { return actionGen; };     
  
  G4double GetParticleEnergy(){return energy;};
  // Get number of particles produced for ascii file stuff
  G4int GetNumberOfParticles(){return nParticles;};       
    
private:

  G4ThreeVector            fCurrentPosition;
  G4ParticleGun*           fParticleGun;
  G4ParticleGun*           particleGunTwoGamma;
  G4ParticleGun*           particleGunThreeGamma;
  G4ParticleGun*           particleGunGammaSource;
  G4ParticleGun*           particleGunGammaDist;
  G4ParticleGun*           particleGun208Tl;
  G4ParticleGun*           particleGun228Ac;
  G4ParticleGun*           particleGun214Bi;
  G4ParticleGun*           particleGunFlat;
  G4ParticleGun*           particleGun40K;
  G4ParticleGun*           particleGunMuon_plus;
  G4ParticleGun*           particleGunMuon_minus;
  G4ParticleGun*           particleGunCRY;
  G4ParticleGun*           particleGunMeteorite;
  G4ParticleGun*           particleGun14N;
  G4ParticleGun*           particleGunGasJet;
  G4ParticleGun*           particleGunGen;
  G4GeneralParticleSource *particleGunGPS;
  PrimaryGeneratorActionTwoGamma* actionTwoGamma;
  PrimaryGeneratorActionOPS* actionThreeGamma;
  PrimaryGeneratorActionThreeGamma* actionThreeGammaP1;
  PrimaryGeneratorActionThreeGamma* actionThreeGammaM1;
  PrimaryGeneratorActionGammaSource* actionGammaSource;
  PrimaryGeneratorActionGammaDist* actionGammaDist;
  PrimaryGeneratorAction208Tl* action208Tl;
  PrimaryGeneratorAction228Ac* action228Ac;
  PrimaryGeneratorAction214Bi* action214Bi;
  PrimaryGeneratorActionFlat* actionFlat;
  PrimaryGeneratorAction40K* action40K;
  PrimaryGeneratorActionMuon_plus* actionMuon_plus;
  PrimaryGeneratorActionMuon_minus* actionMuon_minus;
  PrimaryGeneratorActionCRY* actionCRY;
  PrimaryGeneratorActionMeteorite* actionMeteorite;
  PrimaryGeneratorAction14N* action14N;
  PrimaryGeneratorActionGasJet* actionGasJet;
  PrimaryGeneratorActionGen* actionGen;
  //  ThetaDistribution* thetadistr;
  DetectorConstruction *Detector;
  G4int                    selectedAction;
  //G4double userEnergy;
  std::string caseString;
  PrimaryGeneratorMessenger* gunMessenger;   
  G4double energy;
  G4int nParticles;

  //For Three Gamma Decay
  G4ThreeVector rotatedCoordinates;
  double randomNum;

  std::ofstream outputFile;
  int twoGamma;
  int threeGamma;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
