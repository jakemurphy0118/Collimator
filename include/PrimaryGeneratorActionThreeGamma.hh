//--------------------------------------------------------------------
// PrimaryGeneratorAction.hh
//
// Description: Sorts out generation, events etc. inlcude file
// Changes: 7/18/05 None yet
//--------------------------------------------------------------------

#ifndef PrimaryGeneratorActionThreeGamma_h
#define PrimaryGeneratorActionThreeGamma_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "PrimaryGeneratorActionThreeGamma.hh"
//#include "ThetaDistribution.hh"
#include "globals.hh"
#include <vector>
#include "TH1F.h"
#include <fstream>

class G4ParticleGun;
class G4Event;
class DetectorConstruction;
class PrimaryGeneratorMessengerGammaSource;

using namespace CLHEP;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorActionThreeGamma : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorActionThreeGamma(G4ParticleGun*,DetectorConstruction*,TH1F*);    
  ~PrimaryGeneratorActionThreeGamma();
  
public:
  void GeneratePrimaries(G4Event*);
  void SelectEnergy(G4double ene) { userEnergy = ene; };    
  G4double GetUserEnergy()  { return userEnergy; };  

  // Use this for ascii file naming
  G4double GetParticleEnergy(){return energy;};
  // Get number of particles produced for ascii file stuff
  G4int GetNumberOfParticles(){return nParticles;};
  void setThetaHist(TH1F*);
  
public:
  G4ParticleGun* GetParticleGun() { return particleGun; };
  
private:

  TH1F* thetaHist1;
  G4ParticleGun*  particleGun;
  DetectorConstruction *Detector;
  PrimaryGeneratorMessengerGammaSource* gunMessenger;   
  G4double userEnergy;
        
  G4double energy;
  G4int nParticles;
  Double_t theta;
  G4ThreeVector rotatedCoordinates_1;
  G4ThreeVector rotatedCoordinates_2;
  G4ThreeVector fCurrentPosition;

  std::ofstream planeFile;
  std::ofstream normFile;
  std::ofstream angleFile;
  std::ofstream fooey;
  std::ofstream allFile;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
