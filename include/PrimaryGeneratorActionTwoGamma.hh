//--------------------------------------------------------------------
// PrimaryGeneratorAction.hh
//
// Description: Sorts out generation, events etc. inlcude file
// Changes: 7/18/05 None yet
//--------------------------------------------------------------------

#ifndef PrimaryGeneratorActionTwoGamma_h
#define PrimaryGeneratorActionTwoGamma_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "PrimaryGeneratorActionTwoGamma.hh"
#include "globals.hh"
#include <vector>

class G4ParticleGun;
class G4Event;
class DetectorConstruction;
class PrimaryGeneratorMessengerGammaSource;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorActionTwoGamma : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorActionTwoGamma(G4ParticleGun*,DetectorConstruction*);    
  ~PrimaryGeneratorActionTwoGamma();
  
public:
  void GeneratePrimaries(G4Event*);
  void SelectEnergy(G4double ene) { userEnergy = ene; };    
  G4double GetUserEnergy()  { return userEnergy; };  

  // Use this for ascii file naming
  G4double GetParticleEnergy(){return energy;};
  // Get number of particles produced for ascii file stuff
  G4int GetNumberOfParticles(){return nParticles;};
  
public:
  G4ParticleGun* GetParticleGun() { return particleGun; };
  
private:
  G4ParticleGun*  particleGun;
  DetectorConstruction *Detector;
  PrimaryGeneratorMessengerGammaSource* gunMessenger;   
  G4double userEnergy;
        
  G4double energy;
  G4int nParticles;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
