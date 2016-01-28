#include <cstdlib>
#include <ctime>
#include <fstream>
//#include "ThreeGammaCoordinates.hh"
#include "G4Point3D.hh"
#include "Visualize.hh"
//#include "ThetaDistribution.hh"
#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorActionThreeGamma.hh"
#include "PrimaryGeneratorActionOPS.hh"
#include "PrimaryGeneratorMessenger.hh"
#include "G4Gamma.hh"
#include "G4Polyline.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4Geantino.hh"
#include "G4Electron.hh"
#include "DetectorConstruction.hh"
#include "G4RotationMatrix.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4GeneralParticleSource.hh"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction *DC) 
 :Detector(DC)
{

  twoGamma=0;
  threeGamma=0;
  // default particle kinematic

  G4int n_particle = 1;
  
  fParticleGun = new G4ParticleGun();   

  particleGunGPS = new G4GeneralParticleSource();
 
  std::srand(std::time(0)); 

  //create messenger for this class
  selectedAction=1;
  gunMessenger = new PrimaryGeneratorMessenger(this);

  //Two Gamma
  particleGunTwoGamma = new G4ParticleGun();
  //  actionThreeGamma = new PrimaryGeneratorActionOPS(particleGunTwoGamma,DC);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particleGunGPS;  
  delete gunMessenger;      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  switch(selectedAction)
    {
    case 0:
      {	
	break;
      }
    case 1:
      //Generate an ordinary GPS event
      particleGunGPS->GeneratePrimaryVertex(anEvent);
      break;
    case 2:
      //Generate Back-to-Back Gamma Ray Event
      //      actionThreeGamma->GeneratePrimaries(anEvent);
      break;

    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
