//--------------------------
// PrimaryGeneratorActionTwoGamma.cc
// Description:
// Yak yak yak
// 
//--------------------------

#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorActionOPS.hh"
#include "PsGammaDecayGenerator.hh"
#include "DetectorConstruction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4Gamma.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
//#include "ThreeGammaCoordinates.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include <fstream>
#include "G4Gamma.hh"

using namespace CLHEP;

PrimaryGeneratorActionOPS::PrimaryGeneratorActionOPS(G4ParticleGun *gun, DetectorConstruction *DC) :particleGun(gun),Detector(DC)
{
  gen = new PsGammaDecayGenerator();
  gen->Initialize(); // This is important!
}
PrimaryGeneratorActionOPS::~PrimaryGeneratorActionOPS()
{
  //delete particleGun;
  delete gunMessenger;
}
void PrimaryGeneratorActionOPS::GeneratePrimaries(G4Event* anEvent)
{
  // Generate back-to-back Gamma Rays

  TVector3 k1, k2, k3, normal, sprojection;
  Double_t theta, psi, phi, asymmetryParameter;


  delete gRandom;
  gRandom = new TRandom3(0); // Use UUID as seed. 
  double rand=gRandom->Rndm();

  if(rand<0.33){
    //    gen->GenerateoPs(&k1, &k2, &k3, 0); // Bernreuther m=0                                                                                                           
  }
  else{
    //    gen->GenerateoPs(&k1, &k2, &k3, 1); // Bernreuther m=+-1   
  }
  //  gen->GenerateoPs_CP(&k1, &k2, &k3); // CP-violating case                                                                                               

  G4ThreeVector kappa1(k1.X(),k1.Y(),k1.Z());
  G4ThreeVector kappa2(k2.X(),k2.Y(),k2.Z());
  G4ThreeVector kappa3(k3.X(),k3.Y(),k3.Z());
  kappa1=kappa1.unit();
  kappa2=kappa2.unit();
  kappa3=kappa3.unit();

  double k1_Energy=k1.Mag()*1022*keV;
  double k2_Energy=k2.Mag()*1022*keV;
  double k3_Energy=k3.Mag()*1022*keV;

  /*  particleGun->SetParticleDefinition(G4Gamma::Gamma());
  particleGun->SetParticleEnergy(k1.Mag());
  particleGun->SetParticleMomentumDirection(kappa1);
  particleGun->GeneratePrimaryVertex(anEvent);

  particleGun->SetParticleDefinition(G4Gamma::Gamma());
  particleGun->SetParticleEnergy(k2.Mag());
  particleGun->SetParticleMomentumDirection(kappa2);
  particleGun->GeneratePrimaryVertex(anEvent);

  particleGun->SetParticleDefinition(G4Gamma::Gamma());
  particleGun->SetParticleEnergy(k3.Mag());
  particleGun->SetParticleMomentumDirection(kappa3);
  particleGun->GeneratePrimaryVertex(anEvent);*/

}
