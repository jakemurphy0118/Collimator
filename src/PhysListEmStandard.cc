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
// $Id: PhysListEmStandard.cc,v 1.15 2009/11/20 20:27:21 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-03 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "PhysListEmStandard.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4UserSpecialCuts.hh"
#include "G4PSGeneration.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eMultipleScattering.hh"
#include "G4hMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4MuMultipleScattering.hh"
#include "G4MuonMinusCapture.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"

#include "G4EmProcessOptions.hh"
#include "G4MscStepLimitType.hh"

#include "G4RadioactiveDecay.hh"

using namespace CLHEP;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmStandard::PhysListEmStandard(const G4String& name,G4double BFieldInput)
   :  G4VPhysicsConstructor(name)
{
  BField=BFieldInput;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmStandard::~PhysListEmStandard()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysListEmStandard::ConstructProcess()
{
  // Add standard EM Processes
  G4cout << " <PhysListEmStandard::ConstructProcess> Built Physics List" << G4endl;
  //G4MuMultipleScattering *mums = new G4MuMultipleScattering();
  //G4MuIonisation *muion = new G4MuIonisation();
  //G4MuBremsstrahlung *mubrem = new G4MuBremsstrahlung();
  //G4MuPairProduction *mupair = new G4MuPairProduction();

  aParticleIterator->reset();
  while( (*aParticleIterator)() ){
    G4ParticleDefinition* particle = aParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {
     
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect);
      pmanager->AddDiscreteProcess(new G4ComptonScattering);
      pmanager->AddDiscreteProcess(new G4GammaConversion);
      
    } else if (particleName == "e-") {
  
      pmanager->AddProcess(new G4eMultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4eIonisation(),       -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung,     -1, 3, 3);

      //pmanager->AddDiscreteProcess(new G4UserSpecialCuts());

	    
    } else if (particleName == "e+") {

      pmanager->AddProcess(new G4eMultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4eIonisation(),       -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung,     -1, 3, 3);
      //      pmanager->AddProcess(new G4eplusAnnihilation,    0,-1, 4);
      G4cout << "PS STUFF IS WORKING" << G4endl;
      pmanager->AddProcess(new G4PSGeneration("PSDecay",BField), 1,-1, -1);

      
    } else if (particleName == "mu+" || 
               particleName == "mu-"    ) {

      //G4MuBremsstrahlung *mubrem = new G4MuBremsstrahlung;
    
      pmanager->AddProcess(new G4MuMultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4MuIonisation,        -1, 2, 2);//knock-on/delta electrons
      pmanager->AddProcess(new G4MuBremsstrahlung,    -1, 3, 3);
      pmanager->AddProcess(new G4MuPairProduction,    -1, 4, 4);    

      if( particleName == "mu-" )
	pmanager->AddProcess(new G4MuonMinusCapture(), 0,-1,-1);
   
      //pmanager->AddProcess(mums, -1, 1, 1);
      //pmanager->AddProcess(muion,        -1, 2, 2);//knock-on/delta electron
      //pmanager->AddProcess(mubrem,    -1, 3, 3);
      //pmanager->AddProcess(mupair,    -1, 4, 4);


      //pmanager->RemoveProcess(mums);
      //pmanager->RemoveProcess(muion);
      //pmanager->RemoveProcess(mubrem);
      //pmanager->RemoveProcess(mupair);
    } 
    //The following particles are needed for radioactive decay (eg. 60Co GPS), must be defined. 
    else if (particleName == "proton" || 
               particleName == "pi-" ||  
	       particleName == "pi+") {

      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4hIonisation,         -1, 2, 2);
      pmanager->AddProcess(new G4hBremsstrahlung,     -1, 3, 3);
      pmanager->AddProcess(new G4hPairProduction,     -1, 4, 4);       

    } else if (particleName == "alpha" ||
	       particleName == "He3" || 
	       particleName == "GenericIon") {

      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4ionIonisation,       -1, 2, 2);
     
    } else if ((!particle->IsShortLived()) &&
	       (particle->GetPDGCharge() != 0.0) && 
	       (particle->GetParticleName() != "chargedgeantino")) {

      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4hIonisation,         -1, 2, 2);
    }

  }
  
  // Em options
  //
  // Main options and setting parameters are shown here.
  // Several of them have default values.
  //  
  G4EmProcessOptions emOptions;
  
  //physics tables
  //
  emOptions.SetMinEnergy(100*eV);	//default    
  emOptions.SetMaxEnergy(100*TeV);	//default  
  emOptions.SetDEDXBinning(12*20);	//default=12*7  
  emOptions.SetLambdaBinning(12*20);	//default=12*7
  emOptions.SetSplineFlag(true);	//default  
    
  //coulomb scattering
  //
  emOptions.SetMscStepLimitation(fUseDistanceToBoundary);  //default=fUseSafety
  emOptions.SetMscRangeFactor(0.04);	//default
  emOptions.SetMscGeomFactor (2.5);	//default       
  emOptions.SetSkin(3.);		//default
          
  //energy loss
  //
  emOptions.SetStepFunction(0.2, 100*um);	//default=(0.2, 1*mm)   
  emOptions.SetLinearLossLimit(1.e-2);		//default
   
  //ionization
  //
  emOptions.SetSubCutoff(true);		//default=false  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

