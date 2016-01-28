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
// $Id: PrimaryGeneratorMessenger.cc,v 1.1 2010/07/16 07:37:48 maire Exp $
// GEANT4 tag $Name: geant4-09-04 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorMessenger.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(
                                             PrimaryGeneratorAction* Gun)
:Action(Gun)
{
  Dir = new G4UIdirectory("/gunAction/");
  Dir->SetGuidance("Modular sources for LENA");
    
  selectActionCmd = new G4UIcmdWithAnInteger("/gunAction/selectGunAction",this);
  selectActionCmd->SetGuidance("Select primary generator action");
  selectActionCmd->SetGuidance(" id = 0 : General Particle Source");
  selectActionCmd->SetGuidance(" id = 1 : Gamma point source, emitted isotropically");
  selectActionCmd->SetGuidance(" id = 2 : Gamma source -- distributed adn emitted isotropically");
  selectActionCmd->SetGuidance(" id = 3 : 208Tl room background");
  selectActionCmd->SetGuidance(" id = 4 : 228Ac room background");
  selectActionCmd->SetGuidance(" id = 5 : 214Bi room background");
  selectActionCmd->SetGuidance(" id = 6 : flat room background (0<E<50 MeV)");
  selectActionCmd->SetGuidance(" id = 7 : 40K room background");  
  selectActionCmd->SetGuidance(" id = 8 : All room background");  
  selectActionCmd->SetGuidance(" id = 9 : Muons +");  
  selectActionCmd->SetGuidance(" id = 10 : Muons -");  
  selectActionCmd->SetGuidance(" id = 11 : Cosmic ray library (CRY)");  
  selectActionCmd->SetGuidance(" id = 12 : Meteorite");  
  selectActionCmd->SetGuidance(" id = 13 : 14Na Source");  
  selectActionCmd->SetGuidance(" id = 14 : Gas Jet");  
  selectActionCmd->SetGuidance(" id = 15 : General Reaction");  
  selectActionCmd->SetParameterName("id",false);
  selectActionCmd->SetRange("id>=0 && id<16");
  selectActionCmd->SetDefaultValue(0);

  //selectEnergyCmd = new G4UIcmdWithADoubleAndUnit("/gunAction/selectEnergy",this);
  //selectEnergyCmd->SetParameterName("id",false);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete selectActionCmd;
  //delete selectEnergyCmd;
  delete Dir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command,
                                               G4String newValue)
{ 
  if (command == selectActionCmd)
    Action->SelectAction(selectActionCmd->GetNewIntValue(newValue)); 

  //if (command == selectEnergyCmd)
  //Action->SelectEnergy(selectEnergyCmd->GetNewDoubleValue(newValue));      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

