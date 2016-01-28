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
#include "G4VUserEventInformation.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include <set>
#include <math.h>

#ifndef UserEventInformation_h
#define UserEventInformation_h 1

class UserEventInformation : public G4VUserEventInformation
{
public:
  UserEventInformation();
  ~UserEventInformation();
  
  inline void Print()const{};
 
  void SetMuonVeto(){muonveto=1;}
  G4bool IsMuonVeto(){return muonveto;}

  void SetCopyNo(G4int copyNumber);
  G4int GetCopyNo(){return copyNo;}
  void AddTotA(G4double a1, G4double a2, G4int copyNumber);
  void AddToSetOfHitBars(G4int);
  std::set<G4int> GetSetOfHitBars();
  G4double GetZRecon();
  void AddTotZ(G4double, G4int);//{totZ+=zPos;}
  G4double GetTotZ(G4int);//{return totZ;}
  void AddCounter(){counter++;}
  G4int GetCounter(){return counter;}
  void AddEdep(G4double,G4int);
  G4double GetEnergyDeposition(G4int);
  G4double CalculateZRecon(G4int);

  void AddHitTime(G4double,G4int);
  G4double GetHitTime(G4int);
  G4int GetNumStepsInBar(G4int);
  // Define array for the bars
  // Should this be const?
  
private:

  std::set<G4int> hitBarsSet;
  //  extern G4double totalTotalZ[24];
  G4double* attenuationCoefficients;
  G4double* totalAmplitude1Array;
  G4double* totalAmplitude2Array;
  G4double* totalEDeposition;
  G4double* barHitTime;
  G4int* stepCounter;
  G4double* zReconArray;
  G4double* zPosition;
  //  extern G4double eDeposition[24];
  G4double zRecon6;
  G4bool muonveto;
  G4int copyNo;
  G4double totA1;
  G4double totA2;
  G4double totZ;
  G4int counter;
  G4double edepTot;
};

#endif





