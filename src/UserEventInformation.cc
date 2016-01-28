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
#include "UserEventInformation.hh"
#include "AnalysisManager.hh"
#include <set>
//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
UserEventInformation::UserEventInformation()
  :muonveto(0)//mupairprod(0),muion(0),mumsc(0),mubrem(0),mudecay(0)
{
  // Initialize attenuation coefficients
  attenuationCoefficients = new G4double[24];
  attenuationCoefficients[0] = 0.00447;
  attenuationCoefficients[1] = 0.00454;
  attenuationCoefficients[2] = 0.00413;
  attenuationCoefficients[3] = 0.00454;
  attenuationCoefficients[4] = 0.00441;
  attenuationCoefficients[5] = 0.00465;
  attenuationCoefficients[6] = 0.00411;
  attenuationCoefficients[7] = 0.00425;
  attenuationCoefficients[8] = 0.00350;
  attenuationCoefficients[9] = 0.00359;
  attenuationCoefficients[10] = 0.00404;
  attenuationCoefficients[11] = 0.00381;
  attenuationCoefficients[12] = 0.00329;
  attenuationCoefficients[13] = 0.00368;
  attenuationCoefficients[14] = 0.00355;
  attenuationCoefficients[15] = 0.00372;
  attenuationCoefficients[16] = 0.00368;
  attenuationCoefficients[17] = 0.00367;
  attenuationCoefficients[18] = 0.00363;
  attenuationCoefficients[19] = 0.00425;
  attenuationCoefficients[20] = 0.00390;
  attenuationCoefficients[21] = 0.00398;
  attenuationCoefficients[22] = 0.00380;
  attenuationCoefficients[23] = 0.00405;
  
  
  //  totalEDeposition
  totalAmplitude1Array = new G4double[24];
  totalAmplitude2Array = new G4double[24];
  totalEDeposition = new G4double[24];
  barHitTime = new G4double[24];
  stepCounter = new G4int[24];
  zPosition = new G4double[24];
  for(int i=0;i<24;i++){
    totalAmplitude1Array[i]=0;
    totalAmplitude2Array[i]=0;
    totalEDeposition[i]=0;
    barHitTime[i]=0;
    stepCounter[i]=0;
    zPosition[i]=0;
  }

  zReconArray = new G4double[24];
  muonveto=0;
  copyNo=0;
  counter=0;
  totZ=0;
  totA1=0;
  totA2=0;
  edepTot=0;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
UserEventInformation::~UserEventInformation()
{
  delete totalAmplitude1Array;
  delete totalAmplitude2Array;
  delete totalEDeposition;
  delete barHitTime;
  delete stepCounter;
  delete zPosition;
  delete zReconArray;
}

void UserEventInformation::SetCopyNo(G4int copyNumber)
{
    copyNo = copyNumber;
/*    totA1=0;
    totA2=0;
    totZ=0;
    counter=0;
*/
}

void UserEventInformation::AddTotA(G4double a1, G4double a2, G4int copyNumber)
{
  totalAmplitude1Array[copyNumber] += a1;
  totalAmplitude2Array[copyNumber] += a2;  
  //  totA1 += a1;
  //  totA2 += a2;
}

G4double UserEventInformation::GetZRecon()
{
  totA1 = 0;
  totA2 = 0;
  return (1/(2*0.00447))*log(totA2/totA1); 
}

G4double UserEventInformation::CalculateZRecon(G4int copyNumber)
{
  // reconZ = (1/(2*0.00447))*log(totA2/totA1);
  G4double mu=attenuationCoefficients[copyNumber];
  zReconArray[copyNumber] = (1/(2*mu))*log(totalAmplitude2Array[copyNumber]/totalAmplitude1Array[copyNumber]);
  return zReconArray[copyNumber];
}
void UserEventInformation::AddToSetOfHitBars(G4int copyNumber)
{
  hitBarsSet.insert(copyNumber);
}
std::set<G4int> UserEventInformation::GetSetOfHitBars()
{
  return hitBarsSet;
}
void UserEventInformation::AddEdep(G4double edep,G4int copyNumber){
  totalEDeposition[copyNumber] += edep;
}
G4double UserEventInformation::GetEnergyDeposition(G4int copyNumber)
{
  return totalEDeposition[copyNumber];
}
void UserEventInformation::AddHitTime(G4double hitTime,G4int copyNumber)
{
  barHitTime[copyNumber]+=hitTime;
  stepCounter[copyNumber]+=1;
}
G4double UserEventInformation::GetHitTime(G4int copyNumber)
{
  return barHitTime[copyNumber];
}
G4int UserEventInformation::GetNumStepsInBar(G4int copyNumber)
{
  return stepCounter[copyNumber];
}
void UserEventInformation::AddTotZ(G4double zPos,G4int copyNumber)
{
  zPosition[copyNumber]+=zPos;
}
G4double UserEventInformation::GetTotZ(G4int copyNumber)
{
  return zPosition[copyNumber];
}
