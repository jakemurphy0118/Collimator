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
// $Id: StackingAction.cc,v 1.3 2010/11/24 22:47:23 asaim Exp $
// GEANT4 tag $Name: geant4-09-04 $
//
#include <TROOT.h>
#include <TH1F.h>
#include "StackingAction.hh"
#include "G4VProcess.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ios.hh"
#include <fstream>
#include "TrackInformation.hh"
//#include "TimeDistribution.hh"
//#include "CalorimeterHit.hh"
using namespace std;

StackingAction::StackingAction()
{ 

}

StackingAction::~StackingAction()
{ 
  //  delete timedistr;
}

G4ClassificationOfNewTrack 
StackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
  
  G4ClassificationOfNewTrack classification = fUrgent;

  TrackInformation* trackInfo;
  if(aTrack->GetParentID()==0){
    trackInfo = new TrackInformation(aTrack);
    G4Track* theTrack = (G4Track*)aTrack;
    theTrack->SetUserInformation(trackInfo);
  }else{
    trackInfo = (TrackInformation*)(aTrack->GetUserInformation());
   
  }

  const G4VProcess* creatorProcess = aTrack->GetCreatorProcess();
  
  if(creatorProcess!=NULL){
  G4String processName = creatorProcess->GetProcessName();
  G4int subType = creatorProcess->GetProcessSubType();
  
  //  timedistr = new TimeDistribution();

  if(subType==26){
    //TH1F* timeHist = timedistr->buildHist("oPSDecay");
    //Double_t timeRoot = timeHist->GetRandom();
    //    G4double time = timeRoot*ns;
    G4double time = sameTime;
    (const_cast<G4Track *>(aTrack))->SetGlobalTime(aTrack->GetGlobalTime() - trigger + time + aTrack->GetGlobalTime());
    //    testFile << aTrack->GetGlobalTime() << std::endl;

  }

  //0.125ns
  if(subType==25){
    //    TH1F* timeHist = timedistr->buildHist("pPSDecay");
    //Double_t timeRoot = timeHist->GetRandom();
    //G4double time = timeRoot*ns;
    G4double time = sameTime;
    (const_cast<G4Track *>(aTrack))->SetGlobalTime(aTrack->GetGlobalTime() - trigger + time + aTrack->GetGlobalTime());
    //    testFile << aTrack->GetGlobalTime() << std::endl;
  }
  }

  return classification;
}

void StackingAction::NewStage()
{ ; }
    
void StackingAction::PrepareNewEvent()
{ ; }
void StackingAction::SetTrigger(G4double trig)
{ 
  trigger=trig;
}
void StackingAction::SetSameTime(G4double same)
{ 
  sameTime=same;
}

