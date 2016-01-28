#include "TrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "TrackInformation.hh"
#include "G4VProcess.hh"
#include "G4ParticleDefinition.hh"
#include "G4RunManager.hh"
#include <fstream>
#include <iostream>
#include "G4ios.hh"

TrackingAction::TrackingAction() : fTrackID(0), fParentID(0), fPreX(0), fPreY(0), fPreZ(0), fPreLT(0), fPreGT(0), fPrePT(0), fPostX(0), fPostY(0), fPostZ(0), fPostLT(0), fPostGT(0), fPostPT(0), fPDGMass(0), fPDGWidth(0), fPDGCharge(0), fPDGSpin(0), fPDGiSpin(0), fPDGiParity(0), fPDGiConjugation(0), fPDGIsospin(0), fPDGIsospin3(0), fPDGiIsospin(0), fPDGiIsospin3(0), fPDGiGParity(0), fPDGMagneticMoment(0), fLeptonNumber(0), fBaryonNumber(0), fPDGEncoding(0), fAtomicNumber(0), fAtomicMass(0),fVolume(0), fNextVolume(0)
{}


void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{ 
  nGammas=0; 
  // Reset products of radioactive decay to start at time zero.
  // Otherwise every track may start at 10^10 years in the future.
  //  if( fParentID == 1 && aTrack->GetGlobalTime() > 1 * second )
  
  /*  if( (fParentID==0) && (aTrack->GetGlobalTime() > 1 * second ))
    {

      (const_cast<G4Track *>(aTrack))->SetGlobalTime( 0. );
      }*/
    
  //  testFile.open("test",ios::app);

}

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  
}


