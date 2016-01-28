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
// $Id: TrackInformation.cc,v 1.2 2006/06/29 17:44:25 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
//

#include "TrackInformation.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"

G4Allocator<TrackInformation> aTrackInformationAllocator;

using namespace CLHEP;

TrackInformation::TrackInformation()
{
    originalTrackID = 0;
    particleDefinition = 0;
    originalPosition = G4ThreeVector(0.,0.,0.);
    originalMomentum = G4ThreeVector(0.,0.,0.);
    originalEnergy = 0.;
    originalTime = 0.;
    sourceTrackID = -1;
    sourceTrackID = -1;
    sourceDefinition = 0;
    sourcePosition = G4ThreeVector(0.,0.,0.);
    sourceMomentum = G4ThreeVector(0.,0.,0.);
    sourceEnergy = 0.;
    sourceTime = 0.;

    sourceProcess = -1;
}

TrackInformation::TrackInformation(const G4Track* aTrack)
{
    originalTrackID = aTrack->GetTrackID();
    particleDefinition = aTrack->GetDefinition();
    originalPosition = aTrack->GetPosition();
    originalMomentum = aTrack->GetMomentum();
    originalEnergy = aTrack->GetTotalEnergy();
    originalTime = aTrack->GetGlobalTime();
    sourceTrackID = -1;
    sourceDefinition = 0;
    sourcePosition = G4ThreeVector(0.,0.,0.);
    sourceMomentum = G4ThreeVector(0.,0.,0.);
    sourceEnergy = 0.;
    sourceTime = 0.;

    sourceProcess = -1;
}

TrackInformation::TrackInformation(const TrackInformation* aTrackInfo)
{
    originalTrackID = aTrackInfo->originalTrackID;
    particleDefinition = aTrackInfo->particleDefinition;
    originalPosition = aTrackInfo->originalPosition;
    originalMomentum = aTrackInfo->originalMomentum;
    originalEnergy = aTrackInfo->originalEnergy;
    originalTime = aTrackInfo->originalTime;
    sourceTrackID = aTrackInfo->sourceTrackID;
    sourceDefinition = aTrackInfo->sourceDefinition;
    sourcePosition = aTrackInfo->sourcePosition;
    sourceMomentum = aTrackInfo->sourceMomentum;
    sourceEnergy = aTrackInfo->sourceEnergy;
    sourceTime = aTrackInfo->sourceTime;

    sourceProcess = aTrackInfo->sourceProcess;
}

TrackInformation::~TrackInformation()
{ ; }

TrackInformation& TrackInformation::operator =(const TrackInformation& aTrackInfo)
{

    return *this;
}

void TrackInformation::SetSourceTrackInformation(const G4Track* aTrack)
{
}

void TrackInformation::Print() const
{

}

