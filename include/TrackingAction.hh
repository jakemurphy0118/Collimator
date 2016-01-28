#ifndef TrackingAction_h
#define TrackingAction_h 1

#include "G4UserTrackingAction.hh"
#include "globals.hh"
#include "TFile.h"
#include "TTree.h"
#include <fstream>

class TrackingAction : public G4UserTrackingAction {

  public:
    TrackingAction();
    virtual ~TrackingAction(){};
   
    virtual void PreUserTrackingAction(const G4Track*);
    virtual void PostUserTrackingAction(const G4Track*);

  
private:
  G4int nGammas;
  Double_t fTrackID ;
  Double_t fParentID ;
      
  Double_t fPreX ;
  Double_t fPreY ;
  Double_t fPreZ ;
  Double_t fPreLT ;
  Double_t fPreGT ;
  Double_t fPrePT ;
      
  Double_t fPostX ;
  Double_t fPostY ;
  Double_t fPostZ ;
  Double_t fPostLT ;
  Double_t fPostGT ;
  Double_t fPostPT ;

  Double_t fPDGMass    ;
  Double_t fPDGWidth    ;
  Double_t fPDGCharge   ;
  Double_t fPDGSpin   ;

  Double_t fPDGiSpin    ;
  Double_t fPDGiParity  ;
  Double_t fPDGiConjugation;
  Double_t fPDGIsospin  ;
  Double_t fPDGIsospin3 ;
  Double_t fPDGiIsospin ;
  Double_t fPDGiIsospin3  ;

  Double_t fPDGiGParity  ;
  Double_t fPDGMagneticMoment   ;
  Double_t fLeptonNumber;
  Double_t fBaryonNumber  ;
  Int_t fPDGEncoding ;

  Double_t fAtomicNumber;
  Double_t fAtomicMass  ;

  Double_t fVolume;
  Double_t fNextVolume;

  std::ofstream testFile;

};

#endif
