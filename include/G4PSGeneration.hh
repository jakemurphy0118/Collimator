#ifndef G4PSGeneration_h
#define G4PSGeneration_h 1

#include "G4VEmProcess.hh"
#include "PsGammaDecayGenerator.hh"
#include "G4Positron.hh"
#include "G4VEmModel.hh"
#include <fstream>
//#include "ThetaDistribution.hh"
#include <TH3F.h>
#include <TH1F.h>
#include <TFile.h>
#include "globals.hh"
//#include "G4Positronium.hh"

class G4ParticleDefinition;

class G4PSGeneration : public G4VEmProcess
{

public:

  G4PSGeneration(const G4String& name = "PSDecay",double H=0);

  virtual ~G4PSGeneration();
  
  virtual G4bool IsApplicable(const G4ParticleDefinition& p);

  virtual G4VParticleChange* AtRestDoIt(const G4Track& track,
					const G4Step& stepData);

  virtual G4double AtRestGetPhysicalInteractionLength(const G4Track& track,
						      G4ForceCondition* condition);

  //Print out of the class parameters
  virtual void PrintInfo();

  G4VParticleChange* Do2GammaDecay(const G4Track& track, const G4Step& stepData, G4int, G4bool, G4String );
  G4VParticleChange* Do3GammaDecay(const G4Track& track, const G4Step& stepData, G4int, G4bool, G4String);

protected:

  virtual void InitialiseProcess(const G4ParticleDefinition*);
  
private:

  double H;
  G4bool isInitialised;
  const G4ParticleDefinition* theGamma;
  std::ofstream check;
  //  ThetaDistribution* thetadistr;
  PsGammaDecayGenerator* gen;
  std::ofstream checkDirFile;
  TH3F* spinHist;
  ofstream outfile;
  double BR_pPS;
  double BR_oPS;

};
#endif
