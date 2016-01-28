//--------------------------------------------------------------------
// RunAction.hh
//
// Description: Defines the run itself. include file
// Changes: 7/18/05 None yet
//--------------------------------------------------------------------

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "TimeDistribution.hh"
#include <TH1F.h>
#include "PhysicsList.hh"
#include "G4VUserPhysicsList.hh"
#include "RunActionMessenger.hh"
class G4Run;
class PrimaryGeneratorAction;
class EventAction;

class RunAction : public G4UserRunAction
{
public:
  RunAction(PrimaryGeneratorAction*,PhysicsList*);
  //  RunAction(PrimaryGeneratorAction*,RootAnalysis*);
  ~RunAction();
  TH1F* GetoPSTimingHistogramPtr();
  TH1F* GetpPSTimingHistogramPtr();
  TH1F* GetoPS_2GammaTimingHistogramPtr();
  TH1F* GetpPS_2GammaTimingHistogramPtr();
  TH1F* GetoPS_3GammaTimingHistogramPtr();
  TH1F* GetpPS_3GammaTimingHistogramPtr();
  Double_t GetBranchingRatio();

public:
  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);
  void UpdateRunAction();
  void SetOutputFile(G4String);
  
private:
  PrimaryGeneratorAction* genAction;
  EventAction* evtAction;
  TimeDistribution* timedistr;
  TH1F* timingHistoPS;
  TH1F* timingHistpPS;
  TH1F* timingHistoPS_2Gamma;
  TH1F* timingHistpPS_2Gamma;
  TH1F* timingHistoPS_3Gamma;
  TH1F* timingHistpPS_3Gamma;

  //  TH1F* timingHist;
  Double_t BRatio;
  G4double BField;
  PhysicsList* myPhysListUpdate;
  const G4VUserPhysicsList* physListUpdate;
  //G4VModularPhysicsList* physListUpdate;
  RunActionMessenger* runActMessenger;
  G4String myPath;
};

#endif
