//--------------------------------------------------------------------
// EventAction.hh
//
// Description: Events! Include File.
// Changes: 7/18/05 None yet
//--------------------------------------------------------------------

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4Step.hh"
#include <vector>
#include <fstream>
#include "TimeDistribution.hh"
#include <TH1F.h>
class RunAction;
class PrimaryGeneratorAction;


class EventAction : public G4UserEventAction
{
public:
  EventAction(RunAction*,PrimaryGeneratorAction*);
  //EventAction(RunAction*,PrimaryGeneratorAction*,RootAnalysis*);
  ~EventAction();

public:
  virtual void BeginOfEventAction(const G4Event*);
  virtual void EndOfEventAction(const G4Event*);
  std::vector< std::vector<double> >*  GetGammaRayVector_WithSource();
  std::vector< std::vector<double> >*  GetGammaRayVector_NoSource();
  std::vector< std::vector<double> >*  GetGammaRayCoordinates();
  G4double* GetPsSubTypePtr();
  void SetPsSubTypePtr(G4double*);
  void SetGammaRayVector_WithSource(  std::vector< std::vector<double> >*);
  void SetGammaRayVector_NoSource(  std::vector< std::vector<double> >*);
  void SetGammaRayVector(int);
  void AlternativeReconstruction(std::vector< std::vector<double> >*,std::vector< std::vector<double> >*,std::string);//,int);
  void APEXReconstruction();
  void YamazakiReconstruction(std::vector<G4double>,std::vector<G4double>,std::vector<G4double>,std::vector<G4int>,G4double);
  void RyanReconstruction(std::vector<G4double>,std::vector<G4double>,std::vector<G4double>,std::vector<G4int>,G4double);
  TH1F* GetoPSTimingHistogramPtr();
  TH1F* GetpPSTimingHistogramPtr();
  TH1F* GetoPS_2GammaTimingHistogramPtr();
  TH1F* GetpPS_2GammaTimingHistogramPtr();
  TH1F* GetoPS_3GammaTimingHistogramPtr();
  TH1F* GetpPS_3GammaTimingHistogramPtr();  

  
private:
  G4int event_id;
  std::vector< std::vector<double> > gammaRayStorage_withSource;
  std::vector< std::vector<double> > gammaRayStorage_noSource;
  std::vector< std::vector<double> > gammaRayCoordinates;

  RunAction* runaction;
  PrimaryGeneratorAction* genaction;
  int num;

  //Angle Recon Information
  double n_x, n_y, n_z;
  double dir_x, dir_y, dir_z;
  double dir_x_1, dir_y_1, dir_z_1;
  double dir_x_2, dir_y_2, dir_z_2;
  double dir_x_3, dir_y_3, dir_z_3;
  double r;
  double theta;
  //  int N;
  TH1F* timingHistoPS;
  TH1F* timingHistpPS;
  TH1F* timingHistoPS_2Gamma;
  TH1F* timingHistpPS_2Gamma;
  TH1F* timingHistoPS_3Gamma;
  TH1F* timingHistpPS_3Gamma;
  std::ofstream eventFile;
  double zVertex;
  double PhiVertex;
  Double_t BRatio;
  double IDNum;
  G4double pSSubtype;
  public:


};

#endif
