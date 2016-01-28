//----------------------------------------------------------------------
// AnalysisManager.cc
//
// Description: Sorts out the analysis
//                    At the moment, this is mainly used to handle user commands
//
// Version 0.1, 8/3/05
// Changes: Nothing!
//----------------------------------------------------------------------
//#ifdef G4ANALYSIS_USE

#include <fstream>

#include "G4ios.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4Track.hh"
#include "G4VVisManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"

#include "AnalysisManager.hh"
#include "AnalysisMessenger.hh"

using namespace CLHEP;

AnalysisManager* AnalysisManager::instance = 0;

AnalysisManager::AnalysisManager()
  :fileName("HPGe.root"),asciiFileBool(false),
   binWidth(1.0*CLHEP::keV),detailBool(false),resolutionBool(true),GeResolution(3.*CLHEP::keV),NaIResolution(0.07),SciResolution(0.14)
{
  // Define the Analysis Manager Messenger
  analysisMessenger = new AnalysisMessenger(this);

}

AnalysisManager::~AnalysisManager()
{ 
  delete analysisMessenger;
}

AnalysisManager* AnalysisManager::getInstance()
{
  if (instance == 0) instance = new AnalysisManager();
  return instance;
}


void AnalysisManager::dispose()
{
  if (instance != 0)
  {
    delete instance;
    instance = 0;
  }
}

//#endif // G4ANALYSIS_USE
