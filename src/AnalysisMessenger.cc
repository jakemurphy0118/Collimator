//----------------------------------------------------------------------
// AnalysisMessenger.cc
//
// Description: Sends info to the Analysis Manager
//                     Most generic comands are here too
//
// Version 0.1, 1/18/06
// Changes: Nothing!
//----------------------------------------------------------------------
//#ifdef G4ANALYSIS_USE

//#include <AIDA/AIDA.h>

#include "AnalysisMessenger.hh"
#include "AnalysisManager.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"

#include "G4ios.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4Track.hh"
#include "G4VVisManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4String.hh"

AnalysisMessenger::AnalysisMessenger(AnalysisManager* analysisManager)
  :Analysis(analysisManager)

{
  AnalysisDir = new G4UIdirectory("/analysis/");

  //  G4String analysiscntrl="analysis control"
  //AnalysisDir->SetGuidance(analysiscntrl);

  FileNameCmd = new G4UIcmdWithAString("/analysis/filename",this);
  //  FileNameCmd->SetGuidance("Input the name for the AIDA output file.");
  FileNameCmd->SetParameterName("filename",true,true);
  FileNameCmd->SetDefaultValue("LENAGe.aida");

  AsciiFileCmd = new G4UIcmdWithABool("/analysis/asciifile",this);
  //  AsciiFileCmd->SetGuidance("Set to true if you want an ascii spectrum exported.");
  AsciiFileCmd->SetParameterName("asciiFileBool",true,false);
  AsciiFileCmd->SetDefaultValue(false);

  BinWidthCmd = new G4UIcmdWithADouble("/analysis/binwidth",this);
  //  BinWidthCmd->SetGuidance("Set the width for ascii spectrum bins.");
  BinWidthCmd->SetParameterName("binwidth",true,false);
  BinWidthCmd->SetDefaultValue(1.*(CLHEP::keV));

  DetailCmd = new G4UIcmdWithABool("/analysis/detail",this);
  //  DetailCmd->SetGuidance("Set to true if you want detailed aida stepping information.");
  DetailCmd->SetParameterName("detailBool",true,false);
  DetailCmd->SetDefaultValue(false);

  ResolutionCmd = new G4UIcmdWithABool("/analysis/resolution",this);
  //  ResolutionCmd->SetGuidance("Set to true if you want to use detector resolution.");
  ResolutionCmd->SetParameterName("resolutionBool",true,false);
  ResolutionCmd->SetDefaultValue(false);

  GeResolutionCmd = new G4UIcmdWithADouble("/analysis/geresolution",this);
  //  GeResolutionCmd->SetGuidance("Set the HPGe resolution (full width half max) in keV.");
  GeResolutionCmd->SetParameterName("GeResolution",true,false);
  GeResolutionCmd->SetDefaultValue(3.);

  NaIResolutionCmd = new G4UIcmdWithADouble("/analysis/nairesolution",this);
  //  NaIResolutionCmd->SetGuidance("Set the NaI resolution (%E full width half max)");
  NaIResolutionCmd->SetParameterName("NaIResolution",true,false);
  NaIResolutionCmd->SetDefaultValue(0.07);

  NaIResolutionCmd = new G4UIcmdWithADouble("/analysis/sciresolution",this);
  //  NaIResolutionCmd->SetGuidance("Set the Scintillator resolution (%E full width half max)");
  NaIResolutionCmd->SetParameterName("SciResolution",true,false);
  NaIResolutionCmd->SetDefaultValue(0.14);

  NaIPosResolutionCmd = new G4UIcmdWithADouble("/analysis/naiposresolution",this);
  //  NaIPosResolutionCmd->SetGuidance("Set the Position resolution (%E full width half max)");
  NaIPosResolutionCmd->SetParameterName("PositionResolution",true,false);
  NaIPosResolutionCmd->SetDefaultValue(0.14);
}

AnalysisMessenger::~AnalysisMessenger()
{
  delete FileNameCmd; 
  delete AsciiFileCmd;
  delete BinWidthCmd;
  delete DetailCmd;
  delete ResolutionCmd;
  delete GeResolutionCmd;
  delete NaIResolutionCmd;
}

void AnalysisMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 

  // Where the commands are executed. 
  if( command == FileNameCmd ) { 
    Analysis->SetFileName(newValue);
  }
  if( command == AsciiFileCmd ) {
    Analysis->SetasciiFile(AsciiFileCmd->GetNewBoolValue(newValue));
  }
  if( command == BinWidthCmd ) {
    Analysis->SetBinWidth(BinWidthCmd->GetNewDoubleValue(newValue));
  }
  if( command == DetailCmd ) {
    Analysis->SetDetailBool(DetailCmd->GetNewBoolValue(newValue));
  }
  if( command == ResolutionCmd ) {
    Analysis->SetResolutionBool(ResolutionCmd->GetNewBoolValue(newValue));
  }
  if( command == GeResolutionCmd ) {
    Analysis->SetGeResolution(GeResolutionCmd->GetNewDoubleValue(newValue));
  }
  if( command == NaIResolutionCmd ) {
    Analysis->SetNaIResolution(NaIResolutionCmd->GetNewDoubleValue(newValue));
  }
  if( command == SciResolutionCmd ) {
    Analysis->SetSciResolution(SciResolutionCmd->GetNewDoubleValue(newValue));
  }
}

//#endif // G4ANALYSIS_USE
