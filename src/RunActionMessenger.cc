#include "RunActionMessenger.hh"

#include "RunAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"

RunActionMessenger::RunActionMessenger(RunAction * run) : RunActionPtr(run)
{
  runActionDir = new G4UIdirectory("/runaction/");
  runActionDir->SetGuidance("run action commands");
  UpdateCmd = new G4UIcmdWithoutParameter("/runaction/update",this);
  UpdateCmd->SetGuidance("Update the run action");
  outputFileCmd = new G4UIcmdWithAString("/runaction/setOutputPath",this);
  outputFileCmd->SetGuidance("Update the output path");

}
RunActionMessenger::~RunActionMessenger()
{
  delete UpdateCmd;
  delete outputFileCmd;

}
void RunActionMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{

  if( command == UpdateCmd )
    {
      RunActionPtr->UpdateRunAction();
    }
  if( command == outputFileCmd )
    {
      RunActionPtr->SetOutputFile(newValue);
    }
}

