#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det) : Detector(Det)
{
  detDir = new G4UIdirectory("/det/");
  detDir->SetGuidance("Switch up detector components");
//  defineSourceHolder = new G4UIcmdWithAString("/det/defineSourceHolder",this);  
  UpdateCmd = new G4UIcmdWithoutParameter("/det/update",this);
  UpdateCmd->SetGuidance("Update geometry");


}
DetectorMessenger::~DetectorMessenger()
{
 // delete defineSourceHolder;
  delete UpdateCmd;

}
void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
 /* if( command == defineSourceHolder )
    { 
      Detector->SetSourceHolder(newValue);
    }
*/
  if( command == UpdateCmd )
    {
      Detector->UpdateGeometry();
    }
}

