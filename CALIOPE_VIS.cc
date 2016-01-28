/* Include Geant libraries and your own custom libraries (found in the include subdirectory for the CALIOPE-Simulation code.*/
#include "G4TrajectoryDrawByParticleID.hh"
#include "G4Geantino.hh"
#include <ctime>
#include "G4Polyline.hh"
#include "G4Colour.hh"
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4RunManager.hh"
#include "DetectorConstruction.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "G4VisExecutive.hh"
#include "G4ParticleGun.hh"
#include "Randomize.hh"
#include "AnalysisManager.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "PhysicsList.hh"
#include "SteppingAction.hh"
#include "TrackingAction.hh"
#include "StackingAction.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VUserPhysicsList.hh"
#include "QGSP_BERT_HP.hh"
#include "QGSP_BERT.hh"
#include "LBE.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4Colour.hh"
//Include a normal C++ library
#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

int main(int argc,char** argv)
{
  int counter=0;

  //Get user input for the spin
  if(argc==2){
    counter=atoi(argv[1]);
  }
  else{
    G4cout << "Please input a string for the output file." << G4endl;
    return 0;
  }

  //  std::cout << spinString << std::endl;
  // choose the Random engine                                                              
  time_t systime = time(NULL);
  long seed = (long) systime;
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  CLHEP::HepRandom::setTheSeed(seed);

  // Construct the run manager                                                             
  G4RunManager* runManager = new G4RunManager;

  // Set needed initialization classes                                                     

  //construct the analysis manager

  AnalysisManager* aMgr = AnalysisManager::getInstance();

  //Initialize the detector (geometry, etc)
  DetectorConstruction* detector = new DetectorConstruction;
  runManager->SetUserInitialization(detector);

  PhysicsList* newPhysicsList = new PhysicsList();
  // Intitialize the physics list (tells you what physics processes to run)
  runManager->SetUserInitialization(newPhysicsList);

  //Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  //Create new drawByParticleID model
  //  G4TrajectoryDrawByParticleID* model= new G4TrajectoryDrawByParticleID;
  //configure
  //model->SetDefault("magenta");
  //model->Set("gamma","magenta");
  //  visManager->RegisterModel(model);

  //Set User Action classes
  RunAction* runaction;
  EventAction* eventaction;
  PrimaryGeneratorAction* genaction;
  SteppingAction* steppingaction;
  TrackingAction* trackingaction;
  StackingAction* stackingaction;

  // runManager->SetUserAction(new PrimaryGeneratorAction(detector));
  //causes abort trap 6
  runManager->SetUserAction(genaction=new PrimaryGeneratorAction(detector));
  //00000000000

  runManager->SetUserAction(runaction = new RunAction(genaction,newPhysicsList));
  //  runManager->SetUserAction(runaction = new RunAction(genaction,rootanalysis));
  runManager->SetUserAction(eventaction =
			      new EventAction(runaction,genaction));
    //runManager->SetUserAction(eventaction =
    //		    new EventAction(runaction,genaction,rootanalysis));
  runManager->SetUserAction(stackingaction =
                            new StackingAction());
  runManager->SetUserAction(trackingaction =
                            new TrackingAction());
  runManager->SetUserAction(steppingaction =
			    new SteppingAction(detector, eventaction));

  //  runManager->SetUserInitialization(new PhysicsList);
  //Initialize the G4 kernel (General Initialization Step)
  runManager->Initialize();

  //****** VISUALIZATION ********//                                                                
  /*  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  G4Polyline x_axis;
  G4Polyline y_axis;
  G4Polyline z_axis;
  //Set red line color                                               
  G4Colour red(1.0,0.0,0.0);
  G4VisAttributes att(red);
  x_axis.push_back( G4Point3D (0.,0.,0.));
  x_axis.push_back( G4Point3D (0.*cm,500.*cm,0.));
  x_axis.SetVisAttributes(att);

  y_axis.push_back( G4Point3D (0.,0.,0.));
  y_axis.push_back( G4Point3D (500.*cm,0.*cm,0.));
  y_axis.SetVisAttributes(att);

  z_axis.push_back( G4Point3D (0.,0.,0.));
  z_axis.push_back( G4Point3D (0.*cm,0.*cm,500.*cm));
  z_axis.SetVisAttributes(att);*/

  //=====Visualization====//

  //Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  //Let's make it interactive mode for now
  G4UIExecutive* ui = new G4UIExecutive(argc, argv);

  //  UImanager->ApplyCommand("/vis~/clear/view");
  //  UImanager->ApplyCommand("/vis~/draw/current");
  UImanager->ApplyCommand("/control/execute vis.mac");
  //  pVVisManager->Draw(x_axis);
  //  pVVisManager->Draw(y_axis);
  //  pVVisManager->Draw(z_axis);
  //  UImanager->ApplyCommand("/vis~/show/view");

  //Need to pass a new filename each time
  std::stringstream ss;
  ss << "/analysis/filename " << counter << ".root";

  //  UImanager->ApplyCommand("/vis~/clear/view");
  //  UImanager->ApplyCommand("/vis~/draw/current");
  //  UImanager->ApplyCommand("/control/execute vis.mac");
  //  UImanager->ApplyCommand("/vis~/show/view");
  //  UImanager->ApplyCommand(ss.str().c_str());
  ui->SessionStart();

  //  UImanager->ApplyCommand("/control/execute macros/ion.mac");
  //ui->SessionStart();
  delete ui;
 
  //Job termination

  delete visManager;
  delete aMgr;
  //delete rootanalysis;
  delete detector;
  //delete model;
  delete genaction;
  delete runaction;
  delete eventaction;
  delete stackingaction;
  delete trackingaction;
  delete steppingaction;
  //delete runManager;
  return 0;

}
