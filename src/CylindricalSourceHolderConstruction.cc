//--------------------------------------------------------------------
// DetectorConstruction.cc
//
// Description: The detector definitions, materials etc.
// Changes: 7/15/05 None yet
//--------------------------------------------------------------------
#include "G4VVisManager.hh"
#include "G4Polyline.hh"
#include "G4Point3D.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Torus.hh"
#include "G4Polycone.hh"
#include "G4Ellipsoid.hh"
#include "G4Orb.hh"
#include "G4Trap.hh"
#include "G4Track.hh"
#include "G4Paraboloid.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4PVReplica.hh"
#include "G4RotationMatrix.hh"
#include "G4Trap.hh"
#include "G4Cons.hh"
#include "G4UnitsTable.hh"
#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "G4Navigator.hh"
#include "DetectorConstruction.hh"

using namespace CLHEP;

CylindricalSourceHolderConstruction::CylindricalSourceHolderConstruction()
{
  DefineMaterials();
}

CylindricalSourceHolderConstruction::~CylindricalSourceHolderConstruction()
{
}



void CylindricalSourceHolderConstruction::DefineMaterials()
{
  //---------------Materials-----------------

  G4NistManager* NISTman = G4NistManager::Instance();

  Al = new G4Material("Aluminum", z= 13., a= 26.98*g/mole, density= 2.7*g/cm3);

  // the elements
  a = 14.01*g/mole;
  G4Element* elN = new G4Element("Nitrogen","N",z=7.,a);
  a = 16.00*g/mole;
  G4Element* elO = new G4Element("Oxygen","O",z=8.,a);
  a = 12.011*g/mole;
  G4Element *elC = new G4Element("Carbon","C",z=6.,a);
  a = 1.00794*g/mole;
  G4Element* elH = new G4Element("Hydrogen","H",z=1.,a);
  a = 22.989768*g/mole;
  G4Element* elNa = new G4Element("Sodium","Na",z=11.,a);
  a = 126.90447*g/mole;
  G4Element* elI  = new G4Element("Iodine","I",z=53.,a);
  a = 204.3833*g/mole;
  G4Element* elTl = new G4Element("Thallium","Tl",z=81.,a);
  a = 58.9332*g/mole;
  G4Element* elAl = new G4Element("Aluminium","Al",z=13.,a);
  a = 58.69*g/mole;
  G4Element *elNi = new G4Element("Nickel","Ni",z=28,a);
  a = 95.96*g/mole;
  G4Element *elSi = new G4Element("Silicon","Si",z=14,a);
  a = 30.9738*g/mole;

  // Air
  density = 1.290*mg/cm3;
  Air = new G4Material("Air",density,ncomponents=2);
  Air->AddElement(elN, fractionmass=0.7);
  Air->AddElement(elO, fractionmass=0.3);

  // Water
  Water = new G4Material("Water", density= 1.0*g/cm3,ncomponents=2);
  Water->AddElement(elH, 2);
  Water->AddElement(elO, 1);
  
  // Scintillator BC-408 Material
  Sci = new G4Material("Scintillator", density= 1.032*g/cm3, ncomponents=2);
  Sci->AddElement(elC, natoms=10);
  Sci->AddElement(elH, natoms=11);

  //Silicon Dioxide (SiO2)
  density = 2.648*g/cm3;
  SiO2 = new G4Material("silicon dioxide",density,ncomponents=2);
  SiO2->AddElement(elSi,0.47);
  SiO2->AddElement(elO,0.53);

    // define a material from elements and/or others materials (mixture of mixtures)
  density = 0.010*g/cm3;
  Aerog = new G4Material(name="Aerogel", density, ncomponents=3);
  Aerog->AddMaterial(SiO2, fractionmass=62.5*perCent);
  Aerog->AddMaterial(Water , fractionmass=37.4*perCent);
  Aerog->AddElement (elC , fractionmass= 0.1*perCent);
  
  // Kapton 
  Kapton  = NISTman->FindOrBuildMaterial("G4_KAPTON");

  // Delrin
  delrin = NISTman->FindOrBuildMaterial("G4_POLYOXYMETHYLENE");

  // Anthracene
  anthracene  = NISTman->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

  // Vacuum
  G4double atomicNumber = 1.;
  G4double massOfMole = 1.008*g/mole;
  density = 1.e-25*g/cm3;
  G4double temperature = 77.*kelvin;
  G4double pressure = 0.00133*pascal;
  vacuum =
    new G4Material("Vacuum", atomicNumber,
		   massOfMole, density, kStateGas,
		   temperature, pressure);

}

void CylindricalSourceHolderConstruction::BuildSource(G4LogicalVolume* LENALabVolume){

  G4LogicalVolume* vacuumVolume_log=BuildVacuum(LENALabVolume);

  G4LogicalVolume* aerogelVolume_log=BuildAerogel(vacuumVolume_log);

  //  G4LogicalVolume* kapton_log=BuildKapton(vacuumVolume_log);

  //  G4LogicalVolume* scint_log=BuildScintillator(vacuumVolume_log);

}

G4LogicalVolume* CylindricalSourceHolderConstruction::BuildVacuum(G4LogicalVolume* LENALabVolume){

  //===================================
  // Define Parameters
  //===================================
  density=1.e-25*g/cm3;
  G4double atomicNumber=1.;
  G4double massOfMole=1.008*g/mole;
  G4double temperature=77.*kelvin;
  G4double pressure=0.00133*pascal;
  vacuum=new G4Material("Vacuum",atomicNumber,massOfMole,density,kStateGas,temperature,pressure);

  //===================================
  // Define Dimensions
  //===================================

  G4double vacuumVolumeRadius=4.4501*cm;
  //G4double vacuumVolumeRadius=7.0*cm;
  // This is half the length of the cylinder
  G4double vacuumVolumeLength=5.051*cm;
  //G4double vacuumVolumeLength=12.0*cm;
  G4Tubs* vacuumVolume = new G4Tubs("vacuumVolumeExterior",0,vacuumVolumeRadius,vacuumVolumeLength,0,2*pi);
  G4LogicalVolume* vacuumVolume_log = new G4LogicalVolume(vacuumVolume,vacuum, "vacuumVolumeExt");


  //===================================
  // Place Vacuum
  //===================================

  new G4PVPlacement(0,
		   G4ThreeVector(0,0,0),
		   vacuumVolume_log,
		   "Source",
		   LENALabVolume,
		   false,
		   0);

  //===================================
  // Set Visibility of Vacuum
  //===================================

  G4VisAttributes* VacVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,0.8));
  VacVisAtt->SetForceSolid(1);
  VacVisAtt->SetVisibility(1);
  vacuumVolume_log->SetVisAttributes(VacVisAtt);

  return vacuumVolume_log;

}

G4LogicalVolume* CylindricalSourceHolderConstruction::BuildAerogel(G4LogicalVolume* vacuumVolume_log){


  //Definition of Aerogel

  G4Tubs* aeroDisk = new G4Tubs("aeroDiskSolid",0.0*cm,4.21*cm,2.3*cm,0,2*pi);

  //  G4Orb* aeroSphere = new G4Orb("aeroSphere",2.3*cm);

  G4LogicalVolume* aeroDisk_log = new G4LogicalVolume(aeroDisk,Aerog, "Aerog log");

  //  G4LogicalVolume* aeroSphere_log =
  //new G4LogicalVolume(aeroSphere,Aerog, "Aerog log");

  double halfAerogelLength=2.3001*cm;
  double halfAerogelGap=1*mm;
  double aerogelShift=halfAerogelLength+halfAerogelGap;

  /*  new G4PVPlacement(0,
		    G4ThreeVector(0,0,0),
		    aeroSphere_log,
		    "Source",
		    vacuumVolume_log,
		    false,
		    0);*/

  new G4PVPlacement(0,
		    G4ThreeVector(0,0,aerogelShift),
		    aeroDisk_log,
		    "Source",
		    vacuumVolume_log,
		    false,
		    0);

  new G4PVPlacement(0,
		    G4ThreeVector(0,0,-aerogelShift),
		    aeroDisk_log,
		    "Source",
		    vacuumVolume_log,
		    false,
		    0);
  //Aero Color
  G4VisAttributes* AeroColor = new G4VisAttributes(G4Colour(0.0,0.8,0.8));
  AeroColor->SetForceSolid(1);
  AeroColor->SetVisibility(1);
  aeroDisk_log->SetVisAttributes(AeroColor);

}

G4LogicalVolume* CylindricalSourceHolderConstruction::BuildKapton(G4LogicalVolume* vacuumVolume_log){

  double halfKaptonWidth=0.4*mm;
  double innerFoilRadius=0*mm;
  double outerFoilRadius=31*mm;
  G4Tubs* kaptonFoil = new G4Tubs("kaptonFoilSolid",innerFoilRadius,outerFoilRadius,halfKaptonWidth,0,2*pi);

  G4LogicalVolume* kaptonFoil_log =
    new G4LogicalVolume(kaptonFoil,Kapton, "kapton log");

  new G4PVPlacement(0,
		    G4ThreeVector(0,0,0),
		    kaptonFoil_log,
		    "Source",
		    vacuumVolume_log,
		    false,
		    0);
  
  // Kapton Color is Green
  G4VisAttributes* KapVisAtt = new G4VisAttributes(G4Colour(0.8,0.0,0.8));
  KapVisAtt->SetForceSolid(1);
  kaptonFoil_log->SetVisAttributes(KapVisAtt);

  return vacuumVolume_log;

}

G4LogicalVolume* CylindricalSourceHolderConstruction::BuildScintillator(G4LogicalVolume* vacuumVolume_log){

  //Scintillator
  double halfScintWidth=0.15*mm;
  double halfScintGap=0.2*mm;
  double scintShift=halfScintWidth+halfScintGap;

  G4double outerScintRadius = 4.21*cm;
  G4double innerScintRadius = 0.0*mm;

  G4Tubs* scintDisk = new G4Tubs("scintDiskSolid",innerScintRadius,outerScintRadius,halfScintWidth,0,2*pi);

  G4LogicalVolume* scintDisk_log =
    new G4LogicalVolume(scintDisk,anthracene, "scint log");

  new G4PVPlacement(0,
		    G4ThreeVector(0,0,scintShift),
		    scintDisk_log,
		    "Source",
		    vacuumVolume_log,
		    false,
		    0);
  
  new G4PVPlacement(0,
		    G4ThreeVector(0,0,-scintShift),
		    scintDisk_log,
		    "Source",
		    vacuumVolume_log,
		    false,
		    0);

  // Scintillator Color is Green
  G4VisAttributes* ScintVisAtt = new G4VisAttributes(G4Colour(0.0,0.8,0.0));
  ScintVisAtt->SetForceSolid(1);
  scintDisk_log->SetVisAttributes(ScintVisAtt);
  return vacuumVolume_log;

}

void CylindricalSourceHolderConstruction::BuildShell(G4LogicalVolume* vacuumVolume_log){

  double pipeLength=30.0*cm;
  double halfGap=30.0*cm;
  double halfPipeLength=pipeLength/2.0;
  double shiftPipe=0;
  double pipeRadius=1.0*cm;

  G4double fractionmass;
  G4int ncomponents,natoms;


  G4Tubs* sourcePipeSolid = new G4Tubs("sourcePipeSolid",0,pipeRadius,halfPipeLength,0,360);

  G4LogicalVolume* sourcePipeSolid_log =
    new G4LogicalVolume(sourcePipeSolid,Al, "sourcePipe log");

  /*  new G4PVPlacement(0,
		    G4ThreeVector(0,0,shiftPipe),
		    sourcePipeSolid_log,
		    "Source",
		    LENALabVolume,
		    false,
		    0);


  new G4PVPlacement(0,
		    G4ThreeVector(0,0,-shiftPipe),
		    sourcePipeSolid_log,
		    "Source",
		    LENALabVolume,
		    false,
		    0);*/

  // Rod color is blue
  G4VisAttributes* RodVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,0.8));
  RodVisAtt->SetForceSolid(1);
  sourcePipeSolid_log->SetVisAttributes(RodVisAtt);

  G4RotationMatrix noRotate;

  //==================================
  // Source Tube Parts
  //==================================

  //Below is the original inner source tube radius
  //  G4double innerSourceTubeRadius=4.2*cm;
  G4double innerSourceTubeRadius=0*cm;
  G4double outerSourceTubeRadius=4.3*cm;
  G4double sourceTubeLength=1.0*cm;
  //Hollow source container
  G4Tubs* sourceTube_1 = new G4Tubs("sourcePipeSolid",innerSourceTubeRadius,outerSourceTubeRadius,sourceTubeLength,0,2*pi);
  G4Tubs* sourceTube_0 = new G4Tubs("sourcePipeSolid",innerSourceTubeRadius,4.2*cm,sourceTubeLength,0,2*pi);
  G4Tubs* sourceTube_2 = new G4Tubs("sourcePipeSolid",innerSourceTubeRadius,outerSourceTubeRadius,sourceTubeLength,0,2*pi);

  G4RotationMatrix noRotation;
  noRotation.rotateY(0);
  G4ThreeVector translation(0,0,50);
  G4SubtractionSolid sourceTube("sourceTube_1-sourceTube_0",sourceTube_1,sourceTube_0,&noRotation,translation);

  G4LogicalVolume* sourceTube_log_1 =
    new G4LogicalVolume(sourceTube_1,Al, "sourcePipe log");

  G4LogicalVolume* sourceTube_log =
    new G4LogicalVolume(&sourceTube,Al, "sourcePipe log");

  G4LogicalVolume* sourceTube_log_2 =
    new G4LogicalVolume(sourceTube_2,Al, "sourcePipe log");


  new G4PVPlacement(0,
		    G4ThreeVector(0,0,-4*cm),
		    sourceTube_log_1,
		    "Source",
		    vacuumVolume_log,
		    false,
		    0);


  new G4PVPlacement(0,
		    G4ThreeVector(0,0,4*cm),
		    sourceTube_log_2,
		    "Source",
		    vacuumVolume_log,
		    false,
		    0);

  // Aluminium Tube Color is red
  G4VisAttributes* TubeVisAtt = new G4VisAttributes(G4Colour(0.8,0.0,0.0));
  TubeVisAtt->SetForceSolid(1);
  TubeVisAtt->SetVisibility(1);
  sourceTube_log_1->SetVisAttributes(TubeVisAtt);

  sourceTube_log_2->SetVisAttributes(TubeVisAtt);
  //=============================================

  //=============================================
  // Source Endcaps
  //=============================================

  G4double innerSourceEndcapRadius=0*cm;
  G4double outerSourceEndcapRadius=outerSourceTubeRadius;
  G4double sourceEndcapLength=0.5*mm;

  G4Tubs* sourceTubeEndcap_1 = new G4Tubs("sourcePipeSolid",innerSourceEndcapRadius,outerSourceEndcapRadius,sourceEndcapLength,0,2*pi);
  G4Tubs* sourceTubeEndcap_2 = new G4Tubs("sourcePipeSolid",innerSourceEndcapRadius,outerSourceEndcapRadius,sourceEndcapLength,0,2*pi);

  G4LogicalVolume* sourceTubeEndcap_log_1 =
    new G4LogicalVolume(sourceTubeEndcap_1,Al, "sourcePipe log");

  G4LogicalVolume* sourceTubeEndcap_log_2 =
    new G4LogicalVolume(sourceTubeEndcap_2,Al, "sourcePipe log");

  new G4PVPlacement(0,
		    G4ThreeVector(0,0,-5.0*cm),
		    sourceTubeEndcap_log_1,
		    "Source",
		    vacuumVolume_log,
		    false,
		    0);

  new G4PVPlacement(0,
		    G4ThreeVector(0,0,5.0*cm),
		    sourceTubeEndcap_log_2,
		    "Source",
		    vacuumVolume_log,
		    false,
		    0);

  sourceTubeEndcap_log_1->SetVisAttributes(TubeVisAtt);

  sourceTubeEndcap_log_2->SetVisAttributes(TubeVisAtt);
  
  //===========================================================

  //===========================================================
  // Source Tube Raised Section (to match with threaded section)
  //===========================================================

  //Set the inside radius to the same as the rest
  G4double innerTubeStepRadius=innerSourceTubeRadius;

  //Wider radius so that it can attach to the threaded section
  G4double outerTubeStepRadius=outerSourceTubeRadius+1.5*mm;

  //Length of the raised step
  G4double tubeStepLength=1.0*cm;


  G4Tubs* sourceTubeMiddle_1 = new G4Tubs("sourcePipeSolid",innerTubeStepRadius,outerTubeStepRadius,tubeStepLength,0,2*pi);
  G4Tubs* sourceTubeMiddle_2 = new G4Tubs("sourcePipeSolid",innerSourceTubeRadius,outerTubeStepRadius,tubeStepLength,0,2*pi);

  G4LogicalVolume* sourceTubeMiddle_log_1 =
    new G4LogicalVolume(sourceTubeMiddle_1,Al, "sourcePipe log");

  G4LogicalVolume* sourceTubeMiddle_log_2 =
    new G4LogicalVolume(sourceTubeMiddle_2,Al, "sourcePipe log");

  new G4PVPlacement(0,
		    G4ThreeVector(0,0,-2*cm),
		    sourceTubeMiddle_log_1,
		    "Source",
		    vacuumVolume_log,
		    false,
		    0);
  new G4PVPlacement(0,
		    G4ThreeVector(0,0,2*cm),
		    sourceTubeMiddle_log_2,
		    "Source",
		    vacuumVolume_log,
		    false,
		    0);


  sourceTubeMiddle_log_1->SetVisAttributes(TubeVisAtt);
  sourceTubeMiddle_log_2->SetVisAttributes(TubeVisAtt);

  //============================================================

  //============================================================
  // Source Tube Threading
  //============================================================

  G4double threadingWidth=(outerSourceTubeRadius-innerSourceTubeRadius)/2;

  G4double innerThreadingRadius_1=innerSourceTubeRadius+threadingWidth;

  G4double outerThreadingRadius_1=outerTubeStepRadius;

  G4double innerThreadingRadius_2=innerSourceTubeRadius;

  G4double outerThreadingRadius_2=innerSourceTubeRadius+threadingWidth;

  G4double threadingLength=1.0*cm;

  G4Tubs* sourceTubeThreading_1 = new G4Tubs("sourcePipeSolid",innerThreadingRadius_1,outerThreadingRadius_1,threadingLength+0.25*cm,0,2*pi);

  G4Tubs* sourceTubeThreading_2 = new G4Tubs("sourcePipeSolid",innerThreadingRadius_2,outerThreadingRadius_2,threadingLength,0,2*pi);

  G4LogicalVolume* sourceTubeThreading_log_1 =
    new G4LogicalVolume(sourceTubeThreading_1,Al, "sourcePipe log");

  G4LogicalVolume* sourceTubeThreading_log_2 =
    new G4LogicalVolume(sourceTubeThreading_2,Al, "sourcePipe log");

  //The threaded sections sit on top of each other

  //Scoot this part back a bit for the oring
  new G4PVPlacement(0,
		    G4ThreeVector(0,0,0*cm),
		    sourceTubeThreading_log_1,
		    "Source",
		    vacuumVolume_log,
		    false,
		    0);

  new G4PVPlacement(0,
		    G4ThreeVector(0,0,0*cm),
		    sourceTubeThreading_log_2,
		    "Source",
		    vacuumVolume_log,
		    false,
		    0);

  G4VisAttributes* Thread1VisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  G4VisAttributes* Thread2VisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  Thread1VisAtt->SetForceSolid(1);
  Thread2VisAtt->SetForceSolid(1);
  Thread1VisAtt->SetVisibility(1);
  Thread2VisAtt->SetVisibility(1);
  sourceTubeThreading_log_1->SetVisAttributes(Thread1VisAtt);
  sourceTubeThreading_log_2->SetVisAttributes(Thread2VisAtt);

  //==============================

  G4NistManager* NISTman = G4NistManager::Instance();
  rubber  = NISTman->FindOrBuildMaterial("G4_RUBBER_NATURAL");

  G4double innerOringRadius=innerThreadingRadius_2;
  G4double outerOringRadius=outerThreadingRadius_2;
  G4double oringLength=0.25*cm;

  G4Tubs* sourceTube_Oring = new G4Tubs("sourcePipeSolid",innerOringRadius,outerOringRadius,oringLength,0,2*pi);

  G4LogicalVolume* sourceTube_Oring_log =
    new G4LogicalVolume(sourceTube_Oring,rubber, "sourcePipe log");

  /*  new G4PVPlacement(0,
		    G4ThreeVector(0,0,-1*cm),
		    sourceTube_Oring_log,
		    "Source",
		    vacuumVolume_log,
		    false,
		    0);*/

  G4VisAttributes* OringVisAtt = new G4VisAttributes(G4Colour(0.0,0.8,0.0));
  OringVisAtt->SetForceSolid(1);

  sourceTube_Oring_log->SetVisAttributes(OringVisAtt);


}
