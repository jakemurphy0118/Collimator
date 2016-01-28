#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Trap.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4Color.hh"
#include "G4VisAttributes.hh"
#include "CylindricalSourceHolderConstruction.hh"
#include "SphericalSourceHolderConstruction.hh"
#include "AsymmetricalCylindricalSourceHolderConstruction.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction()
{
//  cylindricalHolder = new CylindricalSourceHolderConstruction();
  //  sphericalHolder = new SphericalSourceHolderConstruction();
  //  asymmetricalCylindricalHolder = new AsymmetricalCylindricalSourceHolderConstruction();
//  electromagnet = new Electromagnet();
  DefineMaterials();

  detectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}

void DetectorConstruction::DefineMaterials()
{

  G4double a; //atomic mass
  G4double z; //atomic number
  G4double density, fractionmass;
  G4int ncomponents,natoms;

  a = 22.989768*g/mole;
  G4Element* elNa = new G4Element("Sodium","Na",z=11.,a);
  a = 126.90447*g/mole;
  G4Element* elI  = new G4Element("Iodine","I",z=53.,a);
  a = 204.3833*g/mole;
  G4Element* elTl = new G4Element("Thallium","Tl",z=81.,a);
  a = 55.85*g/mole;
  G4Element* elFe = new G4Element("Iron" , "Fe" , z= 26., a);
  a = 12.011*g/mole;
  G4Element *elC = new G4Element("Carbon","C",z=6.,a);
  a = 58.9332*g/mole;
  G4Element* elCo = new G4Element( "Cobalt", "Co", 27. ,a);
  a = 16.00*g/mole;
  G4Element* elO = new G4Element("Oxygen","O",z=8.,a);
  a = 28.0855*g/mole;
  G4Element *elSi = new G4Element("Silicon","Si",z=14,a);
  a = 14.01*g/mole;
  G4Element* elN = new G4Element("Nitrogen","N",z=7.,a);

// NaI(Tl) material
  density=3.67*g/cm3;
  G4Material* NaI = new G4Material("NaI",density,ncomponents=2);
  NaI->AddElement(elNa,natoms=1);
  NaI->AddElement(elI,natoms=1);
  NaITl = new G4Material("NaI(Tl)",density,ncomponents=2);
  NaITl->AddMaterial(NaI,fractionmass=99.*perCent);
  NaITl->AddElement(elTl,fractionmass=1.*perCent);

  // Air
  density = 1.290*mg/cm3;
  Air = new G4Material("Air",density,ncomponents=2);
  Air->AddElement(elN, fractionmass=0.7);
  Air->AddElement(elO, fractionmass=0.3);

  // Stainless Steel
  ssteel = new G4Material
    ("Steel", density=7.7*g/cm3, ncomponents=3);
  ssteel->AddElement(elC, 0.04);
  ssteel->AddElement(elFe, 0.88);
  ssteel->AddElement(elCo, 0.08);

  // Quartz
  Quartz = new G4Material("Quartz", density= 2.66*g/cm3,
                                     ncomponents=2);
  Quartz->AddElement(elSi, 1);
  Quartz->AddElement(elO, 2);

}

//=========================================================

G4VPhysicalVolume* DetectorConstruction::Construct()
{
/*
  // First clean out the old geom**********************************************^^
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
*/
/*
  // Get nist material manager and SDManager
  G4NistManager* nist = G4NistManager::Instance();
*/

//-------------------------------------------------
  //------------------Volumes------------------------
  //-------------------------------------------------
  // z - Along Beam
  // x - Horizontal plane
  // y - vertical plane

  G4RotationMatrix noRotate;

  // Define some lengths
  LENALabz = 7200*mm;//beam axis
  LENALabx = 7200*mm;//perpindicular to beam axis
  LENALaby = 9800*mm;//vertical
  //Source_Detectorz = 10.525*mm; // from the drawings for the target holder
  // new source detector distances
  // Source_Detectorz = 5.525*mm;
// ----------WORLD-------------------------------------------------------
//-----------------------------------------------------------------------

  G4double HalfLENALabz = LENALabz*0.5;
  G4double HalfLENALabx = LENALabx*0.5;
  G4double HalfLENALaby = LENALaby*0.5;

  G4Box* LENALab_solid
    = new G4Box("LENALab",HalfLENALabx,HalfLENALaby,HalfLENALabz);
  LENALab_log =
    new G4LogicalVolume(LENALab_solid, Air, "LENALab_log", 0,0,0);
    //new G4LogicalVolume(LENALab_solid, vacuum, "LENALab_log", 0,0,0);
  //G4VPhysicalVolume* LENALab_phys = 
  LENALab_phys =
    new G4PVPlacement(0,                // no rotation
                      G4ThreeVector(),  // at (0,0,0)
                      LENALab_log,      // Logical Volume
                      "LENALab",        // Name
                      0,                // Mother Volume
                      false,            // no boolean operations
                      0);               // Copy number

//----------APPARATUS--------------------------------------------------
//---------------------------------------------------------------------

  //==============================================
  // Testing Different Source Holders
  //==============================================
/*
  //  G4cout << "SOURCEHOLDER CALLED THIS: " << sourceHolderName << G4endl;

  if(sourceHolderName=="Cylinder"){
    G4cout << "Source holder is cylinder" << G4endl;
    //    cylindricalHolder->BuildSource(LENALab_log);
  }
  else if(sourceHolderName=="Spherical"){
    G4cout << "Source holder is spherical" << G4endl;
    //    sphericalHolder->BuildSource(LENALab_log);
  }
  else if(sourceHolderName=="Cylinder_Asymmetrical"){
    G4cout << "Asymmetrical Cylindrical Source Holder" << G4endl;
    //    asymmetricalCylindricalHolder->BuildSource(LENALab_log);
  }
  else{
    G4cout << "Source holder not specified; using default" << G4endl;
    cylindricalHolder->BuildSource(LENALab_log);
    //    asymmetricalCylindricalHolder->BuildSource(LENALab_log);
  }
  */
  BuildAPEX();
  //  electromagnet->BuildCALIOPEMagnetLog(LENALab_log);
  //  new G4PVPlacement(0, G4ThreeVector(), magnet, "Magnet", LENALab_log, false, 0, true);
  LENALab_log ->SetVisAttributes(G4VisAttributes::Invisible);

//##############
  BuildCollimator();
//trying to put source container at middle 
	G4Box* container_box = new G4Box("container",5.*mm,5.*mm,0.45*mm);
	G4LogicalVolume* container_log = new G4LogicalVolume(container_box, Quartz, "container",0,0,0);
	new G4PVPlacement(0, G4ThreeVector(), container_log, "container", LENALab_log, false, 0);
//##############
  return LENALab_phys;

}


void DetectorConstruction::BuildAPEX()
{


//-----------Dimensions---------------------------------------
  G4double SlabLength = 57.0*cm;
  G4double SlabHeight = 0.8*cm;
  G4double SlabWidth = 7.1*cm;

  G4double HalfSlabLength = SlabLength*0.5;
  G4double HalfSlabHeight = SlabHeight*0.5;
  G4double HalfSlabWidth = SlabWidth*0.5;

  G4double TrapLength = 55.0*cm;
  G4double TrapHeight = 6.0*cm;
  G4double TrapWidthShort = 5.5*cm;
  G4double TrapWidthLong = 7.0*cm;

  G4double HalfTrapLength = TrapLength*0.5;
  G4double HalfTrapHeight = TrapHeight*0.5;
  G4double HalfTrapWidthShort = TrapWidthShort*0.5;
  G4double HalfTrapWidthLong = TrapWidthLong*0.5;

  G4double QuartzDia = 4.4*cm;
  G4double QuartzWidth = 1.1*cm;

  G4double QuartzRad = QuartzDia*0.5;
  G4double QuartzThick = QuartzWidth*0.5;

  G4double Radius = 24.5*cm;
  G4double NaIsegRotate;

  G4double Theta = 0.0*deg;
  G4double Phi = 0.0*deg;
  G4double Alpha = 0.0*deg;

  G4ThreeVector zAxis(0,0,1.);
  G4ThreeVector SSteelSlab(0,3.44*cm,0);
  G4ThreeVector QuartzHoleL(0,0,-28.05*cm);
  G4ThreeVector QuartzHoleR(0,0,28.05*cm);
  G4ThreeVector APEXradius(0,-24.5*cm,0);
  G4ThreeVector AngCorradius(0,0,24.5*cm);
  G4RotationMatrix noRotate;

//--------Geometry---------------------------------------------

G4Box* scintcasing_box
    = new G4Box("Slab Volume",HalfSlabWidth,HalfSlabHeight,HalfSlabLength+2.0*mm);

  G4Tubs* cyl_hole
    = new G4Tubs("Cylindrical Hole",0.*mm,QuartzRad,QuartzThick,0,twopi);

  G4Trap* ssteelcasing_trap
    = new G4Trap("Trap Volume",HalfTrapLength+12.0*mm,Theta,Phi,HalfTrapHeight+0.4*mm,
                HalfTrapWidthShort+0.4*mm,HalfTrapWidthLong+0.5*mm,
                Alpha,HalfTrapHeight+0.4*mm,HalfTrapWidthShort+0.4*mm,
                HalfTrapWidthLong+0.5*mm,Alpha);

  G4Box* NaIvolume_box
    = new G4Box("SSteel Slab",HalfSlabWidth,HalfSlabHeight,HalfSlabLength);

  G4Trap* NaIvolume_trap
    = new G4Trap("NaI Casing",HalfTrapLength+10.0*mm,Theta,Phi,HalfTrapHeight+0.4*mm,
                HalfTrapWidthShort+0.4*mm,HalfTrapWidthLong+0.5*mm,
                Alpha,HalfTrapHeight+0.4*mm,HalfTrapWidthShort+0.4*mm,
                HalfTrapWidthLong+0.5*mm,Alpha);

  // subtract material for quartz windows
  G4SubtractionSolid* ssteelcasing_trap_sub
    = new G4SubtractionSolid("Quartz Window Hole L",ssteelcasing_trap,cyl_hole,
                             &noRotate,QuartzHoleL);
  G4SubtractionSolid* scintcasing_trap
    = new G4SubtractionSolid("Quartz Window Hole R",ssteelcasing_trap_sub,cyl_hole,
                             &noRotate,QuartzHoleR);

  // now join the SSteelSlab to the main NaIcasing
  G4UnionSolid* scintcasing
    = new G4UnionSolid("Scintillator Casing",scintcasing_trap,scintcasing_box,
                       &noRotate,SSteelSlab);

  // now join the SlabVolume to the main TrapVolume
  G4UnionSolid* scintvolume
    = new G4UnionSolid("Scintillator Volume",NaIvolume_trap,NaIvolume_box,
                       &noRotate,SSteelSlab);



  G4LogicalVolume* ScintCasing_log
    = new G4LogicalVolume(scintcasing, ssteel, "SSteel NaI Casing", 0, 0, 0);

  G4LogicalVolume* ScintVolume_log
    = new G4LogicalVolume(scintvolume, NaITl /*was air*/, "NaI(Tl) Volume", 0, 0, 0);

  // place SSteel NaI Casing inside the Segment Volume
  new G4PVPlacement(0,
                        G4ThreeVector(),                // at (0,0,0)
                        ScintVolume_log,                // Logical Volume
                        "Segment Volume",               // Name
                        ScintCasing_log,                // Mother Volume
                        false,                          // no boolean operations
                        0);                             // Copy number
//FOR SINGLE BAR PLACEMENT
/*
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(180*deg,90*deg,90*deg),
                        AngCorradius),                  // at radius of angular correlation measurement
                        ScintCasing_log,                // Logical Volume
                        "NaI Segment 1",                // Name
                        LENALab_log,                    // Mother Volume
                        false,                          // no boolean operations
                        0);                             // Copy number
*/
//FOR REPEATED BARS
  for(G4int i=0;i<24;i++){
    NaIsegRotate=(i*(360/24))*deg;
    new G4PVPlacement(G4Transform3D(G4RotationMatrix(zAxis,-NaIsegRotate),
                      G4ThreeVector(Radius*std::sin(NaIsegRotate),Radius*std::cos(NaIsegRotate),0)),
                      ScintVolume_log,          // Logical Volume
                      "APEX Segment ",          // Name
                      LENALab_log,              // Mother Volume
                      false,                    // no boolean operations
                      i);                       // Copy number
  }

  //------------------Quartz Windows-----------------  

  G4Tubs* quartz_cyl
    = new G4Tubs("Quartz Cylinder",0.*mm,QuartzRad,QuartzThick,0,twopi);

  G4LogicalVolume* Quartz_log
    = new G4LogicalVolume(quartz_cyl, Quartz,"Quartz Window",0,0,0);

  new G4PVPlacement(0,                                          // no rotation
                        QuartzHoleL,                            // at (0,0,-28.05*cm)
                        Quartz_log,                             // Logical Volume
                        "Quartz Window L",                      // Name
                        ScintVolume_log,                        // Mother Volume
                        false,                                  // no boolean operations
                        0);                                     // Copy number

  new G4PVPlacement(0,                                          // no rotation
                        QuartzHoleR,                            // at (0,0,28.05*cm)
                        Quartz_log,                             // Logical Volume
                        "Quartz Window R",                      // Name
                        ScintVolume_log,                        // Mother Volume
                        false,                                  // no boolean operations
                        1);                                     // Copy number

  //------------------------------------------------------------
  //-------------- Visualization Stuff -----------------------
  //------------------------------------------------------------

  // The scintillator volume is invisible
//  ScintVolume_log ->SetVisAttributes(G4VisAttributes::Invisible);
  G4VisAttributes* NaIVisAtt = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
  ScintVolume_log->SetVisAttributes(NaIVisAtt);
  NaIVisAtt->SetForceSolid(true);

  // Scintillator casing is dark gray
  G4VisAttributes* ssteelVisAtt = new G4VisAttributes(G4Colour(0.4,0.4,0.4));
  ScintCasing_log->SetVisAttributes(ssteelVisAtt);
//  ssteelVisAtt->SetForceSolid(true);

/*//DOES NOT EXIST
  // NaI(Tl) scintillator is yellow
  G4VisAttributes* detVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  APEXNaI_log->SetVisAttributes(detVisAtt);
*/
  // Quartz window is white
  G4VisAttributes* quartzVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  Quartz_log->SetVisAttributes(quartzVisAtt);
//  quartzVisAtt->SetForceSolid(true);

}
/*
void DetectorConstruction::SetSourceHolder(G4String newValue)
{
  G4cout << "THE SOURCE HOLDER IS: " << newValue << G4endl;
  sourceHolderName=newValue;
  G4cout << "THE SOURCE HOLDER IS: " << sourceHolderName << G4endl;
}

*/
