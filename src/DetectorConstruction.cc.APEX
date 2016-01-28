//--------------------------------------------------------------------
// DetectorConstruction.cc
//
// Description: The detector definitions, materials etc.
// Changes: 7/15/05 None yet
//--------------------------------------------------------------------

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Torus.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4PVReplica.hh"
#include "G4RotationMatrix.hh"
#include "G4Trap.hh"
#include "G4Cons.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4NistManager.hh"

#include "GermSD.hh"
#include "NaIannSD.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"

#include "globals.hh"

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

DetectorConstruction::DetectorConstruction()
  : ActGeCrys_log(0),GeCrys_log(0)
{
  // Defaults for messenger variable values
  GeHLength = 93./2.*mm;
  GeRadius = 89.6/2.*mm;
  GeOffset = 6.*mm;
  GeEndRad = 6.9*mm;
  GeHoleRad = 8.7/2.*mm;
  GeHoleHLength = 79./2.*mm;
  GeFingerRad = 6.9/2.*mm;
  GeVertDisp = 0.*mm;
  GeHorisDisp = 0.*mm;
  GeDeadLayerHThick = 700./2.*um;
  Visuals = false;

  SourceDetectorDist = 10.525*mm;  // distance from detector face to target
  DetectorAngle = 0.*deg;

  ContactPinType = 3;   // Complicated geometry by default

  // The Detectormessenger
  detectorMessenger = new DetectorMessenger(this);

  DefineMaterials();
  SDflag=false;
}

DetectorConstruction::~DetectorConstruction()
{ }

void DetectorConstruction::DefineMaterials()
{
  //---------------Materials-----------------

  G4double a; //atomic mass
  G4double z; //atomic number
  G4double density, fractionmass;
  G4int ncomponents,natoms;

  Al = 
     new G4Material("Aluminum", z= 13., a= 26.98*g/mole, density= 2.7*g/cm3);

  Pb =
     new G4Material("Lead", z=82.,a=207.2*g/mole, density=11340.*kg/m3);

  Cu = 
     new G4Material("Copper", z=29.,a=63.546*g/mole, density=8920.*kg/m3);

  Ta =
     new G4Material("Tantalum",z=73.,a=180.9479*g/mole,density=16650*kg/m3);

  Ge =
     new G4Material("Germanium",z=32.,a=72.61*g/mole,density=5323*kg/m3);

  Li = 
     new G4Material("Lithium",z=3.,a=6.941*g/mole,density=0.534*g/cm3);

  B =
     new G4Material("Boron",z=5.,a=10.811*g/mole,density=2.31*g/cm3);

  // the elements
  a = 14.01*g/mole;
  G4Element* elN = new G4Element("Nitrogen","N",z=7.,a);
  a = 16.00*g/mole;
  G4Element* elO = new G4Element("Oxygen","O",z=8.,a);
  a = 12.011*g/mole;
  G4Element* elC = new G4Element("Carbon","C",z=6.,a);
  a = 1.00794*g/mole;
  G4Element* elH = new G4Element("Hydrogen","H",z=1.,a);
  a = 22.989768*g/mole;
  G4Element* elNa = new G4Element("Sodium","Na",z=11.,a);
  a = 126.90447*g/mole;
  G4Element* elI  = new G4Element("Iodine","I",z=53.,a);
  a = 204.3833*g/mole;
  G4Element* elTl = new G4Element("Thallium","Tl",z=81.,a);
  a = 58.9332*g/mole;
  G4Element* elCo = new G4Element("Cobalt", "Co", 27. ,a);
  a = 55.85*g/mole;
  G4Element* elFe = new G4Element("Iron" , "Fe" , z= 26., a);
  a = 63.546*g/mole;
  G4Element* elCu = new G4Element("Copper","Cu",z=29., a);
  a = 65.39*g/mole;
  G4Element* elZn = new G4Element("Zinc","Zn",z=30., a);
  a = 118.710*g/mole;
  G4Element* elSn = new G4Element("Tin","Sn",z=50.,a);
  a = 28.0855*g/mole;
  G4Element* elSi = new G4Element("Silicon","Si",z=14.,a);
  a = 26.98*g/mole;
  G4Element* elAl = new G4Element("Aluminium","Al",z=13.,a);


  // Air
  density = 1.290*mg/cm3;
  Air = new G4Material("Air",density,ncomponents=2);
  Air->AddElement(elN, fractionmass=0.7);
  Air->AddElement(elO, fractionmass=0.3);

  // Water
  Water = new G4Material("Water", density= 1.0*g/cm3,
  				     ncomponents=2);
  Water->AddElement(elH, 2);
  Water->AddElement(elO, 1);

  // Quartz
  Quartz = new G4Material("Quartz", density= 2.66*g/cm3,
  				     ncomponents=2);
  Quartz->AddElement(elSi, 1);
  Quartz->AddElement(elO, 2);

  // Aluminium Oxide ceramic
  density = 4.0 *g/cm3;
  ceramic = new G4Material("Aluminium Oxide Ceramic",density,ncomponents=2);
  ceramic -> AddElement(elAl,2./5.);
  ceramic -> AddElement(elO,3./5.);

  // NaI(Tl) material
  density=3.67*g/cm3;
  G4Material* NaI = new G4Material("NaI",density,ncomponents=2);
  NaI->AddElement(elNa,natoms=1);
  NaI->AddElement(elI,natoms=1);
  NaITl = new G4Material("NaI(Tl)",density,ncomponents=2);
  NaITl->AddMaterial(NaI,fractionmass=99.*perCent);
  NaITl->AddElement(elTl,fractionmass=1.*perCent);

  // Scintillator BC-408 Material
  Sci = 
    new G4Material("Scintillator", density= 1.032*g/cm3, ncomponents=2);
  Sci->AddElement(elC, natoms=10);
  Sci->AddElement(elH, natoms=11);

  // Stainless Steel
  ssteel = new G4Material
    ("Steel", density=7.7*g/cm3, ncomponents=3);
  ssteel->AddElement(elC, 0.04);
  ssteel->AddElement(elFe, 0.88);
  ssteel->AddElement(elCo, 0.08);

  // Mylar (alenia spazio)
  density= 1.4 *g/cm3;
  mylar = new G4Material("mylar",density,ncomponents=3);
  mylar -> AddElement(elH,0.042);
  mylar -> AddElement(elC,0.625);
  mylar -> AddElement(elO,0.333);

  // Naval Brass for cold finger
  density = 8.442 *g/cm3;
  brass = new G4Material("Naval Brass",density,ncomponents=3);
  brass -> AddElement(elCu,0.59);
  brass -> AddElement(elZn,0.40);
  brass -> AddElement(elSn,0.01);

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

  // NIST
  G4NistManager *man=G4NistManager::Instance();
  man->SetVerbose(1);

  G4cout << G4endl << *(G4Material::GetMaterialTable()) << G4endl;

  // Delrin
  delrin = man->FindOrBuildMaterial("G4_POLYOXYMETHYLENE");

  // Plastic
  plastic = man->FindOrBuildMaterial("G4_POLYSTYRENE");

//
// ------------ Generate & Add Material Properties Table for NaI(Tl) --------
//
  G4MaterialPropertiesTable* MPTNaITl = new G4MaterialPropertiesTable();
  const G4int nEntries = 2;
  G4double PhotonEnergy[nEntries] =
    { 2.26*eV, 3.89*eV };      
  G4double RefractiveIndex1[nEntries] =
    { 1.85, 1.85 };
  MPTNaITl->AddProperty("RINDEX",PhotonEnergy,RefractiveIndex1,nEntries);
  MPTNaITl->AddConstProperty("SCINTILLATIONYIELD",38./keV);
  NaITl->SetMaterialPropertiesTable(MPTNaITl);

  //Print the defined Materials
  G4cout << G4endl << "The Defined Materials are : " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}


G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // First clean out the old geom
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

//**************************************************************************

  //-------------------------------------------------
  //------------------Volumes------------------------
  //-------------------------------------------------
  // z - Along Beam
  // x - Horizontal plane
  // y - vertical plane
  //
  // ----------WORLD-------------

  LENALabz = 1500*mm;
  LENALabx = 1500*mm;
  LENALaby = 1500*mm;

  G4double HalfLENALabz = LENALabz*0.5;
  G4double HalfLENALabx = LENALabx*0.5;
  G4double HalfLENALaby = LENALaby*0.5;

  G4Box* LENALab_solid
    = new G4Box("LENALab",HalfLENALabx,HalfLENALaby,HalfLENALabz);
  LENALab_log = 
    new G4LogicalVolume(LENALab_solid, Air, "LENALab_log", 0,0,0);
  LENALab_phys = 
    new G4PVPlacement(0,                // no rotation
		      G4ThreeVector(),  // at (0,0,0)
		      LENALab_log,      // Logical Volume
		      "LENALab",        // Name
		      0,                // Mother Volume
		      false,            // no boolean operations
		      0);               // Copy number

  // The World is invisible
  LENALab_log ->SetVisAttributes(G4VisAttributes::Invisible);



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  // Build setup

//  BuildHPGe();
//  BuildNaIAnnulus();
//  BuildGammaCube();
  BuildAPEX();
  BuildAPEXcradle();
  BuildTargetHolder();
//  BuildCollimator();
//  BuildRadioactiveSource();

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



  // Sensitive Ge Detector
  if(!SDflag){
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    G4String SDname;
    
    GeSD = new GermSD(SDname="/GeSD", this);
    SDman->AddNewDetector(GeSD);
  
    NaISD = new NaIannSD(SDname="/NaISD", this);
    SDman->AddNewDetector(NaISD);

    SDflag=true;
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// HPGe
//  ActGeCrys_log->SetSensitiveDetector(GeSD);

// Annulus
//  NaISeg_log->SetSensitiveDetector(NaISD);

// APEX
  scintillator_log->SetSensitiveDetector(NaISD);
//  scintillator_log->SetSensitiveDetector(GeSD);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  //-----------------------------------------------------------------
  return LENALab_phys;
}


  //------------------------------------------------------------
  //-------------- Sensitive Detectors -----------------
  //-----------------------------------------------------------


// The commands to change detector geommetry
void DetectorConstruction::SetGeHLength(G4double val)
{
  GeHLength = val;
}
void DetectorConstruction::SetGeRadius(G4double val)
{
  GeRadius = val;
}

void DetectorConstruction::SetGeOffset(G4double val)
{
  GeOffset = val;
}
void DetectorConstruction::SetGeEndRad(G4double val)
{
  GeEndRad = val;
}
void DetectorConstruction::SetGeHoleRad(G4double val)
{
  GeHoleRad = val;
}
void DetectorConstruction::SetGeHoleHLength(G4double val)
{
  GeHoleHLength = val;
}
void DetectorConstruction::SetGeFingerRad(G4double val)
{
  GeFingerRad = val;
}
void DetectorConstruction::SetGeVertDisp(G4double val)
{
  GeVertDisp = val;
}
void DetectorConstruction::SetGeHorisDisp(G4double val)
{
  GeHorisDisp = val;
}
void DetectorConstruction::SetGeDeadLayerHThick(G4double val)
{
  GeDeadLayerHThick = val;
}
void DetectorConstruction::SetVisuals(G4bool val)
{
  Visuals = val;
}
void DetectorConstruction::SetSourceDetectorDist(G4double val)
{
  SourceDetectorDist = val;
}
void DetectorConstruction::SetDetectorAngle(G4double val)
{
  DetectorAngle = val;
}
void DetectorConstruction::SetContactPinType(G4int val)
{
  ContactPinType = val;
}


#include "G4RunManager.hh"
void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}

//#####################################################################
// The Germanium Detector
void DetectorConstruction::BuildHPGe(){

  //Source_Detectorz = 10.525*mm; // from the drawings for the target holder
  Source_Detectorz = 0*mm;
  // new source detector distances
  // Source_Detectorz = 5.525*mm;

  // The Ge Detector Mother Volume
  
  // The crystal housing is the largest thing - so pre-define it here
  G4double houseFrontThickness = 1.5/2.; // the HALF thickness
  G4double houseSideThickness = 1.3/2.*mm;
  G4double houseInnerRadius = 105/2.*mm;
  G4double houseOuterRadius = houseInnerRadius+2*houseSideThickness;
  G4double houseHLength = 100.*mm;
  G4double houseBaseRadius = houseInnerRadius;

  G4double GeMotherRadius = houseOuterRadius;
  G4double GeMotherHLength = houseHLength;
  //G4double GetoTargetLength = 175.0*mm;
  G4double GeMotherZpos = (SourceDetectorDist+GeMotherHLength)
    *cos(DetectorAngle);
  G4double GeMotherXpos = -(SourceDetectorDist+GeMotherHLength)
    *sin(DetectorAngle);
  G4Tubs* GeMotherCylinder_solid
    = new G4Tubs("Ge_Mother",0,
  		 GeMotherRadius,GeMotherHLength,0,twopi);
  G4LogicalVolume* GeMotherCylinder_log =
    new G4LogicalVolume(GeMotherCylinder_solid, vacuum, "GeMother log",
			0,0,0);

  G4RotationMatrix* GeMother_rot = new G4RotationMatrix();
  GeMother_rot->rotateY(DetectorAngle);

  new G4PVPlacement(GeMother_rot,
		    G4ThreeVector(GeMotherXpos,
				  0.*mm,
				  GeMotherZpos),
		      GeMotherCylinder_log,
		      "Germanium Detector Mother Volume",
		      LENALab_log,
		      false,
		      0);

  //---------------------------------------------
  // The HPGe (including dead layer)
  //G4cout << "radius is " << GeRadius << G4endl;
  //GeEndRad*=3.;
  //GeHoleRad*=2.;
  G4double GeCrys_outer = GeRadius;
  G4double GeCrys_Hlen  = GeHLength;
  G4double GeBlock_Hlen = GeCrys_Hlen-GeEndRad/2.;
  G4double GeFinger_outer = GeHoleRad;
  G4double GeFinger_Hlen  = GeHoleHLength-GeHoleRad/2.;
  G4double GeCrys_posz = GeOffset + GeCrys_Hlen + Source_Detectorz;
  G4double GeFinger_posz = -GeFinger_Hlen+GeBlock_Hlen; // relative to crystal
  G4RotationMatrix noRotate;
  G4ThreeVector MoveFinger(0,0,GeFinger_posz);

  G4Tubs* blockGe
    = new G4Tubs("Block of Ge (no finger)", 0.*mm,GeCrys_outer,
		 GeBlock_Hlen,0,twopi);

  // Curved end of crystal
  G4double endGe_torusRad = GeCrys_outer-GeEndRad;
  G4double endGe_posz = -GeBlock_Hlen;
  G4ThreeVector MoveendGe(0,0,endGe_posz);
  G4ThreeVector MoveendGeVis(GeHorisDisp,
			     GeVertDisp,
			     endGe_posz+GeEndRad/2.+GeCrys_posz-GeMotherHLength);
  G4Torus* endGe_solid
    = new G4Torus("Curved end of Ge",   // name
		  0.,                   // Inner Radius
		  GeEndRad,             // outer radius
		  endGe_torusRad,       // The torus radius
		  0,twopi);             // start and end angle

  // This is purely for visualisation - take out for real simulations
  G4LogicalVolume* endGeVIS_log =
      new G4LogicalVolume(endGe_solid,Ge,"End Curve",0,0,0);
  Visuals=0;
  if(Visuals){
    //    endGeVIS_log =
    //      new G4LogicalVolume(endGe_solid,Ge,"End Curve",0,0,0);
    // G4PhysicalVolume* endGe =
    new G4PVPlacement(0,
		      MoveendGeVis,
		      endGeVIS_log,
		      "End Curve",
		      GeMotherCylinder_log,
		      false,
		      0);
  }
  G4UnionSolid* blockge2
    = new G4UnionSolid("blockge2",blockGe,endGe_solid,&noRotate,MoveendGe);

  // Now add in the extra material in the end
  G4double FillPiece_HLenX = 2 *mm;//add extra length to overlap with other HPGe (if they just touch, Geant4 will a seam (or skin)
  G4double FillPiece_shiftX = 2 *mm;
  G4double FillPiece_HLen = GeEndRad/2.;
  G4double FillPiece_Rad = endGe_torusRad;
  G4double FillPiece_posz = endGe_posz-GeEndRad/2.;
  G4ThreeVector MoveFillPiece(0,0,FillPiece_posz+FillPiece_shiftX);
  G4Tubs* FillPiece_solid
    = new G4Tubs("Fill Piece for end of HPGe",
		 0.*mm,                  // inner radius
		 FillPiece_Rad,          // outer radius
		 FillPiece_HLen+FillPiece_HLenX,         // half length
		 0,twopi);               // start and end angle
  G4UnionSolid* blockge3
    = new G4UnionSolid("HPGe Crystal",blockge2,FillPiece_solid,
    		       &noRotate,MoveFillPiece);
  //= new G4UnionSolid("HPGe Crystal",blockGe,FillPiece_solid,
  //	       &noRotate,MoveFillPiece);


  G4Tubs* fingerGe
    = new G4Tubs("Ge Cold Finger", 0.*mm,GeFinger_outer,GeFinger_Hlen,
		 0,twopi);
  G4SubtractionSolid* blockge4
    = new G4SubtractionSolid("Ge Crystal 3",blockge3,fingerGe,
			     &noRotate,MoveFinger);

  // Curved end in Cold finger
  G4double gefingercurvePosz = GeFinger_posz - GeFinger_Hlen;
  G4ThreeVector Movefingercurve(0,0,gefingercurvePosz);
  G4Sphere* gefingercurve =
    new G4Sphere("Curve in cold finger", 0.*mm,GeFinger_outer-0.05*mm,0,twopi,0,pi);
  //GeFinger_outer,0,twopi,0,pi);

  G4SubtractionSolid* GeCrys_solid
    = new G4SubtractionSolid("Ge Crystal",blockge4,gefingercurve,&noRotate,
			     Movefingercurve);

  // The final HPGe crystal logical volume
  GeCrys_log = 
    new G4LogicalVolume(GeCrys_solid,Ge,"Ge Crystal log",0,0,0);
  //new G4LogicalVolume(GeCrys_solid,Ge,"Ge Crystal log",0,0,0);
  // G4VPhysicalVolume* GeCrys_phys = 
  new G4PVPlacement(0, 
		    G4ThreeVector(GeHorisDisp,
				  GeVertDisp,
				  GeCrys_posz+FillPiece_HLen-GeMotherHLength),
		    GeCrys_log,
		    "Ge Crystal",
		    GeMotherCylinder_log,
		    false,
		    0);

  //--------------------------------------------------------------
  // Now, just the active volume of crystal
  G4double GeActCrys_outer = GeRadius-GeDeadLayerHThick*2.;
  G4double GeActCrys_Hlen  = GeHLength-GeDeadLayerHThick;
  G4double GeActEndRad     = GeEndRad; //-GeDeadLayerHThick*2.;
  G4double GeActBlock_Hlen = GeActCrys_Hlen-GeActEndRad/2.;
  G4double GeActFinger_outer = GeHoleRad;
  G4double GeActFinger_Hlen  = GeHoleHLength-GeHoleRad/2.;
  // Positions are with respect to whole HPGe material
  //G4double GeActCrys_posz = 0.;
  G4double GeActFinger_posz = -GeActFinger_Hlen+GeActBlock_Hlen;
  G4ThreeVector MoveActFinger(0,0,GeActFinger_posz);

  G4Tubs* ActblockGe
    = new G4Tubs("Active Block of Ge (no finger)", 0.*mm,GeActCrys_outer,
		 GeActBlock_Hlen,0,twopi);

  // Curved end of crystal
  G4double endGeAct_torusRad = GeActCrys_outer-GeActEndRad;
  G4double endGeAct_posz = -GeActBlock_Hlen;
  G4ThreeVector MoveActendGe(0,0,endGeAct_posz);
  G4ThreeVector MoveActendGeVis(0,
				0,
				endGeAct_posz); //-GeActEndRad+GeActCrys_posz);
  G4Torus* ActendGe_solid
    = new G4Torus("Active Curved end of Ge",   // name
		  0.,                   // Inner Radius
		  GeActEndRad,             // outer radius
		  endGeAct_torusRad,       // The torus radius
		  0,twopi);             // start and end angle




  G4UnionSolid* Actblockge2
    = new G4UnionSolid("Actblockge2",ActblockGe,ActendGe_solid,&noRotate,MoveActendGe);

  // Now add in the extra material in the end
  G4double ActFillPiece_HLen = GeEndRad/2.;
  G4double ActFillPiece_Rad = endGeAct_torusRad;
  G4double ActFillPiece_posz = endGeAct_posz-GeEndRad/2.;
  G4ThreeVector MoveActFillPiece(0,0,ActFillPiece_posz+FillPiece_shiftX);

  G4Tubs* ActFillPiece_solid
    = new G4Tubs("Active Fill Piece for end of HPGe",
		 0.*mm,                  // inner radius
		 ActFillPiece_Rad,          // outer radius
		 ActFillPiece_HLen+FillPiece_HLenX,         // half length
		 0,twopi);               // start and end angle

  G4UnionSolid* Actblockge3
    = new G4UnionSolid("Active HPGe Crystal",Actblockge2,ActFillPiece_solid,
		       &noRotate,MoveActFillPiece);

  G4Tubs* ActfingerGe
    = new G4Tubs("Active Ge Cold Finger", 0.*mm,GeActFinger_outer,GeActFinger_Hlen,
		 0,twopi);
  G4SubtractionSolid* Actblockge4
    = new G4SubtractionSolid("Active Ge Crystal 3",Actblockge3,ActfingerGe,
			     &noRotate,MoveActFinger);

  // Curved end in Cold finger
  G4double ActgefingercurvePosz = GeActFinger_posz - GeActFinger_Hlen;
  G4ThreeVector MoveActfingercurve(0,0,ActgefingercurvePosz);
  G4Sphere* Actgefingercurve =
    new G4Sphere("Active Curve in cold finger", 0.*mm,GeActFinger_outer-0.05*mm,
		 0,twopi,0,pi);

  G4SubtractionSolid* ActGeCrys_solid
    = new G4SubtractionSolid("Active Ge Crystal",Actblockge4,Actgefingercurve,&noRotate,
			     MoveActfingercurve);

  // The final HPGe crystal logical volume
  ActGeCrys_log = 
    new G4LogicalVolume(ActGeCrys_solid,Ge,"Active Ge Crystal log",0,0,0);
  //new G4LogicalVolume(GeCrys_solid,Ge,"Ge Crystal log",0,0,0);
  // G4VPhysicalVolume* GeCrys_phys = 
  new G4PVPlacement(0, 
		    G4ThreeVector(0.,
				  0.,
				  GeDeadLayerHThick),
		    ActGeCrys_log,
		    "Active Ge Crystal",
		    GeCrys_log,
		    false,
		    0);

  
  // This is purely for visualisation - take out for real simulations
  G4LogicalVolume* ActendGeVIS_log =
      new G4LogicalVolume(ActendGe_solid,Ge,"Active End Curve",0,0,0);
  Visuals=0;
  if(Visuals){
    //    endGeVIS_log =
    //      new G4LogicalVolume(endGe_solid,Ge,"End Curve",0,0,0);
    // G4PhysicalVolume* endGe =
    new G4PVPlacement(0,
		      MoveActendGeVis,
		      ActendGeVIS_log,
		      "Active End Curve",
		      ActGeCrys_log,
		      false,
		      0);
  }
  // ---------------------------------------------------------


  // Detector Housing can
  //// normal way of making housing
  G4Tubs* HouseTube_solid
    = new G4Tubs("Housing Tube",houseInnerRadius,
  		 houseOuterRadius,houseHLength,0,twopi);
  G4Tubs* HouseBase_solid
    = new G4Tubs("Housing Base",0.*mm,houseBaseRadius,
  		 houseFrontThickness,0,twopi);
  G4LogicalVolume* HouseTube_log = 
    new G4LogicalVolume(HouseTube_solid, Al, "HouseTube log",
			0,0,0);
  G4LogicalVolume* HouseBase_log = 
    new G4LogicalVolume(HouseBase_solid, Al, "HouseBase log",
			0,0,0);
  // G4VPhysicalVolume* HouseTube_phys = 
    new G4PVPlacement(0,
		      G4ThreeVector(0.*mm,
				    0.*mm,
				    Source_Detectorz+houseHLength-GeMotherHLength),
		      HouseTube_log,
		      "Ge Housing 1",
		      GeMotherCylinder_log,
		      false,
		      0);
  // G4VPhysicalVolume* HouseBase_phys = 
    new G4PVPlacement(0,
		      G4ThreeVector(0.*mm,
				    0.*mm,
				    Source_Detectorz+houseFrontThickness
				    -GeMotherHLength),
		      HouseBase_log,
		      "Ge Housing 2",
		      GeMotherCylinder_log,
		      false,
		      0);



  // G4VPhysicalVolume* GePed_phys = 
//     new G4PVPlacement(0,
// 		      G4ThreeVector(0.*mm,
// 				    0.*mm,
// 				    pedPosz-GeMotherHLength),
// 		      GePed_log,
// 		      "Germanium Pedestal",
// 		      GeMotherCylinder_log,
// 		      false,
// 		      0);

  // Lithium Contact around crystal
//   G4double lithHThick = 0.7/2.*mm; // the HALF thickness
//   G4double lithInnerRadius = GeCrys_outer;
//   G4double lithOuterRadius = lithInnerRadius+2*lithHThick;
//   G4double lithHLength = GeCrys_Hlen+2*lithHThick+0.5*mm;
//   G4double lithBaseRadius = lithInnerRadius;
//   G4double lithPosz = GeCrys_posz;

//     G4Tubs* lithTube_solid
//     = new G4Tubs("Lith",lithInnerRadius,
//   		 lithOuterRadius,lithHLength,0,twopi);
//   G4Tubs* lithBase_solid
//     = new G4Tubs("lith Base",0.*mm,lithBaseRadius,
//   		 lithHThick,0,twopi);
//   G4LogicalVolume* lithTube_log = 
//     new G4LogicalVolume(lithTube_solid, Li, "lithTube log",0,0,0);
//   G4LogicalVolume* lithBase_log = 
//     new G4LogicalVolume(lithBase_solid, Li, "lithBase log",0,0,0);
//   // G4VPhysicalVolume* lithTube_phys =
//     new G4PVPlacement(0,
// 		      G4ThreeVector(GeHorisDisp,
// 				    GeVertDisp,
// 				    lithPosz-GeMotherHLength),
// 		      lithTube_log,
// 		      "Ge Lith 1",
// 		      GeMotherCylinder_log,
// 		      false,
// 		      0);

//   // G4VPhysicalVolume* lithBase_phys = 
//     new G4PVPlacement(0,
// 		      G4ThreeVector(GeHorisDisp,
// 				    GeVertDisp,
// 				    lithPosz-lithHLength+lithHThick-GeMotherHLength),
// 		      lithBase_log,
// 		      "Ge Lith 2",
// 		      GeMotherCylinder_log,
// 		      false,
// 		      0);

  // Mylar surrounding crystal
  G4double mylarHThick = 0.03/2.*mm; // the HALF thickness
  G4double mylarInnerRadius = GeCrys_outer;
  G4double mylarOuterRadius = mylarInnerRadius+2*mylarHThick;
  G4double mylarHLength = GeCrys_Hlen+2*mylarHThick+0.05;
  G4double mylarBaseRadius = mylarInnerRadius;
  G4double mylarPosz = GeCrys_posz;

    G4Tubs* mylarTube_solid
    = new G4Tubs("Mylar",mylarInnerRadius,
  		 mylarOuterRadius,mylarHLength,0,twopi);
  G4Tubs* mylarBase_solid
    = new G4Tubs("mylar Base",0.*mm,mylarBaseRadius,
  		 mylarHThick,0,twopi);
  G4LogicalVolume* mylarTube_log = 
    new G4LogicalVolume(mylarTube_solid, mylar, "mylarTube log",0,0,0);
  G4LogicalVolume* mylarBase_log = 
    new G4LogicalVolume(mylarBase_solid, mylar, "mylarBase log",0,0,0);
  // G4VPhysicalVolume* mylarTube_phys = 
  new G4PVPlacement(0,
		    G4ThreeVector(GeHorisDisp,
				  GeVertDisp,
				  mylarPosz-GeMotherHLength),
		    mylarTube_log,
		    "Ge Mylar 1",
		    GeMotherCylinder_log,
		    false,
		    0);
  // G4VPhysicalVolume* mylarBase_phys = 
  new G4PVPlacement(0,
		    G4ThreeVector(GeHorisDisp,
				  GeVertDisp,
				  mylarPosz-mylarHLength+mylarHThick-GeMotherHLength),
		    mylarBase_log,
		    "Ge Mylar Base",
		    GeMotherCylinder_log,
		    false,
		      0);
  
  // Aluminized mylar
  G4double alummylarHThick = mylarHThick;
  G4double alummylarInnerRadius = mylarOuterRadius;
  G4double alummylarOuterRadius = alummylarInnerRadius+2*alummylarHThick;
  G4double alummylarHLength = mylarHLength;
  G4double alummylarBaseRadius = alummylarInnerRadius;
  G4double alummylarPosz = GeCrys_posz;

  G4Tubs* alummylarTube_solid
    = new G4Tubs("Aluminised mylar",alummylarInnerRadius,
  		 alummylarOuterRadius,alummylarHLength,0,twopi);
  G4Tubs* alummylarBase_solid
    = new G4Tubs("Aluminised mylar Base",0.*mm,alummylarBaseRadius,
  		 alummylarHThick,0,twopi);

  G4LogicalVolume* alummylarTube_log = 
    new G4LogicalVolume(alummylarTube_solid, mylar,"alummylarTube log",0,0,0);
  G4LogicalVolume* alummylarBase_log = 
    new G4LogicalVolume(alummylarBase_solid, mylar, "alummylarBase log",0,0,0);

  // G4VPhysicalVolume*almylarTube_phys = 
  new G4PVPlacement(0,
		    G4ThreeVector(GeHorisDisp,
				  GeVertDisp,
				  alummylarPosz-GeMotherHLength),
		    alummylarTube_log,
		    "Ge Aluminised mylar",
		    GeMotherCylinder_log,
		    false,
		    0);
  // G4VPhysicalVolume* alummylarBase_phys = 
  new G4PVPlacement(0,
		    G4ThreeVector(GeHorisDisp,
				  GeVertDisp,
				  alummylarPosz-alummylarHLength-2.*mylarHThick+
				  alummylarHThick-GeMotherHLength),
		    alummylarBase_log,
		    "Ge aluminised Mylar Base",
		    GeMotherCylinder_log,
		    false,
		    0);


  // Aluminium surrounding mylar
  G4double almylarHThick = 0.76/2.*mm; // the HALF thickness
  G4double almylarInnerRadius = alummylarOuterRadius;
  G4double almylarOuterRadius = almylarInnerRadius+2*almylarHThick;
  G4double almylarHLength = 130./2.*mm; //mylarHLength+2.*almylarHThick;
  //G4double almylarBaseRadius = almylarInnerRadius;
  G4double almylarPosz = GeCrys_posz+almylarHLength-GeHLength;

    G4Tubs* almylarTube_solid
    = new G4Tubs("Almylar",almylarInnerRadius,
  		 almylarOuterRadius,almylarHLength,0,twopi);
//   G4Tubs* almylarBase_solid
//     = new G4Tubs("almylar Base",0.*mm,almylarBaseRadius,
//   		 almylarHThick,0,twopi);
  G4LogicalVolume* almylarTube_log = 
    new G4LogicalVolume(almylarTube_solid, Al,"almylarTube log",0,0,0);
//   G4LogicalVolume* almylarBase_log = 
//     new G4LogicalVolume(almylarBase_solid, Al,"almylarBase log",0,0,0);
  // G4VPhysicalVolume*almylarTube_phys = 
    new G4PVPlacement(0,
		      G4ThreeVector(GeHorisDisp,
				    GeVertDisp,
				    almylarPosz-GeMotherHLength),
		      almylarTube_log,
		      "Ge Almylar 1",
		      GeMotherCylinder_log,
		      false,
		      0);

    // Aluminium Base 
    G4double baseRad = almylarInnerRadius;   //89.6/2.*mm;
  G4double baseHThick = 3.2/2.*mm;
  G4double basePosz = GeCrys_posz - GeCrys_Hlen + 130. - baseHThick;
  //G4double pedRad = 44.5/2.*mm;
  //G4double pedHThick = baseHThick;
  //G4double pedPosz = basePosz + 0.4+2*pedHThick;

  G4Tubs* GeBase_solid
    = new G4Tubs("Germanium base platform", 0.*mm,baseRad,baseHThick,0,twopi);
  //G4Tubs* GePed_solid
  //  = new G4Tubs("Germanium Pedestal", 0.*mm,pedRad,pedHThick,0,twopi);
  G4LogicalVolume* GeBase_log = 
    new G4LogicalVolume(GeBase_solid,Al,"Germanium Base log",0,0,0);
  //G4LogicalVolume* GePed_log = 
  //  new G4LogicalVolume(GePed_solid,Al,"Germanium Pedestal log",0,0,0);
  // G4VPhysicalVolume* GeBase_phys = 
    new G4PVPlacement(0,
		      G4ThreeVector(GeHorisDisp,
				    GeVertDisp,
				    basePosz-GeMotherHLength),
		      GeBase_log,
		      "Germanium Base Platform",
		      GeMotherCylinder_log,
		      false,
		      0);
  // G4VPhysicalVolume* almylarBase_phys = 
//     new G4PVPlacement(0,
// 		      G4ThreeVector(GeHorisDisp,
// 				    GeVertDisp,
// 				    almylarPosz-almylarHLength-almylarHThick-GeMotherHLength),
// 		      almylarBase_log,
// 		      "Ge Almylar Base",
// 		      GeMotherCylinder_log,
// 		      false,
// 		      0);

  // Boron coating inside hole
    G4double boron_outer = GeFinger_outer;
  G4double boron_inner = boron_outer-0.0003/2.*mm;
  G4double boron_Hlen  = GeFinger_Hlen-4.35/2.;
  G4double boron_posz = GeCrys_posz+FillPiece_HLen+GeFinger_posz; 
  // relative to crystal

  G4Tubs* boron_solid
    = new G4Tubs("boron",boron_inner ,boron_outer,
		 boron_Hlen,0,twopi);
  G4LogicalVolume* boron_log = 
    new G4LogicalVolume(boron_solid,B,"Boron log",0,0,0);
  // G4VPhysicalVolume* boron_phys =
    new G4PVPlacement(0,
		      G4ThreeVector(GeHorisDisp,
				    GeVertDisp,
				    boron_posz-GeMotherHLength),
		      boron_log,
		      "Boron in cold finger",
		      GeMotherCylinder_log,
		      false,
		      0);


    // The copper cold finger
    G4double fingerNutRad = 10.8*mm;
    G4double fingerNutHLength = 2.*mm;
    G4double fingerOuter = GeFingerRad;   //6.9/2.*mm;
    G4double fingerOffset = 2.175*mm;
    G4double fingerHLength = almylarHLength-(GeCrys_Hlen-GeHoleHLength)-
      fingerNutHLength-fingerOffset/2.-baseHThick;
    G4double fingerPosz = GeCrys_posz+(GeCrys_Hlen+fingerHLength-
				       2.*GeHoleHLength)+fingerOffset;
    G4double fingerNutPosz = fingerPosz+fingerHLength+fingerNutHLength;
    
    G4LogicalVolume* __attribute__ ((unused)) CuFinger_log=0;
    G4LogicalVolume* __attribute__ ((unused)) brassFingerNarrow_log=0;
    G4LogicalVolume* __attribute__ ((unused)) brassFingerWide_log=0;

    if(ContactPinType == 2){
      fingerOffset = 2.175*mm;
      fingerHLength = almylarHLength-(GeCrys_Hlen-GeHoleHLength)-
	fingerNutHLength-fingerOffset/2.-baseHThick;
      fingerPosz = GeCrys_posz+(GeCrys_Hlen+fingerHLength-
				2.*GeHoleHLength)+fingerOffset;
      fingerNutPosz = fingerPosz+fingerHLength+fingerNutHLength;
      G4Tubs* CuFinger_solid 
	= new G4Tubs("Copper Cold Finger",0.*mm,fingerOuter,fingerHLength,0,twopi);
      G4LogicalVolume* CuFinger_log = 
	new G4LogicalVolume(CuFinger_solid,Cu,"CuFinger log",0,0,0);
      // G4VPhysicalVolume* CuFinger_phys =
      new G4PVPlacement(0,
			G4ThreeVector(GeHorisDisp,
				      GeVertDisp,
				      fingerPosz-GeMotherHLength),
			CuFinger_log,
			"Cu Cold Finger",
			GeMotherCylinder_log,
			false,
			0);

      CuFinger_log->SetVisAttributes(G4VisAttributes(G4Colour(0.84,0.49,0.1)));
    }

    if(ContactPinType == 3){
      fingerOffset = 7.5*mm;
      fingerHLength = almylarHLength-(GeCrys_Hlen-GeHoleHLength)-
	fingerNutHLength-fingerOffset/2.-baseHThick;
      fingerPosz = GeCrys_posz+(GeCrys_Hlen+fingerHLength-
				2.*GeHoleHLength)+fingerOffset;
      fingerNutPosz = fingerPosz+fingerHLength+fingerNutHLength;

      // The holes inside the copper
      G4double fingerHoleNarrowRad = 3.0/2.*mm;
      G4double fingerHoleNarrowHLength = 31./2.*mm + 0.5*mm;
      G4double fingerHoleWideRad = 4.75/2.*mm;
      G4double fingerHoleWideHLength = 44.7/2.*mm + 0.1*mm;
      G4double fingerHoleNarrowPosz = -fingerHLength +             // w.r.t finger
	fingerHoleNarrowHLength + 2.*fingerHoleWideHLength - 0.5*mm;
      G4double fingerHoleWidePosz = -fingerHLength +             // w.r.t finger
	fingerHoleWideHLength - 0.1*mm;
      G4ThreeVector MovefingerHoleNarrow(0,0,fingerHoleNarrowPosz);
      G4ThreeVector MovefingerHoleWide(0,0,fingerHoleWidePosz);

      G4Tubs* CuFinger_solid 
	= new G4Tubs("Copper Cold Finger",0.*mm,fingerOuter,fingerHLength,0,twopi);
      G4Tubs* CuFingerHoleNarrow
	= new G4Tubs("Cold Finger hole 1",0.*mm,fingerHoleNarrowRad,
		     fingerHoleNarrowHLength,0,twopi);
      G4Tubs* CuFingerHoleWide
	= new G4Tubs("Cold Finger hole 2",0.*mm,fingerHoleWideRad,
		     fingerHoleWideHLength,0,twopi);

      // subtract wide hole
      G4SubtractionSolid* CuFinger2_solid
	= new G4SubtractionSolid("Cu Finger 2",CuFinger_solid,CuFingerHoleWide,&noRotate,
				 MovefingerHoleWide);
      // Now subtract the narrow hole
      G4SubtractionSolid* CuFinger3_solid
	= new G4SubtractionSolid("Cu Finger 3",CuFinger2_solid,CuFingerHoleNarrow,
				 &noRotate,MovefingerHoleNarrow);

      G4LogicalVolume* CuFinger_log = 
	new G4LogicalVolume(CuFinger3_solid,Cu,"CuFinger log",0,0,0);
      // G4VPhysicalVolume* CuFinger_phys =
      new G4PVPlacement(0,
			G4ThreeVector(GeHorisDisp,
				      GeVertDisp,
				      fingerPosz-GeMotherHLength),
			CuFinger_log,
			"Cu Cold Finger",
			GeMotherCylinder_log,
			false,
			0);
      CuFinger_log->SetVisAttributes(G4VisAttributes(G4Colour(0.84,0.49,0.1)));
    

      // Now, we need to ass the naval brass finger inside the copper finger
      G4double brassFingerNarrowRad = 2.9/2.*mm;
      G4double brassFingerNarrowHLength = 51.6/2.*mm;
      G4double brassFingerWideRad = 3.8/2.*mm;
      G4double brassFingerWideHLength = 25.2/2.*mm;
      G4double brassFingerExtension = fingerOffset-0.5*mm;   // trial and error
      G4double brassFingerNarrowPosz = fingerPosz-fingerHLength+
	brassFingerNarrowHLength+2.*brassFingerWideHLength-brassFingerExtension;
      G4double brassFingerWidePosz = fingerPosz-fingerHLength+
	brassFingerWideHLength-brassFingerExtension;

      G4Tubs* brassFingerNarrow_solid
	= new G4Tubs("Narrow Brass Finger",0.*mm,brassFingerNarrowRad,
		     brassFingerNarrowHLength,0,twopi);
      G4LogicalVolume* brassFingerNarrow_log = 
	new G4LogicalVolume(brassFingerNarrow_solid,brass,"Narrow brass Finger",0,0,0);
      //    G4VPhysicalVolume* brassFingerNarrow_phys =
      new G4PVPlacement(0,
			G4ThreeVector(GeHorisDisp,
				      GeVertDisp,
				      brassFingerNarrowPosz-GeMotherHLength),
			brassFingerNarrow_log,
			"Narrow Brass Finger",
			GeMotherCylinder_log,
			false,
			0);

      G4Tubs* brassFingerWide_solid
	= new G4Tubs("Wide Brass Finger",0.*mm,brassFingerWideRad,
		     brassFingerWideHLength,0,twopi);
      G4LogicalVolume* brassFingerWide_log = 
	new G4LogicalVolume(brassFingerWide_solid,brass,"Wide brass Finger",0,0,0);
      //    G4VPhysicalVolume* brassFingerWide_phys =
      new G4PVPlacement(0,
			G4ThreeVector(GeHorisDisp,
				      GeVertDisp,
				      brassFingerWidePosz-GeMotherHLength),
			brassFingerWide_log,
			"Wide Brass Finger",
			GeMotherCylinder_log,
			false,
			0);

      G4VisAttributes* BrassVisAtt = new G4VisAttributes(G4Colour(0.74,0.39,0.1));
      brassFingerNarrow_log ->SetVisAttributes(BrassVisAtt);
      brassFingerWide_log ->SetVisAttributes(BrassVisAtt);

    }


    // Copper nut
    G4Tubs * CuNut_solid
      = new G4Tubs("Copper Nut",0.*mm,fingerNutRad,fingerNutHLength,0,twopi);
    G4LogicalVolume* CuNut_log =
      new G4LogicalVolume(CuNut_solid,Cu,"CuNut log",0,0,0);
    new G4PVPlacement(0,
		      G4ThreeVector(GeHorisDisp,
				    GeVertDisp,
				    fingerNutPosz-GeMotherHLength),
		      CuNut_log,
		      "Cu Nut",
		      GeMotherCylinder_log,
		      false,
		      0);


  //------------------------------------------------------------
  //-------------- Visualization Stuff -----------------------
  //------------------------------------------------------------

  // The World is invisible
  GeMotherCylinder_log->SetVisAttributes(G4VisAttributes::Invisible);

  // Make any aluminium off White/grey
  //  G4VisAttributes* AlVisAtt = new G4VisAttributes(G4Colour(0.8,0.8,0.8));
  G4VisAttributes* AlVisAtt2 = new G4VisAttributes(G4Colour(0.8,0.8,0.8));
						   //G4VisAttributes::Invisible);
  //AlVisAtt2->SetForceWireframe(true);
  GeBase_log->SetVisAttributes(AlVisAtt2);
  HouseBase_log->SetVisAttributes(AlVisAtt2);
  HouseTube_log->SetVisAttributes(AlVisAtt2);
  almylarTube_log->SetVisAttributes(AlVisAtt2);
  alummylarTube_log->SetVisAttributes(G4VisAttributes::Invisible);
  alummylarBase_log->SetVisAttributes(G4VisAttributes::Invisible);
  
  // Mylar parts are white
  //G4VisAttributes* mylarVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  G4VisAttributes* mylarVisAtt = 
    new G4VisAttributes(G4VisAttributes::Invisible);
  mylarTube_log->SetVisAttributes(mylarVisAtt);
  mylarBase_log->SetVisAttributes(mylarVisAtt);

  // Copper is Copper Colour
  G4VisAttributes* CuVisAtt = new G4VisAttributes(G4Colour(0.84,0.49,0.1));
//   if(ContactPinType == 2 || ContactPinType == 3){ 
//     CuFinger_log->SetVisAttributes(CuVisAtt);
//   }
  CuNut_log->SetVisAttributes(CuVisAtt);

//   if(ContactPinType == 3){
//     // Brass is dark copper
//     G4VisAttributes* BrassVisAtt = new G4VisAttributes(G4Colour(0.74,0.39,0.1));
//     brassFingerNarrow_log ->SetVisAttributes(BrassVisAtt);
//     brassFingerWide_log ->SetVisAttributes(BrassVisAtt);
//   }

  // The Detectors are green
  G4VisAttributes* GeVisAtt = new G4VisAttributes(G4Colour(0.0,0.5,0.0));
  //GeVisAtt->SetForceWireframe(true);
  GeCrys_log->SetVisAttributes(GeVisAtt);
  G4VisAttributes* ActGeVisAtt = new G4VisAttributes(G4Colour(1.,1.,1.));
  //ActGeVisAtt->SetForceWireframe(true);
  ActGeCrys_log->SetVisAttributes(ActGeVisAtt);
  Visuals=0;
  if(Visuals){
    endGeVIS_log->SetVisAttributes(GeVisAtt);
    ActendGeVIS_log->SetVisAttributes(G4VisAttributes::Invisible);
  }

  // Make boron contact invisible
  boron_log->SetVisAttributes(G4VisAttributes::Invisible);

}

//#####################################################################
// The NaI Annulus
void DetectorConstruction::BuildNaIAnnulus(){

//--------------- NaI Annulus --------------------
  // Annulus housing
  // - inner housing
  G4double innerNaIHouseInner = (119.38/2.)-0.4*mm;
  G4double outerNaIHouseInner = (119.38/2.)+0.4*mm;
  G4double HlengthNaIHouseInner = 361.95/2.*mm;
  G4RotationMatrix noRotate;

  G4Tubs* NaIHouseInner_solid
    = new G4Tubs("NaI Inner Housing", innerNaIHouseInner,outerNaIHouseInner,
		 HlengthNaIHouseInner, 0.*deg,360.*deg);
  G4LogicalVolume* NaIHouseInner_log = 
    new G4LogicalVolume(NaIHouseInner_solid, Al,"NaI Inner Housing log",0,0,0);
  // G4VPhysicalVolume* NaIHouseInner_phys = 
    new G4PVPlacement(0,
		      G4ThreeVector(),
		      NaIHouseInner_log,
		      "NaI Inner Housing",
		      LENALab_log,
		      false,
		      0);

  // - outer housing
  G4double innerNaIHouseOuter = (368.3/2.)-0.4*mm;
  G4double outerNaIHouseOuter = (368.3/2.)+0.4*mm;
  G4double HlengthNaIHouseOuter = 311.15/2.*mm;

  G4Tubs* NaIHouseOuter_solid
    = new G4Tubs("NaI Outer Housing", innerNaIHouseOuter,outerNaIHouseOuter,
		 HlengthNaIHouseOuter, 0.*deg,360.*deg);
  G4LogicalVolume* NaIHouseOuter_log =
    new G4LogicalVolume(NaIHouseOuter_solid, Al,"NaI Outer Housing log",0,0,0);
  // G4VPhysicalVolume* NaIHouseOuter_phys =
    new G4PVPlacement(0,
		      G4ThreeVector(),
		      NaIHouseOuter_log,
		      "NaI Outer Housing",
		      LENALab_log,
		      false,
		      0);

  // - flange
  G4double outerNaIFlange = 412.75/2.*mm;
  G4double innerNaIFlange = innerNaIHouseOuter;
  G4double HlengthNaIFlange = 25.4/2.*mm;
  G4double placementNaIFlange = 168.275*mm;

  G4Tubs* NaIFlange_solid
    = new G4Tubs("NaI Housing Flange", innerNaIFlange,outerNaIFlange,
		 HlengthNaIFlange, 0.*deg,360.*deg);
  G4LogicalVolume* NaIFlange_log =
    new G4LogicalVolume(NaIFlange_solid, Al,"NaI Housing Flange log", 0,0,0);
  // G4VPhysicalVolume* NaIFlange1_phys = 
    new G4PVPlacement(0,
		      G4ThreeVector(0,0,placementNaIFlange),
		      NaIFlange_log,
		      "NaI Housing Flange 1",
		      LENALab_log,
		      false,
		      0);
 // G4VPhysicalVolume* NaIFlange2_phys = 
   new G4PVPlacement(0,
		     G4ThreeVector(0,0,-placementNaIFlange),
		     NaIFlange_log,
		     "NaI Housing Flange 2",
		     LENALab_log,
		     false,
		     0);
  // - fill block
  G4double innerNaIFill = outerNaIHouseInner;
  G4double outerNaIFill = innerNaIFill+47.06;
  G4double HlengthNaIFill = 15.875/2.*mm;
  G4double placementNaIFill = placementNaIFlange+(HlengthNaIFlange-
						  HlengthNaIFill);

  G4Tubs* NaIFill_solid
    = new G4Tubs("NaI Fill Piece", innerNaIFill,outerNaIFill,HlengthNaIFill,
		 0.*deg,360.*deg);
  G4LogicalVolume* NaIFill_log = 
    new G4LogicalVolume(NaIFill_solid, Al,"NaI Fill Piece log",0,0,0);
  // G4VPhysicalVolume*NaIFill1_phys = 
  new G4PVPlacement(0,
		    G4ThreeVector(0,0,placementNaIFill),
		    NaIFill_log,
		    "NaI Fill piece 1",
		    LENALab_log,
		    false,
		    0);
  // G4VPhysicalVolume* NaIFill2_phys = 
  new G4PVPlacement(0,
		    G4ThreeVector(0,0,-placementNaIFill),
		    NaIFill_log,
		    "NaI Fill piece 2",
		    LENALab_log,
		    false,
		    0);
  
  //------ the NaI(Tl) Scintillators
  // some values
  //G4int NoSegments = 16;
  //----------------------------
  // The old sizes
  G4double NaIInnerRad = (60.5-0.1+1.0)*mm;
  G4double NaIOuterRad = (177.8-1.0)*mm;
  G4double NaISegHLength = (323.8/4.-0.1-1.0)*mm;
  
  // make one segment first
  G4double NaISegStartAngle = (22.5+0.5)*deg;
  G4double NaISegSpanningAngle = (44.-1.0)*deg;


  //// make mothers for both halfs
  //  G4Tubs* NaIMother1_solid
  //  = new G4Tubs("NaI Mother 1",NaIInnerRad,NaIOuterRad,
  //		 NaISegHLength,0.*deg,360.*deg);
  //G4LogicalVolume* NaIMother1_log = 
  //  new G4LogicalVolume(NaIMother1_solid, Air,"NaI Mother 1",0,0,0);
  //// G4VPhysicalVolume* NaIMother1_phys =
  //  new G4PVPlacement(0,
  //		      G4ThreeVector(0,0,NaISegHLength+3.2),
  //		      NaIMother1_log,
  //		      "NaI Mother 1",
  //		      NaIHouseOuter_log,
  //		      false,
  //		      0);
  //G4Tubs* NaIMother2_solid
  // = new G4Tubs("NaI Mother 2",NaIInnerRad,NaIOuterRad,
  //		 NaISegHLength,0.*deg,360.*deg);
  //G4LogicalVolume* NaIMother2_log = 
  //  new G4LogicalVolume(NaIMother2_solid, Air, "NaI Mother 2",0,0,0);
  //// G4VPhysicalVolume* NaIMother2_phys =
  //  new G4PVPlacement(0,
  //		      G4ThreeVector(0,0,-NaISegHLength-3.2),
  //		      NaIMother2_log,
  //		      "NaI Mother 2",
  //		      NaIHouseOuter_log,
  //		      false,
  //		      0);

  G4Tubs* NaISegment_solid 
    = new G4Tubs("NaI Segment",NaIInnerRad,NaIOuterRad,
		 NaISegHLength,NaISegStartAngle,NaISegSpanningAngle);
  G4LogicalVolume* NaISeg_log = 
    new G4LogicalVolume(NaISegment_solid, NaITl,"NaI Segment Solid", 0,0,0);

  G4double NaIsegRotate;
  G4ThreeVector zAxis(0,0,1.);


  // The segments are now have labels 0-15. 
  // 0-7 One side
  // 8-15 Opposite side
  // Segment 0 is adjacent to seg 8, 1 next to 9 etc.
  for(G4int i=0;i<8;i++){
    NaIsegRotate=(i*45.+22.5+0.5)*deg;
    new G4PVPlacement(G4Transform3D(G4RotationMatrix(zAxis,NaIsegRotate),
				    G4ThreeVector(0,0,-NaISegHLength-3.2)),
		      NaISeg_log,
		      "NaI Segments 1",
		      LENALab_log,
		      false,
		      i);
  }

  for(G4int i=0;i<8;i++){
    NaIsegRotate=(i*45.+22.5+0.5)*deg;
    new G4PVPlacement(G4Transform3D(G4RotationMatrix(zAxis,NaIsegRotate),
				    G4ThreeVector(0,0,NaISegHLength+3.2)),
		      NaISeg_log,
		      "NaI Segments 2",
		      LENALab_log,
		      false,
		      i+8);
  }

  // Aluminium surrounding NaI segments
  G4double AlNaIInnerRad = NaIInnerRad-1.0*mm;
  G4double AlNaIOuterRad = NaIOuterRad+1.0*mm;
  G4double AlNaISegHLength = NaISegHLength+1.0*mm;
  
  // make one segment first
  G4double AlNaISegStartAngle = NaISegStartAngle;
  G4double AlNaISegSpanningAngle = NaISegSpanningAngle+1.0*deg;

  G4Tubs* AlNaISegment_solid 
    = new G4Tubs("Al NaI Segment",AlNaIInnerRad,AlNaIOuterRad,
		 AlNaISegHLength,AlNaISegStartAngle,AlNaISegSpanningAngle);

  G4SubtractionSolid* OpticalIsolationNaI_solid
    = new G4SubtractionSolid("NaI Optical Isolation",AlNaISegment_solid,
			     NaISegment_solid,&noRotate,G4ThreeVector(0,0,0));
    G4LogicalVolume* OpticalIsolationNaI_log = 
    new G4LogicalVolume(OpticalIsolationNaI_solid, Al,"Optical Isolation",
			0,0,0);

  //------------------------------------------------------------
  //-------------- Visualization Stuff -----------------------
  //------------------------------------------------------------

  // Make any aluminium off White/grey
  G4VisAttributes* AlVisAtt = new G4VisAttributes(G4Colour(0.8,0.8,0.8));
  NaIHouseInner_log->SetVisAttributes(AlVisAtt);
  NaIHouseOuter_log->SetVisAttributes(AlVisAtt);
  NaIFlange_log->SetVisAttributes(AlVisAtt);
  NaIFill_log->SetVisAttributes(AlVisAtt);

  // The Detectors are green
  G4VisAttributes* DetVisAtt = new G4VisAttributes(G4Colour(0.0,0.75,0.0));
  NaISeg_log->SetVisAttributes(DetVisAtt);

}


//#####################################################################
// The Gamma Cube
void DetectorConstruction::BuildGammaCube(){

  // ------------ The Lead Shield ------------
  G4double LeadHthick = 25.4/4.*mm;
  G4double LeadBoxHeight = 539.75 + (2.*LeadHthick);
  // Top
  G4double HlenLeadTopz = (930./2.)-2*25.4;
  G4double HlenLeadTopx = 457.2/2.+2*LeadHthick;
  G4double heightLeadTop= 245.872 + LeadHthick;
  G4RotationMatrix noRotate;

  G4Box* LeadBoxTop_solid
    = new G4Box("Lead Box Top", HlenLeadTopx,LeadHthick,HlenLeadTopz);
  G4LogicalVolume* LeadBoxTop_log = 
    new G4LogicalVolume(LeadBoxTop_solid, Pb,"Leab Box Top log", 0,0,0);
  // G4VPhysicalVolume* LeadBoxTop_phys = 
    new G4PVPlacement(0,
		      G4ThreeVector(0,heightLeadTop,0),
		      LeadBoxTop_log,
		      "Lead Box Top",
		      LENALab_log,
		      false,
		      0);

  // sides
  G4double HlenLeadSidez = HlenLeadTopz;
  G4double HlenLeadSidey = LeadBoxHeight/2. - LeadHthick;
  G4double xposLeadSide  = HlenLeadTopx - LeadHthick;
  G4double yposLeadSide  = HlenLeadSidey-heightLeadTop+LeadHthick;

  G4Box* LeadBoxSide_solid
    = new G4Box("Lead Box Side",LeadHthick,HlenLeadSidey,HlenLeadSidez);
  G4LogicalVolume* LeadBoxSide_log = 
    new G4LogicalVolume(LeadBoxSide_solid, Pb,"Leab Box Side log", 0,0,0);
  // G4VPhysicalVolume* LeadBoxSide1_phys =
    new G4PVPlacement(0,
		      G4ThreeVector(xposLeadSide,-yposLeadSide,0),
		      LeadBoxSide_log,
		      "Lead Box Side 1",
		      LENALab_log,
		      false,
		      0);
  // G4VPhysicalVolume* LeadBoxSide2_phys =
    new G4PVPlacement(0,
		      G4ThreeVector(-xposLeadSide,-yposLeadSide,0),
		      LeadBoxSide_log,
		      "Lead Box Side 2",
		      LENALab_log,
		      false,
		      0);

  // The Detector and Beam Plate 
  G4double HlenLeadBeamx = HlenLeadTopx-2*LeadHthick;
  G4double HlenLeadBeamy = HlenLeadSidey;
  G4double HlenLeadBeamz = LeadHthick;
  G4double BeamPlatePosz = HlenLeadTopz - LeadHthick;
  G4double BeampipeRad = 95.25/2.*mm;
  G4double DetectorRad = 120.65/2.*mm;
  G4ThreeVector MoveDown(0,yposLeadSide,0);

  G4Box* LeadBeamSolid
    = new G4Box("End plates with no holes",HlenLeadBeamx,HlenLeadBeamy,
		HlenLeadBeamz);
  G4Tubs* BeamHole
    = new G4Tubs("Hole for Beam Pipe",0.*mm,BeampipeRad,LeadHthick+1.,
		 0,twopi);
  G4Tubs* DetHole
    = new G4Tubs("Hole for Detector Plate",0.*mm,DetectorRad,LeadHthick+1.,
		 0,twopi);

  G4SubtractionSolid* LeadBoxBeam_solid
    = new G4SubtractionSolid("Beam Pipe Plate",LeadBeamSolid,BeamHole,
			     &noRotate,MoveDown);
  G4SubtractionSolid* LeadBoxDet_solid
    = new G4SubtractionSolid("Detector Pipe Plate",LeadBeamSolid,DetHole,
			     &noRotate,MoveDown);

  G4LogicalVolume* LeadBoxBeam_log = 
    new G4LogicalVolume(LeadBoxBeam_solid, Pb,"Beam Pipe Plate log",0,0,0);
  G4LogicalVolume* LeadBoxDet_log = 
    new G4LogicalVolume(LeadBoxDet_solid, Pb,"Detector Pipe Plate log",0,0,0);
  // G4VPhysicalVolume* LeadBoxBeam_phys = 
    new G4PVPlacement(0,
		      G4ThreeVector(0,-yposLeadSide,-BeamPlatePosz),
		      LeadBoxBeam_log,
		      "Lead Box Beam Side",
		      LENALab_log,
		      false,
		      0);
  // G4VPhysicalVolume* LeadBoxDet_phys = 
    new G4PVPlacement(0,
		      G4ThreeVector(0,-yposLeadSide,BeamPlatePosz),
		      LeadBoxDet_log,
		      "Lead Box Detector Side",
		      LENALab_log,
		      false,
		      0);

  // ---------- Aluminium holding up lead -----------
  G4double AlHthick = 5.*mm;
  G4double AlBoxHeight = LeadBoxHeight-2*LeadHthick;

  // Top
  G4double HlenAlTopz = HlenLeadTopz-2*LeadHthick;
  G4double HlenAlTopx = HlenLeadTopx-2*LeadHthick;
  G4double heightAlTop= heightLeadTop-LeadHthick-AlHthick;

  G4Box* AlBoxTop_solid
    = new G4Box("Al Box Top", HlenAlTopx,AlHthick,HlenAlTopz);
  G4LogicalVolume* AlBoxTop_log = 
    new G4LogicalVolume(AlBoxTop_solid, Al,"Al Box Top log", 0,0,0);
  // G4VPhysicalVolume* AlBoxTop_phys = 
    new G4PVPlacement(0,
		      G4ThreeVector(0,heightAlTop,0),
		      AlBoxTop_log,
		      "Al Box Top",
		      LENALab_log,
		      false,
		      0);

  // sides
  G4double HlenAlSidez = HlenAlTopz;
  G4double HlenAlSidey = AlBoxHeight/2. - AlHthick;
  G4double xposAlSide  = HlenAlTopx - AlHthick;
  G4double yposAlSide  = HlenAlSidey-heightAlTop+AlHthick;

  G4Box* AlBoxSide_solid
    = new G4Box("Al Box Side",AlHthick,HlenAlSidey,HlenAlSidez);
  G4LogicalVolume* AlBoxSide_log = 
    new G4LogicalVolume(AlBoxSide_solid, Al,"Al Box Side log", 0,0,0);
  // G4VPhysicalVolume* AlBoxSide1_phys = 
    new G4PVPlacement(0,
		      G4ThreeVector(+xposAlSide,-yposAlSide,0),
		      AlBoxSide_log,
		      "Al Box Side 1",
		      LENALab_log,
		      false,
		      0);
  // G4VPhysicalVolume* AlBoxSide2_phys = 
    new G4PVPlacement(0,
		      G4ThreeVector(-xposAlSide,-yposAlSide,0),
		      AlBoxSide_log,
		      "Al Box Side 2",
		      LENALab_log,
		      false,
		      0);

  // The Detector and Beam Plate 
  G4double HlenAlBeamx = HlenAlTopx-2*AlHthick;
  G4double HlenAlBeamy = HlenAlSidey;
  G4double HlenAlBeamz = AlHthick;
  G4double AlBeamPlatePosz = BeamPlatePosz-LeadHthick-AlHthick;
  G4ThreeVector MoveDownAl(0,yposAlSide,0);

  G4Box* AlBeamSolid
    = new G4Box("End Al plates with no holes",HlenAlBeamx,HlenAlBeamy,
		HlenAlBeamz);

  G4SubtractionSolid* AlBoxBeam_solid
    = new G4SubtractionSolid("Al Beam Pipe Plate",AlBeamSolid,BeamHole,
			     &noRotate,MoveDownAl);
  G4SubtractionSolid* AlBoxDet_solid
    = new G4SubtractionSolid("Al Detector Pipe Plate",AlBeamSolid,DetHole,
			     &noRotate,MoveDownAl);

  G4LogicalVolume* AlBoxBeam_log =
    new G4LogicalVolume(AlBoxBeam_solid, Al,"Beam Pipe Plate log",0,0,0);
  G4LogicalVolume* AlBoxDet_log = 
    new G4LogicalVolume(AlBoxDet_solid, Al,"Detector Pipe Plate log",0,0,0);
  // G4VPhysicalVolume* AlBoxBeam_phys = 
    new G4PVPlacement(0,
		      G4ThreeVector(0,-yposAlSide,AlBeamPlatePosz),
		      AlBoxBeam_log,
		      "Al Box Beam Side",
		      LENALab_log,
		      false,
		      0);
  // G4VPhysicalVolume* AlBoxDet_phys = 
    new G4PVPlacement(0,
		      G4ThreeVector(0,-yposAlSide,-AlBeamPlatePosz),
		      AlBoxDet_log,
		      "Al Box Detector Side",
		      LENALab_log,
		      false,
		      0);

  // The Base of the Table
  G4double baseLenz = 700*mm;
  G4double baseLenx = 600*mm;
  G4double AlbaseHThick = 1.75*25.4/2.*mm;
  G4double AlbasePosy = 350*mm;

G4Box* AlBase_Solid
    = new G4Box("Base Aluminium",baseLenx,AlbaseHThick,baseLenz);
 G4LogicalVolume* AlBase_log = 
   new G4LogicalVolume(AlBase_Solid,Al,"Base Aluminium Log",0,0,0);
 // G4VPhysicalVolume* AlBase_phys = 
   new G4PVPlacement(0,
		     G4ThreeVector(0,-AlbasePosy,0),
		     AlBase_log,
		     "Base Aluminium",
		     LENALab_log,
		     false,
		     0);

 G4double LeadbaseHThick = 0.5*25.4/2.*mm;
 G4double LeadbasePosy = AlbasePosy + AlbaseHThick + LeadbaseHThick;
 G4Box* LeadBase_Solid
    = new G4Box("Base Lead",baseLenx,LeadbaseHThick,baseLenz);
 G4LogicalVolume* LeadBase_log = 
   new G4LogicalVolume(LeadBase_Solid,Pb,"Base Lead Log",0,0,0);
 // G4VPhysicalVolume* LeadBase_phys =
   new G4PVPlacement(0,
		     G4ThreeVector(0,-LeadbasePosy,0),
		     LeadBase_log,
		     "Base Lead",
		     LENALab_log,
		     false,
		     0);
 
  // ------------- Scintillator Shield ---------
  G4double SciHthick = 25.4*mm;
  G4double SciBoxHeight = 558.8 + (2.*SciHthick);

  // Top
  G4double HlenSciTopz = 930./2.*mm;
  G4double HlenSciTopx = 576/2.*mm;
  G4double heightSciTop= 264.922+SciHthick*mm;

  G4Box* SciBoxTop_solid
    = new G4Box("Sci Box Top", HlenSciTopx,SciHthick,HlenSciTopz);
  G4LogicalVolume* SciBoxTop_log = 
    new G4LogicalVolume(SciBoxTop_solid, Sci,"Leab Box Top log", 0,0,0);
  // G4VPhysicalVolume* SciBoxTop_phys = 
    new G4PVPlacement(0,
		      G4ThreeVector(0,heightSciTop,0),
		      SciBoxTop_log,
		      "Sci Box Top",
		      LENALab_log,
		      false,
		      0);

  // sides
  G4double HlenSciSidez = HlenSciTopz-SciHthick;
  G4double HlenSciSidey = SciBoxHeight/2. - SciHthick;
  G4double zposSciSide  = SciHthick;
  G4double xposSciSide  = xposLeadSide + LeadHthick + SciHthick;
  G4double yposSciSide  = HlenSciSidey-heightSciTop+SciHthick;

  G4Box* SciBoxSide_solid
    = new G4Box("Sci Box Side",SciHthick,HlenSciSidey,HlenSciSidez);
  G4LogicalVolume* SciBoxSide_log = 
    new G4LogicalVolume(SciBoxSide_solid, Sci,"Sci Box Side log", 0,0,0);
  // G4VPhysicalVolume* SciBoxSide1_phys = 
    new G4PVPlacement(0,
		      G4ThreeVector(xposSciSide+3.,-yposSciSide,
				    zposSciSide),
		      SciBoxSide_log,
		      "Sci Box Side 1",
		      LENALab_log,
		      false,
		      0);
  // G4VPhysicalVolume* SciBoxSide2_phys = 
    new G4PVPlacement(0,
		      G4ThreeVector(-xposSciSide-3.,-yposSciSide,
				    -zposSciSide),
		      SciBoxSide_log,
		      "Sci Box Side 2",
		      LENALab_log,
		      false,
		      0);

  // The Detector and Beam Plate 
  G4double HlenSciBeamx = HlenSciTopx-SciHthick;
  G4double HlenSciBeamy = HlenSciSidey;
  G4double HlenSciBeamz = SciHthick;
  G4double SciBeamPlatePosz = HlenSciTopz - SciHthick;
  G4double SciBeampipeRad = 98.425/2.*mm;
  G4double SciDetectorRad = 137.16/2.*mm;
  G4ThreeVector MoveSciDet(SciHthick,yposSciSide+8.255,0);
  G4ThreeVector MoveSciBeam(-SciHthick,yposSciSide,0);

  G4Box* SciBeamSolid
    = new G4Box("End plates with no holes",HlenSciBeamx,HlenSciBeamy,
		HlenSciBeamz);
  G4Tubs* SciBeamHole
    = new G4Tubs("Hole for Beam Pipe",0.*mm,SciBeampipeRad,SciHthick+1.,
		 0,twopi);
  G4Tubs* SciDetHole
    = new G4Tubs("Hole for Detector Plate",0.*mm,SciDetectorRad,SciHthick+1.,
		 0,twopi);

  G4SubtractionSolid* SciBoxBeam_solid
    = new G4SubtractionSolid("Beam Pipe Plate",SciBeamSolid,SciBeamHole,
			     &noRotate,MoveSciBeam);
  G4SubtractionSolid* SciBoxDet_solid
    = new G4SubtractionSolid("Detector Pipe Plate",SciBeamSolid,SciDetHole,
			     &noRotate,MoveSciDet);

  G4LogicalVolume* SciBoxBeam_log = 
    new G4LogicalVolume(SciBoxBeam_solid, Sci,"Beam Pipe Plate log",0,0,0);
  G4LogicalVolume* SciBoxDet_log = 
    new G4LogicalVolume(SciBoxDet_solid, Sci,"Detector Pipe Plate log",0,0,0);
  // G4VPhysicalVolume* SciBoxBeam_phys = 
    new G4PVPlacement(0,
		      G4ThreeVector(SciHthick,-yposSciSide,
				    -SciBeamPlatePosz-6.),
		      SciBoxBeam_log,
		      "Sci Box Beam Side",
		      LENALab_log,
		      false,
		      0);
  // G4VPhysicalVolume* SciBoxDet_phys = 
    new G4PVPlacement(0,
		      G4ThreeVector(-SciHthick,-yposSciSide,
				    SciBeamPlatePosz+6.),
		      SciBoxDet_log,
		      "Sci Box Detector Side",
		      LENALab_log,
		      false,
		      0);


  //------------------------------------------------------------
  //-------------- Visualization Stuff -----------------------
  //------------------------------------------------------------

  G4VisAttributes* AlBoxVisAtt = new G4VisAttributes(G4Colour(0.8,0.8,0.8));
  //AlBoxVisAtt->SetForceWireframe(true);
  AlBoxTop_log->SetVisAttributes(AlBoxVisAtt);
  AlBoxSide_log->SetVisAttributes(AlBoxVisAtt);
  AlBoxBeam_log->SetVisAttributes(AlBoxVisAtt);
  //AlBoxBeam_log->SetVisAttributes(G4VisAttributes::Invisible);
  AlBoxDet_log->SetVisAttributes(AlBoxVisAtt);
  //AlBoxDet_log->SetVisAttributes(G4VisAttributes::Invisible);

  // Scintillator Detectors are Black
  G4VisAttributes* SciVisAtt = new G4VisAttributes(G4Colour(0.1,0.1,0.1));
  //  SciVisAtt->SetForceWireframe(true);
  SciBoxTop_log->SetVisAttributes(SciVisAtt);
  SciBoxSide_log->SetVisAttributes(SciVisAtt);
  SciBoxBeam_log->SetVisAttributes(SciVisAtt);
  //SciBoxBeam_log->SetVisAttributes(G4VisAttributes::Invisible);
  SciBoxDet_log->SetVisAttributes(SciVisAtt);
  //SciBoxDet_log->SetVisAttributes(G4VisAttributes::Invisible);

  // Lead is Red for now
  G4VisAttributes* PbVisAtt = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
  //  PbVisAtt->SetForceWireframe(true);
  LeadBoxTop_log->SetVisAttributes(PbVisAtt);
  LeadBoxSide_log->SetVisAttributes(PbVisAtt);
  LeadBoxBeam_log->SetVisAttributes(PbVisAtt);
  //LeadBoxBeam_log->SetVisAttributes(G4VisAttributes::Invisible);
  LeadBoxDet_log->SetVisAttributes(PbVisAtt);
  //LeadBoxDet_log->SetVisAttributes(G4VisAttributes::Invisible);
  LeadBase_log->SetVisAttributes(PbVisAtt);

}


//#####################################################################
// The APEX Detector
void DetectorConstruction::BuildAPEX(){

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
  G4ThreeVector ssteelslab(0,3.44*cm,0);
  G4ThreeVector cyl_holeL(0,0,-28.05*cm);
  G4ThreeVector cyl_holeR(0,0,28.05*cm);
  G4ThreeVector APEXradius(0,-24.5*cm,0);
  G4ThreeVector AngCorradius(0,0,24.5*cm);
  G4RotationMatrix noRotate;

  //----------------Scintillator Casing--------------

  G4Box* scintcasing_box
    = new G4Box("SSteelSlab",HalfSlabWidth,HalfSlabHeight,HalfSlabLength);

  G4Tubs* cyl_hole
    = new G4Tubs("CylindricalHole",0.*mm,QuartzRad,QuartzThick,0,twopi);

  G4Trap* ssteelcasing_trap
    = new G4Trap("NaIcasing",HalfTrapLength+10.0*mm,Theta,Phi,HalfTrapHeight+0.4*mm,
		HalfTrapWidthShort+0.4*mm,HalfTrapWidthLong+0.5*mm,
		Alpha,HalfTrapHeight+0.4*mm,HalfTrapWidthShort+0.4*mm,
		HalfTrapWidthLong+0.5*mm,Alpha);

  G4Box* volume_box
    = new G4Box("SlabVolume",HalfSlabWidth,HalfSlabHeight,HalfSlabLength+2.0*mm);

  G4Trap* volume_trap
    = new G4Trap("TrapVolume",HalfTrapLength+12.0*mm,Theta,Phi,HalfTrapHeight+0.4*mm,
		HalfTrapWidthShort+0.4*mm,HalfTrapWidthLong+0.5*mm,
		Alpha,HalfTrapHeight+0.4*mm,HalfTrapWidthShort+0.4*mm,
		HalfTrapWidthLong+0.5*mm,Alpha);

  // subtract material for quartz windows
  G4SubtractionSolid* ssteelcasing_trap0
    = new G4SubtractionSolid("Quartz window hole L",ssteelcasing_trap,cyl_hole,
			     &noRotate,cyl_holeL);
  G4SubtractionSolid* scintcasing_trap
    = new G4SubtractionSolid("Quartz window hole R",ssteelcasing_trap0,cyl_hole,
			     &noRotate,cyl_holeR);

  // now join the SSteelSlab to the main NaIcasing
  G4UnionSolid* scintcasing
    = new G4UnionSolid("SSteel NaI Casing",scintcasing_trap,scintcasing_box,
		       &noRotate,ssteelslab);

  // now join the SlabVolume to the main TrapVolume
  G4UnionSolid* scintvolume
    = new G4UnionSolid("Casing Volume",volume_trap,volume_box,
		       &noRotate,ssteelslab);

  G4LogicalVolume* scintcasing_log
    = new G4LogicalVolume(scintcasing, ssteel, "casing_log", 0, 0, 0);

  G4LogicalVolume* scintvolume_log
    = new G4LogicalVolume(scintvolume, Air, "volume_log", 0, 0, 0);

  // place SSteel NaI Casing inside the Casing Volume
  new G4PVPlacement(0,
                        G4ThreeVector(),		// at (0,0,0)
                        scintcasing_log,           	// Logical Volume
                        "Segment Volume",           	// Name
                        scintvolume_log,           	// Mother Volume
                        false,                     	// no boolean operations
                        0);                        	// Copy number


  for(G4int i=0;i<24;i++){
    NaIsegRotate=(i*(360/24))*deg;
    new G4PVPlacement(G4Transform3D(G4RotationMatrix(zAxis,-NaIsegRotate),
                      G4ThreeVector(Radius*std::sin(NaIsegRotate),Radius*std::cos(NaIsegRotate),0)),
		      scintvolume_log,		// Logical Volume
		      "NaI Segment ",		// Name
		      LENALab_log,		// Mother Volume
		      false,			// no boolean operations
		      i);			// Copy number
}
/*
  // now place the Casing Volume inside the LENALab mother volume
  // NOTE: rotation and position are important
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(zAxis,180*deg),
                        APEXradius),			// at radius of APEX
                        scintvolume_log,           	// Logical Volume
                        "NaI Segment 1",           	// Name
                        LENALab_log,           	   	// Mother Volume
                        false,                     	// no boolean operations
                        0);                        	// Copy number

  // now place the Casing Volume inside the LENALab mother volume
  // NOTE: rotation and position are important
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(180*deg,90*deg,90*deg),
                        AngCorradius),			// at radius of angular correlation measurement
                        scintvolume_log,           	// Logical Volume
                        "NaI Segment 1",           	// Name
                        LENALab_log,           	   	// Mother Volume
                        false,                     	// no boolean operations
                        0);                        	// Copy number
*/

  //------------------Quartz Windows-----------------  

  G4Tubs* quartz_cyl
    = new G4Tubs("QuartzWindow",0.*mm,QuartzRad,QuartzThick,0,twopi);

  G4LogicalVolume* quartz_log
    = new G4LogicalVolume(quartz_cyl, Quartz,"quartz_log",0,0,0);

  new G4PVPlacement(0,                             		// no rotation
                        cyl_holeL,				// at (0,0,-28.05*cm)
                        quartz_log,          	   		// Logical Volume
                        "quartzwindowL",            		// Name
                        scintvolume_log,	           	// Mother Volume
                        false,                     		// no boolean operations
                        0);                        		// Copy number

  new G4PVPlacement(0,                             		// no rotation
                        cyl_holeR,				// at (0,0,28.05*cm)
                        quartz_log,          	   		// Logical Volume
                        "quartzwindowR",            		// Name
                        scintvolume_log,	           	// Mother Volume
                        false,                     		// no boolean operations
                        1);                        		// Copy number

  //------------------Scintillators------------------

  G4Trap* scintillator_trap
    = new G4Trap("NaIcrystal",HalfTrapLength,Theta,Phi,HalfTrapHeight,
			HalfTrapWidthShort,HalfTrapWidthLong,Alpha,HalfTrapHeight,
			HalfTrapWidthShort,HalfTrapWidthLong,Alpha);

//  G4LogicalVolume* scintillator_log
	scintillator_log
    = new G4LogicalVolume(scintillator_trap, NaITl, "scint_log", 0, 0, 0);

  //G4VPhysicalVolume* scintillator_phys
  new G4PVPlacement(0,                             // no rotation
                        G4ThreeVector(),           // at (0,0,0)
                        scintillator_log,          // Logical Volume
                        "scintillator",            // Name
                        scintcasing_log,           // Mother Volume
                        false,                     // no boolean operations
                        0);                        // Copy number

  //------------------------------------------------------------
  //-------------- Visualization Stuff -----------------------
  //------------------------------------------------------------

  // The scintillator volume is invisible
  scintvolume_log ->SetVisAttributes(G4VisAttributes::Invisible);

  // Scintillator casing is dark gray
  G4VisAttributes* ssteelVisAtt = new G4VisAttributes(G4Colour(0.4,0.4,0.4));
  scintcasing_log->SetVisAttributes(ssteelVisAtt);

  // NaI(Tl) scintillator is yellow
  G4VisAttributes* detVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  scintillator_log->SetVisAttributes(detVisAtt);

  // Quartz window is white
  G4VisAttributes* quartzVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  quartz_log->SetVisAttributes(quartzVisAtt);

  //scintillator_log->SetSensitiveDetector(NaISD);

}

//#####################################################################
// The APEX Cradle
void DetectorConstruction::BuildAPEXcradle(){

  // The Lead Shield
  G4double LeadShield_outer = 25.75*25.4/2.*mm;
  G4double LeadShield_inner = 24.25*25.4/2.*mm;
  G4double LeadShield_Hlen  = 23.5*25.4/2.*mm;

  G4Tubs* LeadShield_solid
    = new G4Tubs("Lead Shield",LeadShield_inner,LeadShield_outer,
		 LeadShield_Hlen,0,twopi);
  G4LogicalVolume* LeadShield_log =
    new G4LogicalVolume(LeadShield_solid, Pb,"Lead Shield log",0,0,0);

    new G4PVPlacement(0,
		      G4ThreeVector(0,0,0),
		      LeadShield_log,
		      "Lead Shield",
		      LENALab_log,
		      false,
		      0);

  // The Al Cradle
  G4double AlCradle_outer = 26.5*25.4/2.*mm;
  G4double AlCradle_inner = 25.75*25.4/2.*mm;
  G4double AlCradle_Hlen  = 25.75*25.4/2.*mm;

  G4double CradleRad = 26.5*25.4/2.*mm;
  G4double CradleThick = 25.75*25.4/2.*mm;

  G4double AlBoxLength = 0.25*25.4/2.*mm;
  G4double AlBoxHeight = 28.625*25.4/2.*mm;
  G4double AlBoxWidth = 28.5*25.4/2.*mm;

  G4ThreeVector origin(0,0,0);
  G4ThreeVector wingL1(0,0,-11.75*25.4*mm);
  G4ThreeVector wingR1(0,0,11.75*25.4*mm);
  G4ThreeVector wingL2(0,0,-7.5*25.4*mm);
  G4ThreeVector wingR2(0,0,7.5*25.4*mm);
  G4RotationMatrix noRotate;

  G4Box* Al_box
    = new G4Box("AlBox",AlBoxWidth,AlBoxHeight,AlBoxLength);

  G4Tubs* sub_cradle
    = new G4Tubs("CradleSub",0.*mm,CradleRad,CradleThick,0,twopi);


  // main cradle
  G4Tubs* AlCradle_solid
    = new G4Tubs("Al Cradle",AlCradle_inner,AlCradle_outer,
		 AlCradle_Hlen,0,twopi);

  G4LogicalVolume* apexcradle_log =
    new G4LogicalVolume(AlCradle_solid, Al,"cradle_log",0,0,0);

    new G4PVPlacement(0,
		      G4ThreeVector(0,0,0),
		      apexcradle_log,
		      "Al Cradle",
		      LENALab_log,
		      false,
		      0);


  // subtract material for al cradle wings
  G4SubtractionSolid* cradlewings_solid
    = new G4SubtractionSolid("Al wings",Al_box,sub_cradle,
			     &noRotate,origin);

  G4LogicalVolume* cradlewings_log =
    new G4LogicalVolume(cradlewings_solid, Al,"Alwings_log",0,0,0);

    new G4PVPlacement(0,
		      wingL1,
		      cradlewings_log,
		      "Wing L1",
		      LENALab_log,
		      false,
		      0);
    new G4PVPlacement(0,
		      wingR1,
		      cradlewings_log,
		      "Wing R1",
		      LENALab_log,
		      false,
		      1);
    new G4PVPlacement(0,
		      wingL2,
		      cradlewings_log,
		      "Wing L2",
		      LENALab_log,
		      false,
		      2);
    new G4PVPlacement(0,
		      wingR2,
		      cradlewings_log,
		      "Wing R2",
		      LENALab_log,
		      false,
		      1);

  //------------------------------------------------------------
  //-------------- Visualization Stuff -----------------------
  //------------------------------------------------------------

  // Al Cradle is Blue
  G4VisAttributes* CradleVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  apexcradle_log->SetVisAttributes(CradleVisAtt);
  cradlewings_log->SetVisAttributes(CradleVisAtt);

  // Lead is Black
  G4VisAttributes* LeadVisAtt = new G4VisAttributes(G4Colour(0.1,0.1,0.1));
  LeadShield_log->SetVisAttributes(LeadVisAtt);

}

//#####################################################################
// The Target Holder
void DetectorConstruction::BuildTargetHolder(){

  //---------------- The Beam Pipe --------------------------

  // The Stainless Steel Beam Pipe
  G4double BeamPipe_outer = 2.375*25.4/2.*mm;
  G4double BeamPipe_inner = 2.05*25.4/2.*mm;
  G4double BeamPipe_Hlen  = 10*25.4*mm;
  G4double BeamPipe_posz  = -BeamPipe_Hlen+0.5*mm; //+ 0.375*25.4;
  G4RotationMatrix noRotate;

  G4Tubs* BeamPipe_solid
    = new G4Tubs("Beam Pipe",BeamPipe_inner,BeamPipe_outer,
		 BeamPipe_Hlen,0,twopi);
  G4LogicalVolume* BeamPipe_log =
    new G4LogicalVolume(BeamPipe_solid, ssteel,"Beam Pipe log",0,0,0);
  // G4VPhysicalVolume* BeamPipe_phys = 
    new G4PVPlacement(0,
		      G4ThreeVector(0,0,BeamPipe_posz),
		      BeamPipe_log,
		      "Beam Pipe",
		      LENALab_log,
		      false,
		      0);

  // The Copper Cold Trap comes in two pieces: main and end
  G4double BeamPipeCuEnd_Hlen  = 36.5/2.*mm;
  G4double BeamPipeCuEnd_outer = 34.9/2.*mm;
  G4double BeamPipeCuEnd_inner = BeamPipeCuEnd_outer - 4.85*mm;

  G4double BeamPipeCuMain_outer = 1.375*25.4/2.*mm;
  G4double BeamPipeCuMain_inner = 1.25*25.4/2.*mm;
  G4double BeamPipeCuMain_Hlen  = BeamPipe_Hlen-BeamPipeCuEnd_Hlen;

  G4double SuppressorCu_Hlen  = 6.5/2.*mm;
  G4double SuppressorCu_outer = 44.76/2.*mm;
  G4double SuppressorCu_inner = 19.3/2.*mm;
  G4double SuppressorCeramic_Hlen = 14.15/2.*mm - SuppressorCu_Hlen;
  G4double SuppressorCeramic_outer = SuppressorCu_outer;
  G4double SuppressorCeramic_inner = 27.97/2.*mm;

  G4double BeamPipeCuMain_posz  = (BeamPipe_posz-
				   BeamPipeCuEnd_Hlen-
				   2.*(SuppressorCu_Hlen+SuppressorCeramic_Hlen)-
				   0.25*25.4*mm);
  G4double BeamPipeCuEnd_posz  = (BeamPipeCuMain_posz+
				  BeamPipeCuMain_Hlen+
				  BeamPipeCuEnd_Hlen);
  G4double SuppressorCeramic_posz = (BeamPipeCuEnd_posz + 
				     BeamPipeCuEnd_Hlen +
				     SuppressorCeramic_Hlen);
  G4double SuppressorCu_posz = (SuppressorCeramic_posz + 
				SuppressorCeramic_Hlen +
				SuppressorCu_Hlen);

  G4Tubs* BeamPipeCuMain_solid
    = new G4Tubs("Main Cu Beam Pipe",BeamPipeCuMain_inner,BeamPipeCuMain_outer,
		 BeamPipeCuMain_Hlen,0,twopi);
  G4LogicalVolume* BeamPipeCuMain_log = 
    new G4LogicalVolume(BeamPipeCuMain_solid, Cu,"Main Cu Beam Pipe log",0,0,0);
  // G4VPhysicalVolume* BeamPipeCuMain_phys =
    new G4PVPlacement(0,
		      G4ThreeVector(0,0,BeamPipeCuMain_posz),
		      BeamPipeCuMain_log,
		      "Main Cu Beam Pipe",
		      LENALab_log,
		      false,
		      0);

  G4Tubs* BeamPipeCuEnd_solid
    = new G4Tubs("End of Cu Beam Pipe",BeamPipeCuEnd_inner,BeamPipeCuEnd_outer,
		 BeamPipeCuEnd_Hlen,0,twopi);
  G4LogicalVolume* BeamPipeCuEnd_log = 
    new G4LogicalVolume(BeamPipeCuEnd_solid, Cu,"End of Cu Beam Pipe log",0,0,0);
  // G4VPhysicalVolume* BeamPipeCuEnd_phys =
  new G4PVPlacement(0,
		      G4ThreeVector(0,0,BeamPipeCuEnd_posz),
		      BeamPipeCuEnd_log,
		      "End of Cu Beam Pipe",
		      LENALab_log,
		      false,
		      0);

  G4Tubs* SuppressorCeramic_solid
    = new G4Tubs("Ceramic in Suppressor",SuppressorCeramic_inner,SuppressorCeramic_outer,
		 SuppressorCeramic_Hlen,0,twopi);
  G4LogicalVolume* SuppressorCeramic_log = 
    new G4LogicalVolume(SuppressorCeramic_solid, ceramic,"Ceramic in Suppressor log",0,0,0);
  // G4VPhysicalVolume* SuppressorCeramic_phys =
    new G4PVPlacement(0,
		      G4ThreeVector(0,0,SuppressorCeramic_posz),
		      SuppressorCeramic_log,
		      "Ceramic in Suppressor",
		      LENALab_log,
		      false,
		      0);

  G4Tubs* SuppressorCu_solid
    = new G4Tubs("Cu in Suppressor",SuppressorCu_inner,SuppressorCu_outer,
		 SuppressorCu_Hlen,0,twopi);

  // Cone to subtract
  G4double SuppressorCone_largerad = SuppressorCeramic_inner;
  G4double SuppressorCone_smallrad = SuppressorCu_inner;
    // Cone
  G4Cons* SuppressorCone_solid
    = new G4Cons("Subtracted part of suppressor",
		 0,                      // inner rad at -Hlength
		 SuppressorCone_largerad,     // outer rad at -Hlength
		 0,                      // inner rad at HLength
		 SuppressorCone_smallrad,     // outer rad at HLength
		 SuppressorCu_Hlen+0.1*mm,    // half length
		 0,twopi);   
  // The subtraction
  G4SubtractionSolid* SuppressorKnife_solid
    = new G4SubtractionSolid("Suppressor Cu knife edge",SuppressorCu_solid,
			     SuppressorCone_solid,
			     &noRotate,G4ThreeVector(0,0,0));
  G4LogicalVolume* SuppressorCu_log = 
    new G4LogicalVolume(SuppressorKnife_solid, Cu,"Cu in Suppressor log",0,0,0);

  // G4VPhysicalVolume* SuppressorCu_phys =
    new G4PVPlacement(0,
		      G4ThreeVector(0,0,SuppressorCu_posz),
		      SuppressorCu_log,
		      "Cu in Suppressor",
		      LENALab_log,
		      false,
		      0);

  // The Tantalum Target Backing
  G4double TantBack_outer = 1.85*25.4/2.*mm;
  G4double TantBack_Hthick = 0.25*mm;
  
  G4Tubs* TantBack_solid
    = new G4Tubs("Tantalum Backing",0.*mm,TantBack_outer,
		 TantBack_Hthick,0,twopi);
  G4LogicalVolume* TantBack_log = 
    new G4LogicalVolume(TantBack_solid, Ta,"Tantalum",0,0,0);
  // G4VPhysicalVolume* TantBack_phys = 
    new G4PVPlacement(0,
		      G4ThreeVector(0,0,TantBack_Hthick),
		      TantBack_log,
		      "Tantalum",
		      LENALab_log,
		      false,
		      0);

  // The Back piece that buts up against the detector
  G4double THolder_outer = 1.285*25.4/2.*mm;
  G4double THolder_Hthick = 0.075*25.4/2.*mm;
  G4double THolder_posz = 0.375*25.4 - THolder_Hthick+2.*TantBack_Hthick;

  G4Tubs* THolder_solid
    = new G4Tubs("Back of Target Holder",0.*mm,THolder_outer,THolder_Hthick,
		 0,twopi);
  G4LogicalVolume* THolder_log = 
    new G4LogicalVolume(THolder_solid,ssteel,"Target Holder log",0,0,0);
  // G4VPhysicalVolume* THolder_phys = 
    new G4PVPlacement(0,
		      G4ThreeVector(0,0,THolder_posz),
		      THolder_log,
		      "Back of Target Holder",
		      LENALab_log,
		      false,
		      0);

  // The Edges (main frame)
  G4double THolder2_outer = 2.375*25.4/2.*mm;
  G4double THolder2_inner = 1.160*25.4/2.*mm;
  G4double THolder2_Hthick = 0.235*25.4/2.*mm;
  G4double THolder2_posz = THolder2_Hthick+2*TantBack_Hthick;

  G4Tubs* THolderFrame_solid
    = new G4Tubs("Frame of Target Holder",THolder2_inner,THolder2_outer,
		 THolder2_Hthick,0,twopi);

  // Wings for water pipes
  G4double WingHThick = THolder2_Hthick;
  G4double WingHLength = 0.375*25.4/2.*mm;
  G4double WingPosy   = THolder2_outer*cos(asin(WingHLength/THolder2_outer))+WingHLength;
    //THolder2_outer+WingHLength;   // wrt THolder2

  G4ThreeVector MoveWings1(0,WingPosy,0);
  G4ThreeVector MoveWings2(0,-WingPosy,0); 

//  G4ThreeVector MoveWings1(WingPosy,0,0);
//  G4ThreeVector MoveWings2(-WingPosy,0,0); 

  G4Box* Wings_solid
    = new G4Box("Wings",WingHLength,WingHLength,WingHThick);
  // now join the wings to the main target holder frame
  G4UnionSolid* THolder2p_solid
    = new G4UnionSolid("Edges of Target Holder",THolderFrame_solid,Wings_solid,
		       &noRotate,MoveWings1);
  G4UnionSolid* THolder2_solid
    = new G4UnionSolid("Edges of Target Holder",THolder2p_solid,Wings_solid,
		       &noRotate,MoveWings2);

  G4LogicalVolume* THolder2_log = 
    new G4LogicalVolume(THolder2_solid,ssteel,"Target Holder Edges log",0,0,0);
  // G4VPhysicalVolume* THolder2_phys = 
  new G4PVPlacement(0,
		    G4ThreeVector(0,0,THolder2_posz),
		    THolder2_log,
		    "Edges of Target Holder",
		    LENALab_log,
		    false,
		    0);
  
//   G4LogicalVolume* Wings_log =
//     new G4LogicalVolume(Wings_solid,ssteel,"Wings",0,0,0);
//   // G4VPhysicalVolume* Wings_phys = 
//   new G4PVPlacement(0,
// 		    G4ThreeVector(WingPosx,0,THolder2_posz),
// 		    Wings_log,
// 		    "Wings of Target Holder",
// 		    LENALab_log,
// 		    false,
// 		    0);

  // The rest of the Target Holder (to fill in water chamber)
  G4double THolder3_outer = THolder_outer;
  G4double THolder3_inner = THolder2_inner;
  G4double THolder3_Hthick = (0.375*25.4 - //2*TantBack_Hthick - 
			      2*THolder2_Hthick - 2*THolder_Hthick)/2.;
  G4double THolder3_posz = THolder_posz - THolder_Hthick - THolder3_Hthick;

  G4Tubs* THolder3_solid
    = new G4Tubs("Rest of Target Holder",THolder3_inner,THolder3_outer,
		 THolder3_Hthick,0,twopi);
  G4LogicalVolume* THolder3_log = 
    new G4LogicalVolume(THolder3_solid,ssteel,"Target Holder rest log",0,0,0);
  // G4VPhysicalVolume* THolder3_phys =
  new G4PVPlacement(0,
		    G4ThreeVector(0,0,THolder3_posz),
		    THolder3_log,
		    "Rest of Target Holder",
		    LENALab_log,
		    false,
		    0);
  
  // The raised ring
  G4double ringHHeight = 0.125*25.4/2.*mm;
  G4double ringHThick = 0.062*25.4/2.*mm;
  G4double ringOuterRad = 2.045*25.4/2.*mm;
  G4double ringInnerRad = ringOuterRad-2*ringHThick;
  G4double ring_posz = THolder2_posz-THolder2_Hthick-ringHHeight;

  G4Tubs* Tring_solid
    = new G4Tubs("Target Holder Raised Ring",ringInnerRad,ringOuterRad,ringHHeight,
		 0,twopi);
  G4LogicalVolume* Tring_log = 
    new G4LogicalVolume(Tring_solid,ssteel,"Raised Ring log",0,0,0);
  // G4PhysicalVolume* Tring_phys =
  new G4PVPlacement(0,
		    G4ThreeVector(0,0,ring_posz),
		    Tring_log,
		    "Raised Ring",
		    LENALab_log,
		    false,
		    0);

  // The frame to hold down the target
  // Made by subtracting a cone from a ring
  G4double frame_outer = 47./2.*mm;
  G4double frame_HThick = 1.*mm;
  G4double frame_posz = -frame_HThick-0.01*mm;
  G4double framecone_largerad = 34.5/2.*mm;
  G4double framecone_smallrad = 31.8/2.*mm;

  // main frame
  G4Tubs* mainframe_solid
    = new G4Tubs("Main part of frame",0.*mm,frame_outer,
		 frame_HThick,0,twopi);
  // Cone
  G4Cons* subframe_solid
    = new G4Cons("Subtracted part",
		 0,                      // inner rad at -Hlength
		 framecone_largerad,     // outer rad at -Hlength
		 0,                      // inner rad at HLength
		 framecone_smallrad,     // outer rad at HLength
		 frame_HThick+0.1*mm,    // half length
		 0,twopi);               // start and end angles
  // The subtraction
  G4SubtractionSolid* frame_solid
    = new G4SubtractionSolid("Frame to hold target",mainframe_solid,subframe_solid,
			     &noRotate,G4ThreeVector(0,0,0));

  // The logical Volume
  G4LogicalVolume* frame_log =
    new G4LogicalVolume(frame_solid,ssteel,"Frame to hold target log",0,0,0);
  // G4VPhysicalVolume* frame_phys = 
  new G4PVPlacement(0, 
		    G4ThreeVector(0.,
				  0.,
				  frame_posz),
		    frame_log,
		    "Target frame",
		    LENALab_log,
		    false,
		    0);


    // water in target holder
  G4double TWater_outer = 1.160*25.4/2.*mm;
  G4double TWater_Hdepth = 0.3*25.4/2*mm;
  G4double TWater_pos = TWater_Hdepth+2.*TantBack_Hthick;
  
  G4Tubs* TWater_solid
    = new G4Tubs("Target Holder Water", 0.*mm,TWater_outer,TWater_Hdepth,
  		 0,twopi);
  G4LogicalVolume* TWater_log = 
    new G4LogicalVolume(TWater_solid,Water,"Water Log",0,0,0);
  // G4VPhysicalVolume* TWater_phys = 
  new G4PVPlacement(0,
		    G4ThreeVector(0,0,TWater_pos),
		    TWater_log,
		    "Water in Target Holder",
		    LENALab_log,
		    false,
		    0);

  //------------------------------------------------------------
  //-------------- Visualization Stuff -----------------------
  //------------------------------------------------------------

  // SSteels are Dark Grey
  G4VisAttributes* SteelVisAtt1 = new G4VisAttributes(G4Colour(0.4,0.4,0.4));
  G4VisAttributes* SteelVisAtt2 = new G4VisAttributes(G4Colour(0.35,0.35,0.35));
  //G4VisAttributes* SteelVisAtt1 = new G4VisAttributes(G4VisAttributes::Invisible);
  //G4VisAttributes* SteelVisAtt2 = new G4VisAttributes(G4VisAttributes::Invisible);
  BeamPipe_log->SetVisAttributes(SteelVisAtt1);
  THolder_log->SetVisAttributes(SteelVisAtt1);
  THolder2_log->SetVisAttributes(SteelVisAtt2);
  THolder3_log->SetVisAttributes(SteelVisAtt2);
  Tring_log->SetVisAttributes(SteelVisAtt2);
  frame_log->SetVisAttributes(SteelVisAtt2);

  // Copper is Copper Color
  G4VisAttributes* CuVisAtt = new G4VisAttributes(G4Colour(0.84,0.49,0.1));
  BeamPipeCuMain_log->SetVisAttributes(CuVisAtt);
  BeamPipeCuEnd_log->SetVisAttributes(CuVisAtt);
  SuppressorCu_log->SetVisAttributes(CuVisAtt);

  // Ceramic is White
  G4VisAttributes* ceramicVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  SuppressorCeramic_log->SetVisAttributes(ceramicVisAtt);

  // Tantalum backing is darker
  G4VisAttributes* TaVisAtt = new G4VisAttributes(G4Colour(0.2,0.2,0.2));
  TantBack_log->SetVisAttributes(TaVisAtt);

  // Water is Blue
  G4VisAttributes* WaterVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  TWater_log->SetVisAttributes(WaterVisAtt);

}

//#####################################################################
// The Collimator
void DetectorConstruction::BuildCollimator(){

  //---------------- The Collimator --------------------------

  // The Aluminum Collimator Pipe
  G4double CollimatorPipe_outer = 4.5*25.4/2.*mm;
  G4double CollimatorPipe_inner = 4.25*25.4/2.*mm;
  G4double CollimatorPipe_Hlen  = 570.0/2.*mm;
  G4RotationMatrix noRotate;

  G4Tubs* CollimatorPipe_solid
    = new G4Tubs("Collimator Pipe",CollimatorPipe_inner,CollimatorPipe_outer,
		 CollimatorPipe_Hlen,0,twopi);
  G4LogicalVolume* CollimatorPipe_log =
    new G4LogicalVolume(CollimatorPipe_solid, Al,"Collimator Pipe log",0,0,0);
  // G4VPhysicalVolume* CollimatorPipe_phys = 
    new G4PVPlacement(0,
		      G4ThreeVector(0,0,0),
		      CollimatorPipe_log,
		      "Collimator Pipe",
		      LENALab_log,
		      false,
		      0);

  // The Delrin Container
  G4double DelrinContainer_outer = 4.25*25.4/2.*mm;
  G4double DelrinContainer_inner = 4.0*25.4/2.*mm;
  G4double DelrinContainer_Hlen  = 4.14*25.4/2.*mm;

  G4Tubs* DelrinContainer_solid
    = new G4Tubs("Delrin Container",DelrinContainer_inner,DelrinContainer_outer,
		 DelrinContainer_Hlen,0,twopi);
  G4LogicalVolume* DelrinContainer_log =
    new G4LogicalVolume(DelrinContainer_solid, delrin,"Delrin Container log",0,0,0);
  // G4VPhysicalVolume* DelrinContainer_phys = 
    new G4PVPlacement(0,
		      G4ThreeVector(0,0,0),
		      DelrinContainer_log,
		      "Delrin Container",
		      LENALab_log,
		      false,
		      0);

  // Lead Cylinder
  G4double LeadRad = 4.0*25.4/2.*mm;
  G4double LeadThick  = 2.0*25.4/2.*mm;

  //G4ThreeVector lead_cylL(0,0,-31.75*mm);
  //G4ThreeVector lead_cylR(0,0,25.401*mm);
  G4ThreeVector lead_cylL(0,0,-27.178*mm);
  G4ThreeVector lead_cylR(0,0,27.178*mm);


  G4Tubs* lead_cyl
    = new G4Tubs("LeadCylinder",0.*mm,LeadRad,LeadThick,0,twopi);

  G4LogicalVolume* leadcyl_log
    = new G4LogicalVolume(lead_cyl, Pb,"Lead Cylinder log",0,0,0);

  new G4PVPlacement(0,                             		// no rotation
                        lead_cylL,				// at (0,0,-31.75*mm)
                        leadcyl_log,          	   		// Logical Volume
                        "leadcylL",            			// Name
                        LENALab_log,	           		// Mother Volume
                        false,                     		// no boolean operations
                        0);                        		// Copy number

  new G4PVPlacement(0,                             		// no rotation
                        lead_cylR,				// at (0,0,25.4*mm)
                        leadcyl_log,          	   		// Logical Volume
                        "leadcylR",            			// Name
                        LENALab_log,	           		// Mother Volume
                        false,                     		// no boolean operations
                        1);                        		// Copy number
/*
  // Gamma Source
  G4double SourceRad = 25.4/2.*mm;
  G4double SourceThick  = 0.14*25.4/2.*mm;
  G4double SourceHoleRad = 2.5*mm;
  G4double SourceHoleHLen = 3.175/2.*mm;
  G4ThreeVector holepos(0,0,0);


  G4Tubs* sourceMain_solid 
    = new G4Tubs("Plastic Source",0.*mm,SourceRad,SourceThick,0,twopi);
  G4Tubs* sourceHole_solid
    = new G4Tubs("Source hole",0.*mm,SourceHoleRad,SourceHoleHLen,0,twopi);

  G4SubtractionSolid* source_solid
    =new G4SubtractionSolid("Plastic Source",sourceMain_solid,
			    sourceHole_solid,&noRotate,holepos);

  G4LogicalVolume* source_log
    = new G4LogicalVolume(source_solid, plastic,"source_log",0,0,0);

    new G4PVPlacement(0,
		      G4ThreeVector(0,0,0),
		      source_log,
		      "Gamma Source",
		      LENALab_log,
		      false,
		      0);

  // Source is Red
  G4VisAttributes* PlasticVisAtt = new G4VisAttributes(G4Colour(1.0,0,0));
  source_log->SetVisAttributes(PlasticVisAtt);
*/

  //------------------------------------------------------------
  //-------------- Visualization Stuff -----------------------
  //------------------------------------------------------------

  // Aluminum is gray
  G4VisAttributes* AlVisAtt = new G4VisAttributes(G4Colour(0.8,0.8,0.8));
  CollimatorPipe_log->SetVisAttributes(AlVisAtt);

  // Delrin is White
  G4VisAttributes* delrinVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  DelrinContainer_log->SetVisAttributes(delrinVisAtt);

  // Lead is Black
  G4VisAttributes* LeadVisAtt = new G4VisAttributes(G4Colour(0.1,0.1,0.1));
  leadcyl_log->SetVisAttributes(LeadVisAtt);

}
//#####################################################################
// The Radioactive Source
void DetectorConstruction::BuildRadioactiveSource(){

  // ---------------- Radioactive Source ------------
  G4double sourceRad = 25.4/2.*mm;
  G4double sourceHLen = 6.35/2.*mm;
  G4double sourceHoleRad = 5.0/2.*mm;
  G4double sourceHoleHLen = 3.175/2.*mm;
  G4double source_posz = -sourceHLen;
  G4double sourceHole_posz = 2.77*mm + sourceHoleHLen - sourceHLen;  // with respect to plastic
  G4ThreeVector holepos(0,0,sourceHole_posz);
  G4RotationMatrix noRotate;

  G4Tubs* sourceMain_solid 
    = new G4Tubs("Plastic Source",0.*mm,sourceRad,sourceHLen,0,twopi);
  G4Tubs* sourceHole_solid
    = new G4Tubs("Source hole",0.*mm,sourceHoleRad,sourceHoleHLen,0,twopi);
  G4SubtractionSolid* source_solid
    = new G4SubtractionSolid("Plastic Source",sourceMain_solid,
			    sourceHole_solid,&noRotate,holepos);

  G4LogicalVolume* source_log =
    new G4LogicalVolume(source_solid, plastic, "source_log",0,0,0);
  //  G4VPhysicalVolume* source_phys =
  new G4PVPlacement(0,
		    G4ThreeVector(0,0,source_posz),
		    source_log,
		    "Source",
		    LENALab_log,
		    false,
		    0);

  //------------------------------------------------------------
  //-------------- Visualization Stuff -----------------------
  //------------------------------------------------------------

  // Source is Red
  G4VisAttributes* PlasticVisAtt = new G4VisAttributes(G4Colour(1.0,0,0));
  source_log->SetVisAttributes(PlasticVisAtt);
}


//####################################################################
//find the name of the material at the coordinates (x,y,z)
G4String DetectorConstruction::GetVolumeMaterial(G4double x, G4double y, G4double z){
  //G4cout << "(x,y,z) " << "(" << x << "," << y << "," << z << ")" << G4endl;
  G4ThreeVector pos(x,y,z);
  aNavigator->SetWorldVolume(LENALab_phys);
  aNavigator->LocateGlobalPointAndSetup(pos);
  
  G4TouchableHistoryHandle aTouchable = aNavigator->CreateTouchableHistoryHandle();

  /*
  G4VPhysicalVolume *aVol = aTouchable->GetHistory()->GetVolume(0);
  G4Material *aMat = aVol->GetLogicalVolume()->GetMaterial();
  G4String theMat = aMat->GetName();
  */

  G4int depth =  aTouchable->GetHistory()->GetDepth();
  G4String theMat = aTouchable->GetHistory()->GetVolume(depth)->GetLogicalVolume()->GetMaterial()->GetName();

  //G4String theMat = "test";

  return theMat;

}
