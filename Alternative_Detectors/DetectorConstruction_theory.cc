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
#include "GermSD.hh"
#include "NaIannSD.hh"
#include "ScintSD.hh"
#include "APEXSD.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"

#include "globals.hh"

#include "G4Navigator.hh"

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

using namespace CLHEP;

DetectorConstruction::DetectorConstruction()
  : ActGeCrys_log(0),GeCrys_log(0),NaISeg_log(0),APEXNaI_log(0)
{
  //Default for printing materials table
  printMaterialTable.assign("no");

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

  elecEkinMin = 0*eV;

  //needed in multiple funcitons
  LeadHthick = 25.4/4.*mm;
  HlenLeadTopx = 457.2/2.+2*LeadHthick;
  xposLeadSide  = HlenLeadTopx - LeadHthick;
  
  aNavigator= new G4Navigator();
 
}

DetectorConstruction::~DetectorConstruction()
{
  delete aNavigator;

  delete detectorMessenger;
 }



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

  //C = 
  //new G4Material("Carbon",z=6.,a=12.011*g/mole,density=2.267*g/cm3);//graphite near room temp.

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
  G4Element* elCo = new G4Element( "Cobalt", "Co", 27. ,a);
  a = 55.85*g/mole;
  G4Element* elFe = new G4Element("Iron" , "Fe" , z= 26., a);
  a = 63.546*g/mole;
  G4Element* elCu = new G4Element("Copper","Cu",z=29., a);
  a = 65.39*g/mole;
  G4Element* elZn = new G4Element("Zinc","Zn",z=30., a);
  a = 118.710*g/mole;
  G4Element* elSn = new G4Element("Tin","Sn",z=50.,a);
  a = 26.98*g/mole;
  G4Element* elAl = new G4Element("Aluminium","Al",z=13.,a);
  a = 58.69*g/mole;
  G4Element *elNi = new G4Element("Nickel","Ni",z=28,a);
  a = 95.96*g/mole;
  G4Element *elMo = new G4Element("Molybdenum","Mo",z=42,a);
  a = 40.078*g/mole;
  G4Element *elCa = new G4Element("Calcium","Ca",z=20,a);
  a = 24.3050*g/mole;
  G4Element *elMg = new G4Element("Magnesium","Mg",z=12,a);
  a = 28.0855*g/mole;
  G4Element *elSi = new G4Element("Silicon","Si",z=14,a);
  a = 30.9738*g/mole;
  G4Element *elP = new G4Element("Phosphorus","P",z=15,a);
  a = 32.065*g/mole;
  G4Element *elS = new G4Element("Sulfur","S",z=16,a);
  a = 54.938*g/mole;
  G4Element *elMn = new G4Element("Manganese","Mn",z=25,a);

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

  //Steel1010A
  steel1010A = new G4Material
    ("Steel1010A", density=7.87*g/cm3, ncomponents=5);
  steel1010A->AddElement(elFe,0.9918);
  steel1010A->AddElement(elMn,0.0030);
  steel1010A->AddElement(elS,0.0005);
  steel1010A->AddElement(elP,0.0004);
  steel1010A->AddElement(elC,0.0008);
  

  //mylar (alenia spazio)
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

  // Aluminium Oxide ceramic
  density = 4.0 *g/cm3;
  ceramic = new G4Material("Aluminium Oxide Ceramic",density,ncomponents=2);
  ceramic -> AddElement(elAl,2./5.);
  ceramic -> AddElement(elO,3./5.);

  // Mu - Metal
  density = 8747*kg/m3;
  mumetal = new G4Material("Mu-Metal",density,ncomponents=4);
  mumetal->AddElement(elNi,0.75);
  mumetal->AddElement(elFe,0.15);
  mumetal->AddElement(elCu,0.05);
  mumetal->AddElement(elMo,0.05);

  //sytrofoam (Expanded polystyrene)
  density = 0.01627*g/cm3;
  styrofoam = new G4Material("styrofoam",density,ncomponents=2);
  styrofoam->AddElement(elH,0.077418);
  styrofoam->AddElement(elC,0.922582);


  //Aluminum Oxide (Al2O3)
  density = 3.95*g/cm3;
  Al2O3 = new G4Material("aluminum oxide",density,ncomponents=2);
  Al2O3->AddElement(elAl,0.53);
  Al2O3->AddElement(elO,0.47);

  //Calcium Oxide (CaO)
  density = 3.35*g/cm3;
  CaO = new G4Material("calcium oxide",density,ncomponents=2);
  CaO->AddElement(elCa,0.71);
  CaO->AddElement(elO,0.29);

  //Iron(II) Oxide (FeO)
  density = 5.745*g/cm3;
  FeO = new G4Material("iron(II) oxide",density,ncomponents=2);
  FeO->AddElement(elFe,0.78);
  FeO->AddElement(elO,0.22);

  //Magnesium Oxide (MgO)
  density = 3.58*g/cm3;
  MgO = new G4Material("magnesium oxide",density,ncomponents=2);
  MgO->AddElement(elMg,0.60);
  MgO->AddElement(elO,0.40);

  //Phosphorus Pentoxide (P2O5)
  //O is the letter, not the number
  density = 2.39*g/cm3;
  P2O5 = new G4Material("phosphorus pentoxide",density,ncomponents=2);
  P2O5->AddElement(elP,0.22);
  P2O5->AddElement(elO,0.78);

  //Silicon Dioxide (SiO2)
  density = 2.648*g/cm3;
  SiO2 = new G4Material("silicon dioxide",density,ncomponents=2);
  SiO2->AddElement(elSi,0.47);
  SiO2->AddElement(elO,0.53);

  //Sodium Oxide (Na2O)
  density = 2.27*g/cm3;
  Na2O = new G4Material("sodium oxide",density,ncomponents=2);
  Na2O->AddElement(elNa,0.74);
  Na2O->AddElement(elO,0.26);

  //Nickel(II) Oxide (NiO)
  density = 6.67*g/cm3;
  NiO = new G4Material("nickel (II) oxide",density,ncomponents=2);
  NiO->AddElement(elNi,0.79);
  NiO->AddElement(elO,0.21);

  // Concrete
  G4NistManager* NISTman = G4NistManager::Instance();
  concrete  = NISTman->FindOrBuildMaterial("G4_CONCRETE");

  // Quartz
  Quartz = new G4Material("Quartz", density= 2.66*g/cm3,
                                     ncomponents=2);
  Quartz->AddElement(elSi, 1);
  Quartz->AddElement(elO, 2);

  // Delrin
  delrin = NISTman->FindOrBuildMaterial("G4_POLYOXYMETHYLENE");

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

  //meteor rock
/*
  density = 3.73*g/cm3;
  meteorite = new G4Material("meteorite",density,ncomponents=8);
  
    //if number provided are mass fraction
  meteorite->AddMaterial(Al2O3,0.08);
  meteorite->AddMaterial(CaO,0.06);
  meteorite->AddMaterial(FeO,0.15);
  meteorite->AddMaterial(MgO,0.24);
  meteorite->AddMaterial(Na2O,0.05);
  meteorite->AddMaterial(NiO,0.01);
  //meteorite->AddMaterial(elP2O5,??);
  meteorite->AddElement(elS,0.01);
  meteorite->AddMaterial(SiO2,0.40);
  */
  /*
  //if number provided are by number (but entered in mass fraction, as Geant4 likes)
  meteorite->AddMaterial(Al2O3,0.136);
  meteorite->AddMaterial(CaO,0.056);
  meteorite->AddMaterial(FeO,0.179);
  meteorite->AddMaterial(MgO,0.161);
  meteorite->AddMaterial(Na2O,0.052);
  meteorite->AddMaterial(NiO,0.012);
  //meteorite->AddMaterial(elP2O5,??);
  meteorite->AddElement(elS,0.005);
  meteorite->AddMaterial(SiO2,0.399);
  */

  //Make meteorite air, for source runs
  density = 1.290*mg/cm3;
  meteorite = new G4Material("meteorite",density,ncomponents=2);
  meteorite->AddElement(elN, fractionmass=0.7);
  meteorite->AddElement(elO, fractionmass=0.3);
  
  
  /*
    //make meteorite material vacuum to calculate self absorption of meteor rock 
    // run as real meteor rock, then as vacuum, then compare. Difference is self absorption
    meteorite =
    new G4Material("meteorite", atomicNumber,
    massOfMole, density, kStateGas,
    temperature, pressure);
  */
  


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

  //Print the defined Materials if desired by user
  if(printMaterialTable.compare("yes")==0){
    //  G4cout << G4endl << "The Defined Materials are : " << G4endl << G4endl;
    //  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  }
}


G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // First clean out the old geom
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  //  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  //****** VISUALIZATION ********//
    /*  if(pVVisManager)
    {
      G4Polyline x_axis;
      //Set red line color                                                                              
      G4Colour red(1.0,0.0,0.0);
      G4VisAttributes att(red);
      att.SetColour(red);
      x_axis.push_back( G4Point3D (0.,0.,0.));
      x_axis.push_back( G4Point3D (0.,500.*cm,0.));
      x_axis.SetVisAttributes(att);
      pVVisManager->Draw(x_axis);
    }*/
//**************************************************************************

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
  // ----------WORLD-------------

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

  ////////////////////////////////////////////////////////////////////
  //
  // Build the Apparatus
  //
  //////////////////////////////////////////////////////////////////////


  // ------------ The Germanium Detector ----------------
  // To turn this off, you must also turn off HPGe sensitive detectors -- this is a pain!
  //BuildHPGe();
  //  G4cout << "Built HPGe yo" << G4endl;

  //  BuildAero();
  //    BuildMagnet();
  // --------------   Build the beam pipe -------------
  // The beam pipe, can be replaced by the NaI(Tl) plug
  //  BuildBeamPipe();  //remove for meteorite
  //  G4cout << "Built Beam Pipe" << G4endl;

  // ---------------- The Target Holder -------------
  //  BuildTargetHolder();  //remove for meteorite
  //  G4cout << "Built Target Holder" << G4endl;

  //------------------ NaI Annulus ------------------
  // To turn this off, you must also turn off NaI sensitive detectors
  //BuildNaIann(); 
  //G4cout << "Built NaI Annulus" << G4endl;
  //BuildPbShield(); 
  //G4cout << "Built Lead Shielding" << G4endl;
  //BuildMuonPanels(); 
  //G4cout << "Built Plastic Scintillators" << G4endl;

  // ---------------- Radioactive Source ------------
  // Plasic (mylar) to hold the 'puck' source -- not to be included in 'beam target' simulations
  BuildSource(); 
  //  G4cout << "Built Source" << G4endl;

  // ---------------- The Meteorite --------------------
  // The rock plus the styrofoam holder
  //BuildMeteorite();
  //G4cout << "Built Metorite Source" << G4endl;

  // ---------------- The Walls --------------------
  // Concrete walls are important for the background spectrum, not source/target
  //BuildWalls();
  //G4cout << "Built Walls" << G4endl;

  // ---------------- The NaI Plug --------------------
  // A NaI(Tl) detector that replaces the beam pipe
  //    To use this, remember to set the NaIPlug sensitive detector (SD) below
  //  BuildNaIPlug();
  //  G4cout << "Built NaI(Tl) Plug" << G4endl;

  // ---------------- The APEX Detector  --------------------
  // A NaI(Tl) detector that replaces the Annulus
  //BuildAPEX();
  //  G4cout << "Built APEX" << G4endl;

  // ---------------- The APEX Cradle  --------------------
  // A lead and aluminum cradle that surrounds the APEX detector
  //BuildAPEXcradle();
  //G4cout << "Built APEX Cradle" << G4endl;

  // ---------------- The APEX Collimator  --------------------
  // A NaI(Tl) detector that replaces the Annulus
  //BuildCollimator();
  //G4cout << "Built Collimator" << G4endl;

  // ---------------- The Gas Target (for DIANA) ---------
  // The Gas Target geometry to be used for the underground accerator
  //BuildGasTarget();
  //G4cout << " Built Gas Target" << G4endl;

  ////////////////////////////////////////////////////////////////////////




  //------------------------------------------------------------
  //-------------- Sensitive Detectors -----------------
  //-----------------------------------------------------------
  
  // Sensitive Ge Detector
  /*  if(!SDflag){
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    G4String SDname;
    
    GeSD = new GermSD(SDname="/GeSD", this);
    SDman->AddNewDetector(GeSD);
  
    NaISD = new NaIannSD(SDname="/NaISD", this);
    SDman->AddNewDetector(NaISD);

    ApexSD = new APEXSD(SDname="/APEXSD", this);
    SDman->AddNewDetector(ApexSD);

    SciSD = new ScintSD(SDname="/ScintSD", this);
    SDman->AddNewDetector(SciSD);

    SDflag=true;
  }
  //comment out unused volume_log:
  ActGeCrys_log->SetSensitiveDetector(GeSD);
  //NaISeg_log->SetSensitiveDetector(NaISD);
  //NaIPlug_log->SetSensitiveDetector(NaISD);
  APEXNaI_log->SetSensitiveDetector(ApexSD);

  //4 volumes for muon veto box
  //SciBoxTop_log->SetSensitiveDetector(SciSD);
  //SciBoxSide_log->SetSensitiveDetector(SciSD);
  //SciBoxBeam_log->SetSensitiveDetector(SciSD);
  //SciBoxDet_log->SetSensitiveDetector(SciSD);
  */
  // The World is invisible
  LENALab_log ->SetVisAttributes(G4VisAttributes::Invisible);
  

  //-----------------------------------------------------------------
  return LENALab_phys;
}



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

void DetectorConstruction::BuildNaIann(){

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
  NaISeg_log = 
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
  //G4double AlNaIInnerRad = NaIInnerRad-1.0*mm;
  //G4double AlNaIOuterRad = NaIOuterRad+1.0*mm;
  //G4double AlNaISegHLength = NaISegHLength+1.0*mm;
  
  // make one segment first
  //G4double AlNaISegStartAngle = NaISegStartAngle;
  //G4double AlNaISegSpanningAngle = NaISegSpanningAngle+1.0*deg;

  //G4Tubs* AlNaISegment_solid 
  //= new G4Tubs("Al NaI Segment",AlNaIInnerRad,AlNaIOuterRad,
  //	 AlNaISegHLength,AlNaISegStartAngle,AlNaISegSpanningAngle);
  
  // G4SubtractionSolid* OpticalIsolationNaI_solid
  //  = new G4SubtractionSolid("NaI Optical Isolation",AlNaISegment_solid,
  //			     NaISegment_solid,&noRotate,G4ThreeVector(0,0,0));
  //G4LogicalVolume* OpticalIsolationNaI_log = 
  //  new G4LogicalVolume(OpticalIsolationNaI_solid, Al,"Optical Isolation",
  //			0,0,0);
  
  
  //   for(G4int i=0;i<8;i++){
  //     NaIsegRotate=(i*45.+22.5)*deg;
  //     new G4PVPlacement(G4Transform3D(G4RotationMatrix(zAxis,NaIsegRotate),
  // 				    G4ThreeVector(0,0,-NaISegHLength-3.2)),
  // 		      OpticalIsolationNaI_log,
  // 		      "NaI Optical Isolation Segments 1",
  // 		      LENALab_log,
  // 		      false,
  // 		      0);
  //   }
  
  //   for(G4int i=0;i<8;i++){
  //     NaIsegRotate=(i*45.+22.5)*deg;
  //     new G4PVPlacement(G4Transform3D(G4RotationMatrix(zAxis,NaIsegRotate),
  // 				    G4ThreeVector(0,0,NaISegHLength+3.2)),
  // 		      OpticalIsolationNaI_log,
  // 		      "NaI Optical Isolation Segments 2",
  // 		      LENALab_log,
  // 		      false,
  // 		      0);
  //   }
  
  //PMTs -- make one, then place 16
  G4double PMTInnerRad = 30.*mm;
  G4double PMTOuterRad = 31.*mm;
  G4double PMTHLength = 167.75/2.*mm;
  G4Tubs *PMTOutertube_solid
    = new G4Tubs("PMT Outer Tube",0.*mm,PMTOuterRad,PMTHLength,0.,twopi);
  
  G4double PMTInnerHLength = PMTHLength-(PMTOuterRad-PMTInnerRad);
  G4Tubs *PMTvolume_solid
    =new G4Tubs("PMT Volume",0.*mm,PMTInnerRad,PMTInnerHLength,0,twopi);

  G4SubtractionSolid *PMTTube_solid
    = new G4SubtractionSolid("PMT Tube sub solid",PMTOutertube_solid,PMTvolume_solid,
			     &noRotate,G4ThreeVector(0,0,0));
  
  G4LogicalVolume *PMTtube_log
    = new G4LogicalVolume(PMTTube_solid,mumetal,"PMT tube",0,0,0);
  
  
  G4double PMTsegRotate;
  G4double PMTx;// = (NaIOuterRad-PMTOuterRad)*cos(NaISegStartAngle);
  G4double PMTy;// = (NaIOuterRad-PMTOuterRad)*sin(NaISegStartAngle);

  for(int i=0;i<8;i++){
    //PMTsegRotate=(i*45.+22.5+0.5)*deg;
    PMTsegRotate=(i*45.)*deg;
    PMTx = (NaIOuterRad-PMTOuterRad)*cos(NaISegStartAngle+PMTsegRotate);
    PMTy = (NaIOuterRad-PMTOuterRad)*sin(NaISegStartAngle+PMTsegRotate);
    new G4PVPlacement(G4Transform3D(G4RotationMatrix(zAxis,PMTsegRotate),
				    G4ThreeVector(PMTx,PMTy,PMTHLength+2.*NaISegHLength+3.2)),
		      PMTtube_log,
		      "PMTs",
		      LENALab_log,
		      false,
		      i);
  }
  for(int i=0;i<8;i++){
    //PMTsegRotate=(i*45.+22.5+0.5)*deg;
    PMTsegRotate=(i*45.)*deg;
    PMTx = (NaIOuterRad-PMTOuterRad)*cos(NaISegStartAngle+PMTsegRotate);
    PMTy = (NaIOuterRad-PMTOuterRad)*sin(NaISegStartAngle+PMTsegRotate);
    new G4PVPlacement(G4Transform3D(G4RotationMatrix(zAxis,PMTsegRotate),
				    G4ThreeVector(PMTx,PMTy,-PMTHLength-2.*NaISegHLength-3.2)),
		      PMTtube_log,
		      "PMTs",
		      LENALab_log,
		      false,
		      i+8);
  }
  PMTtube_log->SetVisAttributes(G4VisAttributes(G4Colour(0.8,0.8,0.8)));

  
  // Make any aluminium off White/grey
  G4VisAttributes* AlVisAtt = new G4VisAttributes(G4Colour(0.8,0.8,0.8));
  NaIHouseInner_log->SetVisAttributes(AlVisAtt);
  NaIHouseOuter_log->SetVisAttributes(AlVisAtt);
  NaIFlange_log->SetVisAttributes(AlVisAtt);
  NaIFill_log->SetVisAttributes(AlVisAtt);

  // Mylar parts are white
  //  G4VisAttributes* mylarVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));

  // The Detectors are green
  //G4VisAttributes* DetVisAtt = new G4VisAttributes(G4VisAttributes::Invisible);
  G4VisAttributes* DetVisAtt = new G4VisAttributes(G4Colour(0.0,0.75,0.0));
  NaISeg_log->SetVisAttributes(DetVisAtt);

}


void DetectorConstruction::BuildPbShield(){  

  // ------------ The Lead Shield ------------
  G4double LeadBoxHeight = 539.75 + (2.*LeadHthick);
  // Top
  G4double HlenLeadTopz = (930./2.)-2*25.4;
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
  G4double BeampipeRad = 118./2.*mm;//95.25/2.*mm;
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

   
  G4VisAttributes* AlBoxVisAtt = new G4VisAttributes(G4Colour(0.8,0.8,0.8));
  //AlBoxVisAtt->SetForceWireframe(true);
  AlBoxTop_log->SetVisAttributes(AlBoxVisAtt);
  AlBoxSide_log->SetVisAttributes(AlBoxVisAtt);
  AlBoxBeam_log->SetVisAttributes(AlBoxVisAtt);
  //AlBoxBeam_log->SetVisAttributes(G4VisAttributes::Invisible);
  AlBoxDet_log->SetVisAttributes(AlBoxVisAtt);
  //AlBoxDet_log->SetVisAttributes(G4VisAttributes::Invisible);


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

void DetectorConstruction::BuildMuonPanels(){

  // ------------- Scintillator Shield ---------
  G4double SciHthick = 25.4*mm;
  G4double SciBoxHeight = 558.8 + (2.*SciHthick);

  // Top
  G4double HlenSciTopz = 930./2.*mm;
  G4double HlenSciTopx = 576/2.*mm;
  G4double heightSciTop= 264.922+SciHthick*mm;

  G4RotationMatrix noRotate;
  G4ThreeVector xAxis(1,0,0);
  G4ThreeVector yAxis(0,1,0);
  G4ThreeVector zAxis(0,0,1);

  G4Box* SciBoxTop_solid
    = new G4Box("Sci Box Top", HlenSciTopx,SciHthick,HlenSciTopz);
  SciBoxTop_log = 
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
  SciBoxSide_log = 
    new G4LogicalVolume(SciBoxSide_solid, Sci,"Sci Box Side log", 0,0,0);

  G4double panelsideRotateZ1 = 0*deg;
  G4double panelsideRotateZ2 = 0*deg;
  G4double panelRotateR = (HlenSciSidey+SciHthick);
  G4double SciBoxSideX1 = xposSciSide+3.+panelRotateR*sin(panelsideRotateZ1);
  G4double SciBoxSideY1 = -yposSciSide+panelRotateR*sin(panelsideRotateZ1);
  G4double SciBoxSideZ1 = zposSciSide;
  // G4VPhysicalVolume* SciBoxSide1_phys = 
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(zAxis, panelsideRotateZ1),
				  G4ThreeVector(SciBoxSideX1,SciBoxSideY1,SciBoxSideZ1)),
		    SciBoxSide_log,
		    "Sci Box Side 1",
		    LENALab_log,
		    false,
		    1);
  
  G4double SciBoxSideX2 = -xposSciSide-3.-panelRotateR*sin(panelsideRotateZ2);
  G4double SciBoxSideY2 = -yposSciSide+panelRotateR*sin(panelsideRotateZ2);
  G4double SciBoxSideZ2 = -zposSciSide;
  // G4VPhysicalVolume* SciBoxSide2_phys = 
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(zAxis, panelsideRotateZ2),
				  G4ThreeVector(SciBoxSideX2,SciBoxSideY2,SciBoxSideZ2)),
		    SciBoxSide_log,
		    "Sci Box Side 2",
		    LENALab_log,
		    false,
		    2);
		    
  // The Detector and Beam Plate 
  G4double HlenSciBeamx = HlenSciTopx-SciHthick;
  G4double HlenSciBeamy = HlenSciSidey;
  G4double HlenSciBeamz = SciHthick;
  G4double SciBeamPlatePosz = HlenSciTopz - SciHthick;
  G4double SciBeampipeRad = 121./2.*mm;//98.425/2.*mm;
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

  SciBoxBeam_log = 
    new G4LogicalVolume(SciBoxBeam_solid, Sci,"Beam Pipe Plate log",0,0,0);
  SciBoxDet_log = 
    new G4LogicalVolume(SciBoxDet_solid, Sci,"Detector Pipe Plate log",0,0,0);
  
  // G4VPhysicalVolume* SciBoxBeam_phys = 
  //This one is removed for meteorite
  new G4PVPlacement(0,
		      G4ThreeVector(SciHthick,-yposSciSide,
				    -SciBeamPlatePosz-6.),
		    SciBoxBeam_log,
		    "Sci Box Beam Side",
		    LENALab_log,
		    false,
		    3);
  
  // G4VPhysicalVolume* SciBoxDet_phys = 
  new G4PVPlacement(0,
		    G4ThreeVector(-SciHthick,-yposSciSide,
				  SciBeamPlatePosz+6.),
		    SciBoxDet_log,
		    "Sci Box Detector Side",
		    LENALab_log,
		    false,
		    4);

  //------------------------------------------------------------
  //-------------- Visualization Stuff -----------------------
  //------------------------------------------------------------

  // Scintillator Detectors are Black
  G4VisAttributes* SciVisAtt = new G4VisAttributes(G4Colour(0.1,0.1,0.1));
  //  SciVisAtt->SetForceWireframe(true);
  SciBoxTop_log->SetVisAttributes(SciVisAtt);
  SciBoxSide_log->SetVisAttributes(SciVisAtt);
  SciBoxBeam_log->SetVisAttributes(SciVisAtt);
  //SciBoxBeam_log->SetVisAttributes(G4VisAttributes::Invisible);
  SciBoxDet_log->SetVisAttributes(SciVisAtt);
  //SciBoxDet_log->SetVisAttributes(G4VisAttributes::Invisible);


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
  G4double houseHLength = 296./2.*mm;
  G4double houseBaseRadius = houseInnerRadius;

  G4double GeMotherRadius = houseOuterRadius;
  G4double GeMotherHLength = houseHLength;
  //G4double GetoTargetLength = 175.0*mm;
  //  G4double GeMotherZpos = (SourceDetectorDist+GeMotherHLength)
  //*cos(DetectorAngle);
  //G4double GeMotherXpos = -(SourceDetectorDist+GeMotherHLength)
  //*sin(DetectorAngle);
  G4Tubs* GeMotherCylinder_solid
    = new G4Tubs("Ge_Mother",0,
  		 GeMotherRadius,GeMotherHLength,0,twopi);
  G4LogicalVolume* GeMotherCylinder_log =
    new G4LogicalVolume(GeMotherCylinder_solid, vacuum, "GeMother log",
			0,0,0);

  G4RotationMatrix* GeMother_rot = new G4RotationMatrix();
  GeMother_rot->rotateY(DetectorAngle);

  /*  new G4PVPlacement(GeMother_rot,
		    G4ThreeVector(GeMotherXpos,
				  0.*mm,
				  GeMotherZpos),
		      GeMotherCylinder_log,
		      "Germanium Detector Mother Volume",
		      LENALab_log,
		      false,
		      0);*/

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

  //varibles for 'curved end of crystal
  G4double endGe_torusRad = GeCrys_outer-GeEndRad;
  G4double endGe_posz = -GeBlock_Hlen;
  G4ThreeVector MoveendGe(0,0,endGe_posz);
  G4ThreeVector MoveendGeVis(GeHorisDisp,
			     GeVertDisp,
			     endGe_posz+GeEndRad/2.+GeCrys_posz-GeMotherHLength);

  //varibles for fill piece in end of crystal
  G4double FillPiece_HLenX = 2 *mm;//add extra length to overlap with other HPGe (if they just touch, Geant4 will a seam (or skin)
  G4double FillPiece_shiftX = 2 *mm;
  G4double FillPiece_HLen = GeEndRad/2.;
  G4double FillPiece_Rad = endGe_torusRad;
  G4double FillPiece_posz = endGe_posz-GeEndRad/2.;
  G4ThreeVector MoveFillPiece(0,0,FillPiece_posz+FillPiece_shiftX);


  
  //comment this out to just use the Act Ge, which is all that's needed (i think)
  // This was here to try and add a 'dead layer'. This didn't work, as the two Ge 
  // crystals formed a 'skin' at the boundary.
  
  G4Tubs* blockGe
    = new G4Tubs("Block of Ge (no finger)", 0.*mm,GeCrys_outer,
		 GeBlock_Hlen,0,twopi);

  // Curved end of crystal 
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
  /*  new G4PVPlacement(0, 
		    G4ThreeVector(GeHorisDisp,
				  GeVertDisp,
				  GeCrys_posz+FillPiece_HLen-GeMotherHLength),
		    GeCrys_log,
		    "Ge Crystal",
		    GeMotherCylinder_log,
		    false,
		    0);*/
  
  //--------------------------------------------------------------
  // Now, just the active volume of crystal
  //G4cout << G4BestUnit(GeDeadLayerHThick,"Length");

  G4double GeActCrys_outer = GeRadius-GeDeadLayerHThick*2.;
  G4double GeActCrys_Hlen  = GeHLength-GeDeadLayerHThick*2.;
  G4double GeActEndRad     = GeEndRad;//-GeDeadLayerHThick*2.;
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
  //G4ThreeVector MoveActendGeVis(0,
  //			0,
  //			endGeAct_posz); //-GeActEndRad+GeActCrys_posz);
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
  
  /*  new G4PVPlacement(0, 
		    G4ThreeVector(0.,
				  0.,
				  GeDeadLayerHThick),
		    ActGeCrys_log,
		    "Active Ge Crystal",
		    GeCrys_log,
		    false,
		    0);*/
  /*
  new G4PVPlacement(0, 
		    G4ThreeVector(GeHorisDisp,
				  GeVertDisp,
				  GeCrys_posz+FillPiece_HLen-GeMotherHLength),
		    ActGeCrys_log,
		    "Active Ge Crystal",
		    GeMotherCylinder_log,
		    false,
		    0);
  */
  
  // This is purely for visualisation - take out for real simulations
  /*
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
  */
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
  /*    new G4PVPlacement(0,
		      G4ThreeVector(0.*mm,
				    0.*mm,
				    Source_Detectorz+houseHLength-GeMotherHLength),
		      HouseTube_log,
		      "Ge Housing 1",
		      GeMotherCylinder_log,
		      false,
		      0);*/
  // G4VPhysicalVolume* HouseBase_phys = 
  /*    new G4PVPlacement(0,
		      G4ThreeVector(0.*mm,
				    0.*mm,
				    Source_Detectorz+houseFrontThickness
				    -GeMotherHLength),
		      HouseBase_log,
		      "Ge Housing 2",
		      GeMotherCylinder_log,
		      false,
		      0);*/



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
  //  G4double mylarPosz = GeCrys_posz;

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
  /*  new G4PVPlacement(0,
		    G4ThreeVector(GeHorisDisp,
				  GeVertDisp,
				  mylarPosz-GeMotherHLength),
		    mylarTube_log,
		    "Ge Mylar 1",
		    GeMotherCylinder_log,
		    false,
		    0);*/
  // G4VPhysicalVolume* mylarBase_phys = 
  /*  new G4PVPlacement(0,
		    G4ThreeVector(GeHorisDisp,
				  GeVertDisp,
				  mylarPosz-mylarHLength+mylarHThick-GeMotherHLength),
		    mylarBase_log,
		    "Ge Mylar Base",
		    GeMotherCylinder_log,
		    false,
		    0);*/
  
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
  /*  new G4PVPlacement(0,
		    G4ThreeVector(GeHorisDisp,
				  GeVertDisp,
				  alummylarPosz-GeMotherHLength),
		    alummylarTube_log,
		    "Ge Aluminised mylar",
		    GeMotherCylinder_log,
		    false,
		    0);*/
  // G4VPhysicalVolume* alummylarBase_phys = 
  /*  new G4PVPlacement(0,
		    G4ThreeVector(GeHorisDisp,
				  GeVertDisp,
				  alummylarPosz-alummylarHLength-2.*mylarHThick+
				  alummylarHThick-GeMotherHLength),
		    alummylarBase_log,
		    "Ge aluminised Mylar Base",
		    GeMotherCylinder_log,
		    false,
		    0);*/


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
  /*    new G4PVPlacement(0,
		      G4ThreeVector(GeHorisDisp,
				    GeVertDisp,
				    almylarPosz-GeMotherHLength),
		      almylarTube_log,
		      "Ge Almylar 1",
		      GeMotherCylinder_log,
		      false,
		      0);*/

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
  /*    new G4PVPlacement(0,
		      G4ThreeVector(GeHorisDisp,
				    GeVertDisp,
				    basePosz-GeMotherHLength),
		      GeBase_log,
		      "Germanium Base Platform",
		      GeMotherCylinder_log,
		      false,
		      0);*/
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
  //  G4double boron_posz = GeCrys_posz+FillPiece_HLen+GeFinger_posz; 
  // relative to crystal

  G4Tubs* boron_solid
    = new G4Tubs("boron",boron_inner ,boron_outer,
		 boron_Hlen,0,twopi);
  G4LogicalVolume* boron_log = 
    new G4LogicalVolume(boron_solid,B,"Boron log",0,0,0);
  // G4VPhysicalVolume* boron_phys =
  /*    new G4PVPlacement(0,
		      G4ThreeVector(GeHorisDisp,
				    GeVertDisp,
				    boron_posz-GeMotherHLength),
		      boron_log,
		      "Boron in cold finger",
		      GeMotherCylinder_log,
		      false,
		      0);*/


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
      /*      new G4PVPlacement(0,
			G4ThreeVector(GeHorisDisp,
				      GeVertDisp,
				      fingerPosz-GeMotherHLength),
			CuFinger_log,
			"Cu Cold Finger",
			GeMotherCylinder_log,
			false,
			0);*/
      
      
      CuFinger_log->SetVisAttributes(G4VisAttributes(G4Colour(0.84,0.49,0.1)));

      // Now, we need to ass the naval brass finger inside the copper finger
      G4double brassFingerNarrowRad = 2.9/2.*mm;
      G4double brassFingerNarrowHLength = 51.6/2.*mm;
      G4double brassFingerWideRad = 3.8/2.*mm;
      G4double brassFingerWideHLength = 25.2/2.*mm;
      G4double brassFingerExtension = fingerOffset-0.5*mm;   // trial and error
      //      G4double brassFingerNarrowPosz = fingerPosz-fingerHLength+
	brassFingerNarrowHLength+2.*brassFingerWideHLength-brassFingerExtension;
	//G4double brassFingerWidePosz = fingerPosz-fingerHLength+
	brassFingerWideHLength-brassFingerExtension;

      G4Tubs* brassFingerNarrow_solid
	= new G4Tubs("Narrow Brass Finger",0.*mm,brassFingerNarrowRad,
		     brassFingerNarrowHLength,0,twopi);
      G4LogicalVolume* brassFingerNarrow_log = 
	new G4LogicalVolume(brassFingerNarrow_solid,brass,"Narrow brass Finger",0,0,0);
      //    G4VPhysicalVolume* brassFingerNarrow_phys =
      /*      new G4PVPlacement(0,
			G4ThreeVector(GeHorisDisp,
				      GeVertDisp,
				      brassFingerNarrowPosz-GeMotherHLength),
			brassFingerNarrow_log,
			"Narrow Brass Finger",
			GeMotherCylinder_log,
			false,
			0);*/

      G4Tubs* brassFingerWide_solid
	= new G4Tubs("Wide Brass Finger",0.*mm,brassFingerWideRad,
		     brassFingerWideHLength,0,twopi);
      G4LogicalVolume* brassFingerWide_log = 
	new G4LogicalVolume(brassFingerWide_solid,brass,"Wide brass Finger",0,0,0);
      //    G4VPhysicalVolume* brassFingerWide_phys =
      /*      new G4PVPlacement(0,
			G4ThreeVector(GeHorisDisp,
				      GeVertDisp,
				      brassFingerWidePosz-GeMotherHLength),
			brassFingerWide_log,
			"Wide Brass Finger",
			GeMotherCylinder_log,
			false,
			0);*/
 
      G4VisAttributes* BrassVisAtt = new G4VisAttributes(G4Colour(0.74,0.39,0.1));
      brassFingerNarrow_log ->SetVisAttributes(BrassVisAtt);
      brassFingerWide_log ->SetVisAttributes(BrassVisAtt);

    }


    // Copper nut
    G4Tubs * CuNut_solid
      = new G4Tubs("Copper Nut",0.*mm,fingerNutRad,fingerNutHLength,0,twopi);
    G4LogicalVolume* CuNut_log =
      new G4LogicalVolume(CuNut_solid,Cu,"CuNut log",0,0,0);
    /*    new G4PVPlacement(0,
		      G4ThreeVector(GeHorisDisp,
				    GeVertDisp,
				    fingerNutPosz-GeMotherHLength),
		      CuNut_log,
		      "Cu Nut",
		      GeMotherCylinder_log,
		      false,
		      0);*/

    //Aluminum Pedestal
    G4double AlPedRad = 44.5/2*mm;
    G4double AlPedHLength = 3.2/2.*mm;

    G4Tubs *Alped_solid
      =new G4Tubs("Aluminum Pedestal",0.*mm,AlPedRad,AlPedHLength,0,twopi);
    G4LogicalVolume *Alped_log = 
      new G4LogicalVolume(Alped_solid,Al,"AlPed log",0,0,0);
    new G4PVPlacement(0,
		    G4ThreeVector(0,
				  0,
				   basePosz-GeMotherHLength+0.4*mm+3.2*mm),
		    Alped_log,
		    "Al Pedestal",
		    GeMotherCylinder_log,
		    false,
		    0);
    G4VisAttributes* alpedVisAtt = 
      new G4VisAttributes(G4Colour(0.8,0.8,0.8));
    Alped_log->SetVisAttributes(alpedVisAtt);
    

    // Cu cold finger behind the crystal
    
    G4double CuRodRad = 32./2.*mm;
    G4double CuRodHLength = houseHLength-135./2.*mm-3.7*mm;
    
    G4Tubs *CuRod_solid
      = new G4Tubs("Copper Rod",0.*mm,CuRodRad,CuRodHLength,0,twopi);
    G4LogicalVolume *CuRod_log = 
      new G4LogicalVolume(CuRod_solid,Cu,"CuRod log",0,0,0);
    /*    new G4PVPlacement(0,
		      G4ThreeVector(0,
				    0,
				    basePosz-GeMotherHLength+0.4*mm+3.2*mm+1.7*mm+CuRodHLength),
		      CuRod_log,
		      "Cu Rod",
		      GeMotherCylinder_log,
		      false,
		      0);*/
    G4VisAttributes* cuRodVisAtt = 
      new G4VisAttributes(G4Colour(0.84,0.49,0.1));
    CuRod_log->SetVisAttributes(cuRodVisAtt);
    
    // Al back plate to the HPGe Housing
    G4double AlGeBaseRad = houseInnerRadius;//houseBaseRadius-houseOuterRadius-houseInnerRadius;
    G4double AlGeBaseHLength = 7.5/2.*mm;

    G4Tubs *AlGeBase_solid
      = new G4Tubs("Al HPGe Base",CuRodRad+1.*mm,AlGeBaseRad,AlGeBaseHLength,0,twopi);
    G4LogicalVolume *AlGeBase_log = 
      new G4LogicalVolume(AlGeBase_solid,Al,"Al HPGe Base log",0,0,0);
    /*    new G4PVPlacement(0,
		      G4ThreeVector(0,
				    0,
				    Source_Detectorz-AlGeBaseHLength
				    +GeMotherHLength),
		      AlGeBase_log,
		      "Al Ge Base",
		      GeMotherCylinder_log,
		      false,
		      0);*/
    G4VisAttributes* algebaseVisAtt = 
      new G4VisAttributes(G4Colour(0.8,0.8,0.8));
    AlGeBase_log->SetVisAttributes(algebaseVisAtt);


    // Al back plate in the Ge housing
    G4double AlGeMidRad = houseInnerRadius;//houseBaseRadius-houseOuterRadius-houseInnerRadius;
    G4double AlGeMidHLength = 13.8/2.*mm;

    G4Tubs *AlGeMid_solid
      = new G4Tubs("Al HPGe Middle Plate",CuRodRad+1.*mm,AlGeMidRad,AlGeMidHLength,0,twopi);
    G4LogicalVolume *AlGeMid_log = 
      new G4LogicalVolume(AlGeMid_solid,Al,"Al HPGe Plate log",0,0,0);
    /*    new G4PVPlacement(0,
		      G4ThreeVector(0,
				    0,
				    Source_Detectorz-AlGeBaseHLength
				    +GeMotherHLength-65.2*mm),
		      AlGeMid_log,
		      "Al Ge Mid Plate",
		      GeMotherCylinder_log,
		      false,
		      0); //plate locate ~65.2 mm behind the back base of the housing*/
    G4VisAttributes* algeplateVisAtt = 
      new G4VisAttributes(G4Colour(0.8,0.8,0.8));
    AlGeMid_log->SetVisAttributes(algeplateVisAtt);

    //stainless tube to the dewer
    G4double StTubeRad = 32./2.*mm;  
    G4double StTubeHLength = 110./2.*mm;
    G4Tubs *StTube_solid
      = new G4Tubs("St Tube",0.*mm,StTubeRad,StTubeHLength,0,twopi);
    G4LogicalVolume *StTube_log
      = new G4LogicalVolume(StTube_solid,ssteel,"Steel Tube",0,0,0);
    /*    new G4PVPlacement(GeMother_rot,
		      G4ThreeVector(GeMotherXpos,
				    0,
				    GeMotherZpos+houseHLength+StTubeHLength),
		      StTube_log,
		      "Steel Tube",
		      LENALab_log,
		      false,
		      0);*/
    G4VisAttributes* sttubeVisAtt = 
      new G4VisAttributes(G4Colour(0.8,0.8,0.8));
    StTube_log->SetVisAttributes(sttubeVisAtt);

    //Cu pipe in the stainless tube to the dewer
    G4double CuTubeHLength = 110./2.*mm;
    G4double CuTubeRad = 30./2.*mm;
    G4Tubs *CuTube_solid
      = new G4Tubs("Cu Tube",0.*mm,CuTubeRad,CuTubeHLength,0,twopi);
    G4LogicalVolume *CuTube_log
      = new G4LogicalVolume(CuTube_solid,Cu,"Cu Tube",0,0,0);
    /*    new G4PVPlacement(GeMother_rot,
		      G4ThreeVector(0,
				    0,
				    0),
		      CuTube_log,
		      "Cu Tube",
		      StTube_log,
		      false,
		      0);*/
    G4VisAttributes* cutubeVisAtt = 
      new G4VisAttributes(G4Colour(0.84,0.49,0.1));
    CuTube_log->SetVisAttributes(cutubeVisAtt);


    //flange where steel tube meets dewer
    G4double StFlangeRad = 92.24/2.*mm;
    G4double StFlangeHLength = 2./2.*mm;
    G4Tubs *StFlange_solid
      = new G4Tubs("St flange",0.*mm,StFlangeRad,StFlangeHLength,0,twopi);
    G4LogicalVolume *StFlange_log
      = new G4LogicalVolume(StFlange_solid,ssteel,"Steel Flange",0,0,0);
    /*    new G4PVPlacement(GeMother_rot,
		      G4ThreeVector(GeMotherXpos,
				    0,
				    GeMotherZpos+houseHLength+2.*StTubeHLength+StFlangeHLength),
		      StFlange_log,
		      "Steel Flange",
		      LENALab_log,
		      false,
		      0);*/
    StFlange_log->SetVisAttributes(sttubeVisAtt);  

    
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
  //G4VisAttributes* GeVisAtt = new G4VisAttributes(G4Colour(0.0,1,0.0));
  //GeVisAtt->SetForceWireframe(true);
  //GeCrys_log->SetVisAttributes(GeVisAtt);
  G4VisAttributes* ActGeVisAtt = new G4VisAttributes(G4Colour(1.,1.,0.));
  //G4VisAttributes* GeVisAtt = new G4VisAttributes(G4Colour(1.,0.,0.));
  //ActGeVisAtt->SetForceWireframe(true);
  ActGeCrys_log->SetVisAttributes(ActGeVisAtt);
  GeCrys_log->SetVisAttributes(ActGeVisAtt);
  /*
  Visuals=0;
  if(Visuals){
    //endGeVIS_log->SetVisAttributes(GeVisAtt);
    ActendGeVIS_log->SetVisAttributes(G4VisAttributes::Invisible);
  }
  */
  // Make boron contact invisible
  boron_log->SetVisAttributes(G4VisAttributes::Invisible);


}

//#####################################################################
// The Target Holder
void DetectorConstruction::BuildTargetHolder(){

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
  G4RotationMatrix noRotate;
  G4ThreeVector MoveWings1(0,WingPosy,0);
  G4ThreeVector MoveWings2(0,-WingPosy,0); 
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
  // Plastic on the end of the target holder for isolation
  G4double isol_outer = THolder_outer;
  G4double isol_HThick = 0.5/2.*mm;
  G4double isol_posz = THolder_posz+THolder_Hthick+isol_HThick;

  G4Tubs* isol_solid
    = new G4Tubs("Target Chamber Isolation",0.*mm,isol_outer,isol_HThick,
		 0,twopi);
  G4LogicalVolume* isol_log = 
    new G4LogicalVolume(isol_solid,mylar,"Target Chamber Isolation log",0,0,0);
  // G4VPhysicalVolume* isol_phys = 
    new G4PVPlacement(0,
		      G4ThreeVector(0,0,isol_posz),
		      isol_log,
		      "Isolation at Back of Target Holder",
		      LENALab_log,
		      false,
		      0);

  //------------------------------------------------------------
  //-------------- Visualization Stuff -----------------------
  //------------------------------------------------------------

  // Tantalum and steel are Dark Grey
  G4VisAttributes* SteelVisAtt1 = new G4VisAttributes(G4Colour(0.4,0.4,0.4));
  G4VisAttributes* SteelVisAtt2 = new G4VisAttributes(G4Colour(0.35,0.35,0.35));
  //G4VisAttributes* SteelVisAtt1 = new G4VisAttributes(G4VisAttributes::Invisible);
  //G4VisAttributes* SteelVisAtt2 = new G4VisAttributes(G4VisAttributes::Invisible);
  THolder_log->SetVisAttributes(SteelVisAtt1);
  THolder2_log->SetVisAttributes(SteelVisAtt2);
  THolder3_log->SetVisAttributes(SteelVisAtt2);
  Tring_log->SetVisAttributes(SteelVisAtt2);
  frame_log->SetVisAttributes(SteelVisAtt2);

  // Tantalum backing is darker
  G4VisAttributes* TaVisAtt = new G4VisAttributes(G4Colour(0.2,0.2,0.2));
  TantBack_log->SetVisAttributes(TaVisAtt);

  // Water is Blue
  G4VisAttributes* WaterVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  TWater_log->SetVisAttributes(WaterVisAtt);

  // isolator is black
  G4VisAttributes* isolVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,0.0));
  isol_log->SetVisAttributes(isolVisAtt);

}


//#####################################################################
// The Beam Pipe
void DetectorConstruction::BuildBeamPipe(){

  G4RotationMatrix noRotate;

  //---------------- The Beam Pipe --------------------------
  // The Al Beam Pipe
  G4double BeamPipeAl_outer = 2.375*25.4/2.*mm;
  G4double BeamPipeAl_inner = 2.05*25.4/2.*mm;
  G4double BeamPipeAl_Hlen  = 10*25.4*mm;
  G4double BeamPipeAl_posz  = -BeamPipeAl_Hlen+0.5*mm; //+ 0.375*25.4;
  
  G4Tubs* BeamPipeAl_solid
    = new G4Tubs("Al Beam Pipe",BeamPipeAl_inner,BeamPipeAl_outer,
		 BeamPipeAl_Hlen,0,twopi);
  G4LogicalVolume* BeamPipeAl_log =
    new G4LogicalVolume(BeamPipeAl_solid, ssteel,"Al Beam Pipe log",0,0,0);
  // G4VPhysicalVolume* BeamPipeAl_phys = 
  new G4PVPlacement(0,
		    G4ThreeVector(0,0,BeamPipeAl_posz),
		    BeamPipeAl_log,
		    "Al Beam Pipe",
		    LENALab_log,
		    false,
		    0);
  // The Copper Cold Trap comes in two pieces: main and end
  G4double BeamPipeCuEnd_Hlen  = 36.5/2.*mm;
  G4double BeamPipeCuEnd_outer = 34.9/2.*mm;
  G4double BeamPipeCuEnd_inner = BeamPipeCuEnd_outer - 4.85*mm;

  G4double BeamPipeCuMain_outer = 1.375*25.4/2.*mm;
  G4double BeamPipeCuMain_inner = 1.25*25.4/2.*mm;
  G4double BeamPipeCuMain_Hlen  = BeamPipeAl_Hlen-BeamPipeCuEnd_Hlen;

  G4double SuppressorCu_Hlen  = 6.5/2.*mm;
  G4double SuppressorCu_outer = 44.76/2.*mm;
  G4double SuppressorCu_inner = 19.3/2.*mm;
  G4double SuppressorCeramic_Hlen = 14.15/2.*mm - SuppressorCu_Hlen;
  G4double SuppressorCeramic_outer = SuppressorCu_outer;
  G4double SuppressorCeramic_inner = 27.97/2.*mm;

  G4double BeamPipeCuMain_posz  = (BeamPipeAl_posz-
				   BeamPipeCuEnd_Hlen-
				   2.*(SuppressorCu_Hlen+SuppressorCeramic_Hlen)-
				   0.5*25.4*mm);
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
    /*
    G4double BeamPipeVac_Hlen = BeamPipeCuMain_Hlen+SuppressorCu_Hlen+SuppressorCeramic_Hlen+BeamPipeCuEnd_Hlen;
    G4double BeamPipeVac_posz = BeamPipeCuMain_posz;
    G4Tubs* BeamPipeVac_solid
      = new G4Tubs("Beam Pipe Vacuum Tube",0.,SuppressorCone_smallrad,
		   BeamPipeVac_Hlen,0,twopi);
    G4LogicalVolume* BeamPipeCac_log = 
      new G4LogicalVolume(BeamPipeVac_solid, vacuum,"Beam Pipe Vac log",0,0,0);
    // G4VPhysicalVolume* BeamPipeCuMain_phys =
    new G4PVPlacement(0,
		      G4ThreeVector(0,0,BeamPipeVac_posz),
		      BeamPipeCac_log,
		      "Beam Pipe Vacuum tube",
		      LENALab_log,
		      false,
		      0);
    */

  // Make any aluminium off White/grey
  G4VisAttributes* AlVisAtt = new G4VisAttributes(G4Colour(0.8,0.8,0.8));

  //AlVisAtt->SetForceWireframe(true);
  BeamPipeAl_log->SetVisAttributes(AlVisAtt);
  
  // Copper is Copper Colour
  G4VisAttributes* CuVisAtt = new G4VisAttributes(G4Colour(0.84,0.49,0.1));
  BeamPipeCuMain_log->SetVisAttributes(CuVisAtt);
  BeamPipeCuEnd_log->SetVisAttributes(CuVisAtt);
  SuppressorCu_log->SetVisAttributes(CuVisAtt);

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

G4ThreeVector DetectorConstruction::GetWorldXYZ(){
 
  G4ThreeVector MotherVolumeXYZ (LENALabx,LENALaby,LENALabz);

  return MotherVolumeXYZ;

}


void DetectorConstruction::BuildWalls(){

  G4double WallHThickness = 300.*mm;
  G4double WallHLengthX = LENALabx*0.5-2.*WallHThickness-1.*mm;
  G4double WallHLengthY = WallHThickness;
  G4double WallHLengthZ = LENALabz*0.5-2.*WallHThickness-1.*mm;

  G4Box* Wall_solid 
    = new G4Box("Walls",WallHLengthX-1*mm,WallHLengthY-1*mm,WallHLengthZ-1*mm);
  G4LogicalVolume* Wall_log = 
    new G4LogicalVolume(Wall_solid,concrete,"Wall Solid", 0,0,0);

  G4Box* Roof_solid 
    = new G4Box("Roof",WallHLengthX-1*mm,WallHLengthY*3-1*mm,WallHLengthZ-1*mm);
  G4LogicalVolume* Roof_log = 
    new G4LogicalVolume(Roof_solid,concrete,"Wall Solid", 0,0,0);

  G4Box* Floor_solid 
    = new G4Box("Floor",WallHLengthX-1*mm,WallHLengthY*3-1*mm,WallHLengthZ-1*mm);
  G4LogicalVolume* Floor_log = 
    new G4LogicalVolume(Floor_solid,concrete,"Floor Solid", 0,0,0);

  G4double WallRotate;
  G4ThreeVector xAxis(1,0,0);
  G4ThreeVector yAxis(0,1,0);
  G4ThreeVector zAxis(0,0,1);


  // The segments are now have labels 0-15. 
  // 0-7 One side
  // 8-15 Opposite side
  // Segment 0 is adjacent to seg 8, 1 next to 9 etc.
 
  //WallRotate=(i*45.+22.5+0.5)*deg;

  //along beampipe
  WallRotate=90.*deg;
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(xAxis,WallRotate),
				  G4ThreeVector(0,0,LENALabz*0.5-WallHThickness-1*mm)),
		    Wall_log,
		    "Walls",
		    LENALab_log,
		    false,
		    0);
  
  //along beampipe
  WallRotate=90.*deg;
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(xAxis,WallRotate),
				  G4ThreeVector(0,0,-LENALabz*0.5+WallHThickness+1*mm)),
		    Wall_log,
		    "Walls",
		    LENALab_log,
		    false,
		    1);
  
  //floor
  /*
  WallRotate=0.*deg;
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(xAxis,WallRotate),
				  G4ThreeVector(0,-LENALaby*0.5+(WallHThickness*3)+1.*mm,0)),
		    Floor_log,
		    "Floor",
		    LENALab_log,
		    false,
		    2);
  */
  
  G4RotationMatrix* Met_rot = new G4RotationMatrix();
  Met_rot->rotateX(0*deg);
  Met_rot->rotateY(0*deg);
  Met_rot->rotateZ(0*deg);
  new G4PVPlacement(Met_rot,
		    G4ThreeVector(0,-1100*mm-WallHThickness*3.,0),
		    Floor_log,
		    "Floor",
		    LENALab_log,
		    false,
		    2);
  
  //ceiling
  
  WallRotate=0.*deg;
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(xAxis,WallRotate),
				  G4ThreeVector(0,LENALaby*0.5-(WallHThickness*3)-1.*mm,0)),
		    Roof_log,
		    "Roof",
		    LENALab_log,
		    false,
		    3);
  /*
  new G4PVPlacement(Met_rot,
		    G4ThreeVector(0,LENALaby*0.5-(WallHThickness*1)-1.*mm,0),
		    Roof_log,
		    "Roof",
		    LENALab_log,
		    false,
		    3);
  */
  WallRotate=90.*deg;
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(zAxis,WallRotate),
				  G4ThreeVector(LENALabx*0.5-WallHThickness-1.*mm,0,0)),
		    Wall_log,
		    "Walls",
		    LENALab_log,
		    false,
		    4);

  WallRotate=90.*deg;
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(zAxis,WallRotate),
				  G4ThreeVector(-LENALabx*0.5+WallHThickness+1.*mm,0,0)),
		    Wall_log,
		    "Walls",
		    LENALab_log,
		    false,
		    5);

  Wall_log ->SetVisAttributes(G4VisAttributes::Invisible);
  Floor_log ->SetVisAttributes(G4VisAttributes::Invisible);
  Roof_log ->SetVisAttributes(G4VisAttributes::Invisible);
}

void DetectorConstruction::BuildNaIPlug(){

  // variables for the mother volume
  G4double NaIPlugMotherHLength = 184.095*mm + 0.5*mm;
  G4double NaIPlugMotherRadius = 59*mm;
  G4double NaIPlugMotherZpos = -(NaIPlugMotherHLength);
  G4double NaIPlugZShift = -50*mm+SourceDetectorDist; //change this to move NaI plug (z-shift: -away, +toward Ge Crystal)

  // variables for the Al can 
  G4double naiplugcanThickness = 1.0*mm;
  //G4double naiplugcanOuterRad = NaIPlugMotherRadius-1.0*mm;
  //G4double naiplugcanInnerRad = naiplugcanOuterRad-naiplugcanThickness;
  //G4double naiplugcanHLen = 150.0*mm;

  G4double naiplugcanLength = NaIPlugMotherHLength*2.-1.*mm;
  G4double naiplugcanKink1 = 110.02*mm;
  G4double naiplugcanKink2 = 129.56*mm;
  G4double naiplugcanKink3a = 184.94*mm;
  G4double naiplugcanKink3b = 184.95*mm;
  G4double naiplugcanKink4a = 203.3*mm;
  G4double naiplugcanKink4b = 203.4*mm;
  //G4double naiplugcanbackOuterRad = (119.38/4.);
  //G4double naiplugcanbackInnerRad = naiplugcanbackOuterRad-naiplugcanThickness;

  G4double naiplugOuterRad1 = 31.155*mm;
  G4double naiplugOuterRad2 = 31.155*mm;
  G4double naiplugOuterRad3 = 31.155*mm;
  G4double naiplugOuterRad4 = 48.475*mm;
  G4double naiplugOuterRad5 = 48.475*mm;
  G4double naiplugOuterRad6 = 58.96*mm;
  G4double naiplugOuterRad7 = 58.96*mm;
  G4double naiplugOuterRad8 = 58.46*mm;
  G4double naiplugOuterRad9 = 58.460*mm;
  G4double naiplugOuterRad10 = 58.460*mm;

  G4double naiplugInnerRad1 = naiplugOuterRad1-naiplugcanThickness;
  G4double naiplugInnerRad2 = naiplugOuterRad2-naiplugcanThickness;
  G4double naiplugInnerRad3 = naiplugOuterRad3-naiplugcanThickness;
  G4double naiplugInnerRad4 = naiplugOuterRad4-naiplugcanThickness;
  G4double naiplugInnerRad5 = naiplugOuterRad5-naiplugcanThickness;
  G4double naiplugInnerRad6 = naiplugOuterRad6-naiplugcanThickness;
  G4double naiplugInnerRad7 = naiplugOuterRad7-naiplugcanThickness;
  G4double naiplugInnerRad8 = naiplugOuterRad8-naiplugcanThickness;
  G4double naiplugInnerRad9 = naiplugOuterRad9-naiplugcanThickness;
  G4double naiplugInnerRad10 = naiplugOuterRad10-naiplugcanThickness;

  // variables for the NaI(Tl) detector
  G4double naiplugRad =  naiplugInnerRad8-0.70*mm;
  G4double naiplugHLen = (177.0/2.)*mm;

  ///// Create a mother volume to build the NaI(Tl) plug in //////
  G4Tubs* NaIPlugMotherCylinder_solid
    = new G4Tubs("NaIPlug_Mother",0,
  		 NaIPlugMotherRadius,NaIPlugMotherHLength,0,twopi);
  G4LogicalVolume* NaIPlugMotherCylinder_log =
    new G4LogicalVolume(NaIPlugMotherCylinder_solid, vacuum, "NaI Plug Mother log",
			0,0,0);
  new G4PVPlacement(0,
		    G4ThreeVector(0.*mm,
				  0.*mm,
				  NaIPlugMotherZpos+NaIPlugZShift),
		      NaIPlugMotherCylinder_log,
		      "NaI Plug Detector Mother Volume",
		      LENALab_log,
		      false,
		      0);

  ////// Build the Al can to hold the NaI(Tl) detector //////////

  const G4int Nnaiplugz = 10;//number of sections 
  //G4double naiplugcanrz[Nnaiplugz];//z position of each "kink"
  //G4double naiplugcanrouter[Nnaiplugz];//inner radius at each "kink"
  //G4double naiplugcanrinner[Nnaiplugz];//outer radius at each "kink"
  
  G4double naiplugcanrz[] = {0*mm,naiplugcanThickness,naiplugcanKink1,naiplugcanKink2,naiplugcanKink3a,naiplugcanKink3b,naiplugcanKink4a,naiplugcanKink4b,naiplugcanLength-naiplugcanThickness,naiplugcanLength};
  G4double naiplugcanrinner[] = {naiplugInnerRad1,naiplugInnerRad2,naiplugInnerRad3,naiplugInnerRad4,naiplugInnerRad5,naiplugInnerRad6,naiplugInnerRad7,naiplugInnerRad8,naiplugInnerRad9,naiplugInnerRad10};
  G4double naiplugcanrouter[] = {naiplugOuterRad1,naiplugOuterRad2,naiplugOuterRad3,naiplugOuterRad4,naiplugOuterRad5,naiplugOuterRad6,naiplugOuterRad7,naiplugOuterRad8,naiplugOuterRad9,naiplugOuterRad10};

  G4Polycone *NaIPlugCan_solid
    = new G4Polycone("NaI Plug Can",0,twopi,Nnaiplugz,naiplugcanrz,naiplugcanrinner,naiplugcanrouter);
  G4LogicalVolume *NaIPlugCan_log = 
    new G4LogicalVolume(NaIPlugCan_solid,Al,"NaI PLug Can log",0,0,0);
  new G4PVPlacement(0,
		    G4ThreeVector(0.*mm,
				  0.*mm,
				  NaIPlugMotherZpos),
		    NaIPlugCan_log,
		    "NaI Plug Can",
		    NaIPlugMotherCylinder_log,
		    false,
		    0);

 

  ///////// Build and place the NaI(Tl) detector ////////////////

  G4Tubs* NaIPlug_solid 
    = new G4Tubs("NaI Plug", 0.*mm, naiplugRad, naiplugHLen, 0,twopi);
  NaIPlug_log = 
    new G4LogicalVolume(NaIPlug_solid,NaITl,"Active NaI Crystal log",0,0,0);
  new G4PVPlacement(0,
		    G4ThreeVector(0.*mm,
				  0.*mm,
				  NaIPlugMotherHLength-naiplugHLen-naiplugcanThickness-1.*mm),//1mm back from front edge of NaI plug can
		    NaIPlug_log,
		    "NaI Plug",
		    NaIPlugMotherCylinder_log,
		    false,
		    16);

  NaIPlugMotherCylinder_log->SetVisAttributes(G4VisAttributes::Invisible);

  G4VisAttributes* naiplugcanVisAtt = new G4VisAttributes(G4Colour(0.8,0.8,0.8));
  NaIPlugCan_log->SetVisAttributes(naiplugcanVisAtt);

  G4VisAttributes* naiplugVisAtt = new G4VisAttributes(G4Colour(0.,0.75,0.));
  NaIPlug_log->SetVisAttributes(naiplugVisAtt);


}


void DetectorConstruction::BuildMeteorite(){

  // variables for the mother volume
  G4double SourceMotherRadius = 116./2.*mm;
  G4double SourceMotherHLength = (48.71/2.)*mm;
  G4double SourceMotherZpos = -SourceMotherHLength+SourceDetectorDist;

  // variables for the styrofoam holder
  G4double sfRad = 115.96/2.*mm;
  G4double sfHLen = (48.7/2.)*mm;

  G4RotationMatrix noRotate;
  
  ///////  Build the mother volume to put the Styrofoam holder and Meteorite in //////

  G4Tubs* SourceMotherCylinder_solid
    = new G4Tubs("Source_Mother",0,
  		 SourceMotherRadius,SourceMotherHLength,0,twopi);
  G4LogicalVolume* SourceMotherCylinder_log =
    new G4LogicalVolume(SourceMotherCylinder_solid, Air, "Meteorite Mother log",
			0,0,0);
  new G4PVPlacement(0,
		    G4ThreeVector(0.*mm,
				  0.*mm,
				  SourceMotherZpos),
		      SourceMotherCylinder_log,
		      "Meteorite Mother Volume",
		      LENALab_log,
		      false,
		      0);

  ///// Build the stryofoam holder /////////////////

  G4Tubs* Styrofoam_solid 
    = new G4Tubs("Styrofoam_holder", 0.*mm, sfRad, sfHLen, 0,twopi);
  G4LogicalVolume *Styrofoam_log = 
    new G4LogicalVolume(Styrofoam_solid,styrofoam,"Stryofoam log",0,0,0);
  new G4PVPlacement(0,
		    G4ThreeVector(0.*mm,
				  0.*mm,
				  0.*mm),
		    Styrofoam_log,
		    "Stryrofoam holder",
		    SourceMotherCylinder_log,
		    false,
		    0);

  /////// Build the Meteorite and place in the holder ///////////

  //Add rock (probably need G4Intersection...sphere+box)

  //rock geometry

  //G4double pxsemiaxis = (10.6/2.)*cm;
  //G4double pysemiaxis = (10.6/2.)*cm;
  //G4double pzsemiaxis = (9.1/2.)*cm;
  //G4double pzBottomCut = 0;
  //G4double pzTopCut = 0;

  G4double R1 = (0.)*cm;//Radius at -Dz
  G4double R2 = (10.6)*cm;//Radius at +Dz greater than R1
  G4double Dz = (4.6/2.)*cm;//Half length Z
  //G4Ellipsoid* MetEll_solid 
  //  = new G4Ellipsoid("Meteorite1", pxsemiaxis, pysemiaxis, pzsemiaxis, pzBottomCut,pzTopCut);
  G4Paraboloid *MetEll_solid
    = new G4Paraboloid("Meteorite1",
		       Dz,
		       R1,
		       R2);
  
  G4double pxsemiaxis = R2;
  G4double pysemiaxis = R2;

  G4double metboxlx = (7)*cm;
  G4double metboxly = (7)*cm;
  //G4double metboxlz = (7)*cm;//thickness of rock

  G4double metboxx = (7.901/2.)*cm;
  G4double metboxy = (8.871/2.)*cm;
  G4double metboxz = (4.620/2.)*cm;//thickness of rock
  G4Box *MetBox_solid
    = new G4Box("Meteorite2",metboxlx,metboxly,metboxz);
 
  G4double mettransz = 0.;//pzsemiaxis-(metboxz*2); 
  G4double mettransx = -((pxsemiaxis+metboxlx)-(metboxx*2));
  G4double mettransy = -((pysemiaxis+metboxly)-(metboxy*2));
 
  G4ThreeVector mettrans(mettransx,mettransy,mettransz);

  G4IntersectionSolid* Meteorite_solid 
    = new G4IntersectionSolid("Meteorite3", MetBox_solid, MetEll_solid,&noRotate,mettrans);  

  G4LogicalVolume *Meteorite_log = 
    new G4LogicalVolume(Meteorite_solid,meteorite,"Metorite log",0,0,0);
 
  G4RotationMatrix* Met_rot = new G4RotationMatrix();
  Met_rot->rotateX(0*deg);
  Met_rot->rotateY(0*deg);

  G4double xoffset = -1.5*cm;
  G4double yoffset = -1.5*cm;
  
  new G4PVPlacement(Met_rot,//rotation
		    G4ThreeVector(-mettransx/2.+xoffset,
				  -mettransy/2.+yoffset,
				  -mettransz/2.),
		    Meteorite_log,
		    "Meteorite",
		    Styrofoam_log,
		    false,
		    0);
 
  G4VisAttributes* metVisAtt = new G4VisAttributes(G4Colour(0.,0.,0.75));
  Meteorite_log->SetVisAttributes(metVisAtt);
  
  G4VisAttributes* sfVisAtt = new G4VisAttributes(G4Colour(0.,1.,1.));
  Styrofoam_log->SetVisAttributes(sfVisAtt);

}

void DetectorConstruction::BuildGasTarget(){

 

  // Target geometry
  /*
  G4Box* Wall_solid 
    = new G4Box("Walls",WallHLengthX,WallHLengthY,WallHLengthZ);
  G4LogicalVolume* Wall_log = 
    new G4LogicalVolume(Wall_solid,concrete,"Wall Solid", 0,0,0);
  new G4PVPlacement(0,
		    G4ThreeVector(0.*mm,
				  0.*mm,
				  0.*mm),
		    Meteorite_log,
		    "Meteorite",
		    Styrofoam_log,
		    false,
		    0);
  */
  // Source geometry

 
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
  G4ThreeVector SSteelSlab(0,3.44*cm,0);
  G4ThreeVector QuartzHoleL(0,0,-28.05*cm);
  G4ThreeVector QuartzHoleR(0,0,28.05*cm);
  G4ThreeVector APEXradius(0,-24.5*cm,0);
  G4ThreeVector AngCorradius(0,0,24.5*cm);
  G4RotationMatrix noRotate;

  //----------------Scintillator Casing--------------

  G4Box* scintcasing_box
    = new G4Box("SSteel Slab",HalfSlabWidth,HalfSlabHeight,HalfSlabLength);

  G4Tubs* cyl_hole
    = new G4Tubs("Cylindrical Hole",0.*mm,QuartzRad,QuartzThick,0,twopi);

  G4Trap* ssteelcasing_trap
    = new G4Trap("NaI Casing",HalfTrapLength+10.0*mm,Theta,Phi,HalfTrapHeight+0.4*mm,
		HalfTrapWidthShort+0.4*mm,HalfTrapWidthLong+0.5*mm,
		Alpha,HalfTrapHeight+0.4*mm,HalfTrapWidthShort+0.4*mm,
		HalfTrapWidthLong+0.5*mm,Alpha);

  G4Box* volume_box
    = new G4Box("Slab Volume",HalfSlabWidth,HalfSlabHeight,HalfSlabLength+2.0*mm);

  G4Trap* volume_trap
    = new G4Trap("Trap Volume",HalfTrapLength+12.0*mm,Theta,Phi,HalfTrapHeight+0.4*mm,
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
    = new G4UnionSolid("Scintillator Volume",volume_trap,volume_box,
		       &noRotate,SSteelSlab);

  G4LogicalVolume* ScintCasing_log
    = new G4LogicalVolume(scintcasing, ssteel, "SSteel NaI Casing", 0, 0, 0);

  G4LogicalVolume* ScintVolume_log
    = new G4LogicalVolume(scintvolume, Air, "Air Volume", 0, 0, 0);

  // place SSteel NaI Casing inside the Segment Volume
  new G4PVPlacement(0,
                        G4ThreeVector(),		// at (0,0,0)
                        ScintCasing_log,           	// Logical Volume
                        "Segment Volume",           	// Name
                        ScintVolume_log,           	// Mother Volume
                        false,                     	// no boolean operations
                        0);                        	// Copy number

  //THIS IS IMPORTANT. I DON'T KNOW WHY THERE IS NO COMMENT EXPLAINING WHAT THIS IS.
  for(G4int i=0;i<24;i++){
    NaIsegRotate=(i*(360/24))*deg;
    new G4PVPlacement(G4Transform3D(G4RotationMatrix(zAxis,-NaIsegRotate),
                      G4ThreeVector(Radius*std::sin(NaIsegRotate),Radius*std::cos(NaIsegRotate),0)),
		      ScintVolume_log,		// Logical Volume
		      "APEX Segment ",		// Name
		      LENALab_log,		// Mother Volume
		      false,			// no boolean operations
		      i);			// Copy number
}
/*
  // now place the Casing Volume inside the LENALab mother volume
  // NOTE: rotation and position are important
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(zAxis,180*deg),
                        APEXradius),			// at radius of APEX
                        ScintVolume_log,           	// Logical Volume
                        "NaI Segment 1",           	// Name
                        LENALab_log,           	   	// Mother Volume
                        false,                     	// no boolean operations
                        0);                        	// Copy number

  // now place the Casing Volume inside the LENALab mother volume
  // NOTE: rotation and position are important
  new G4PVPlacement(G4Transform3D(G4RotationMatrix(180*deg,90*deg,90*deg),
                        AngCorradius),			// at radius of angular correlation measurement
                        ScintVolume_log,           	// Logical Volume
                        "NaI Segment 1",           	// Name
                        LENALab_log,           	   	// Mother Volume
                        false,                     	// no boolean operations
                        0);                        	// Copy number
*/

  //------------------Quartz Windows-----------------  

  G4Tubs* quartz_cyl
    = new G4Tubs("Quartz Cylinder",0.*mm,QuartzRad,QuartzThick,0,twopi);

  G4LogicalVolume* Quartz_log
    = new G4LogicalVolume(quartz_cyl, Quartz,"Quartz Window",0,0,0);

  new G4PVPlacement(0,                             		// no rotation
                        QuartzHoleL,				// at (0,0,-28.05*cm)
                        Quartz_log,          	   		// Logical Volume
                        "Quartz Window L",            		// Name
                        ScintVolume_log,	           	// Mother Volume
                        false,                     		// no boolean operations
                        0);                        		// Copy number

  new G4PVPlacement(0,                             		// no rotation
                        QuartzHoleR,				// at (0,0,28.05*cm)
                        Quartz_log,          	   		// Logical Volume
                        "Quartz Window R",            		// Name
                        ScintVolume_log,	           	// Mother Volume
                        false,                     		// no boolean operations
                        1);                        		// Copy number

  //------------------Scintillators------------------

  G4Trap* scintillator_trap
    = new G4Trap("NaI Crystal",HalfTrapLength,Theta,Phi,HalfTrapHeight,
			HalfTrapWidthShort,HalfTrapWidthLong,Alpha,HalfTrapHeight,
			HalfTrapWidthShort,HalfTrapWidthLong,Alpha);

//  G4LogicalVolume* APEXseg_log
	APEXNaI_log
    = new G4LogicalVolume(scintillator_trap, NaITl, "APEX NaI Scintillator", 0, 0, 0);

  //G4VPhysicalVolume* scintillator_phys
  new G4PVPlacement(0,                             // no rotation
                        G4ThreeVector(),           // at (0,0,0)
                        APEXNaI_log,		   // Logical Volume
                        "NaI Scintillator",        // Name
                        ScintCasing_log,           // Mother Volume
                        false,                     // no boolean operations
                        0);                        // Copy number

  //------------------------------------------------------------
  //-------------- Visualization Stuff -----------------------
  //------------------------------------------------------------

  // The scintillator volume is invisible
  ScintVolume_log ->SetVisAttributes(G4VisAttributes::Invisible);

  // Scintillator casing is dark gray
  G4VisAttributes* ssteelVisAtt = new G4VisAttributes(G4Colour(0.4,0.4,0.4));
  ScintCasing_log->SetVisAttributes(ssteelVisAtt);

  // NaI(Tl) scintillator is yellow
  G4VisAttributes* detVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  APEXNaI_log->SetVisAttributes(detVisAtt);

  // Quartz window is white
  G4VisAttributes* quartzVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  Quartz_log->SetVisAttributes(quartzVisAtt);

}


//#####################################################################
// The APEX Cradle
void DetectorConstruction::BuildAPEXcradle(){

  // The Lead Shield
  G4double LeadShield_outer = 25.75*25.4/2.*mm;
  G4double LeadShield_inner = 24.25*25.4/2.*mm;
  G4double LeadShield_Hlen  = 23.5*25.4/2.*mm;

  G4Tubs* LeadShield_solid
    = new G4Tubs("APEX Shield",LeadShield_inner,LeadShield_outer,
		 LeadShield_Hlen,0,twopi);
  G4LogicalVolume* LeadShield_log =
    new G4LogicalVolume(LeadShield_solid, Pb,"Lead Shield",0,0,0);

    new G4PVPlacement(0,
		      G4ThreeVector(0,0,0),
		      LeadShield_log,
		      "APEX Lead Shield",
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

  G4LogicalVolume* APEXcradle_log =
    new G4LogicalVolume(AlCradle_solid, Al,"Al Cradle",0,0,0);

    new G4PVPlacement(0,
		      G4ThreeVector(0,0,0),
		      APEXcradle_log,
		      "APEX Cradle",
		      LENALab_log,
		      false,
		      0);


  // subtract material for al cradle wings
  G4SubtractionSolid* cradlewings_solid
    = new G4SubtractionSolid("APEX Wings",Al_box,sub_cradle,
			     &noRotate,origin);

  G4LogicalVolume* APEXwings_log =
    new G4LogicalVolume(cradlewings_solid, Al,"Al Wings",0,0,0);

    new G4PVPlacement(0,
		      wingL1,
		      APEXwings_log,
		      "APEX Wing L1",
		      LENALab_log,
		      false,
		      0);
    new G4PVPlacement(0,
		      wingR1,
		      APEXwings_log,
		      "APEX Wing R1",
		      LENALab_log,
		      false,
		      1);
    new G4PVPlacement(0,
		      wingL2,
		      APEXwings_log,
		      "APEX Wing L2",
		      LENALab_log,
		      false,
		      2);
    new G4PVPlacement(0,
		      wingR2,
		      APEXwings_log,
		      "APEX Wing R2",
		      LENALab_log,
		      false,
		      1);

  //------------------------------------------------------------
  //-------------- Visualization Stuff -----------------------
  //------------------------------------------------------------

  // Al Cradle is Blue
  G4VisAttributes* CradleVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  APEXcradle_log->SetVisAttributes(CradleVisAtt);
  APEXwings_log->SetVisAttributes(CradleVisAtt);

  // Lead is Black
  G4VisAttributes* LeadVisAtt = new G4VisAttributes(G4Colour(0.1,0.1,0.1));
  LeadShield_log->SetVisAttributes(LeadVisAtt);

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
    new G4LogicalVolume(CollimatorPipe_solid, Al,"Al Collimator Pipe",0,0,0);
  // G4VPhysicalVolume* CollimatorPipe_phys = 
    new G4PVPlacement(0,
		      G4ThreeVector(0,0,0),
		      CollimatorPipe_log,
		      "APEX Collimator Pipe",
		      LENALab_log,
		      false,
		      0);

  // The Delrin Container
  G4double DelrinContainer_outer = 4.25*25.4/2.*mm;
  G4double DelrinContainer_inner = 4.0*25.4/2.*mm;
  G4double DelrinContainer_Hlen  = 4.14*25.4/2.*mm;

  G4Tubs* DelrinContainer_solid
    = new G4Tubs("Collimator Container",DelrinContainer_inner,DelrinContainer_outer,
		 DelrinContainer_Hlen,0,twopi);
  G4LogicalVolume* DelrinContainer_log =
    new G4LogicalVolume(DelrinContainer_solid, delrin,"Delrin Container",0,0,0);
  // G4VPhysicalVolume* DelrinContainer_phys = 
    new G4PVPlacement(0,
		      G4ThreeVector(0,0,0),
		      DelrinContainer_log,
		      "Delrin Collimator Container",
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
    = new G4Tubs("Collimator Cylinder",0.*mm,LeadRad,LeadThick,0,twopi);

  G4LogicalVolume* LeadCyl_log
    = new G4LogicalVolume(lead_cyl, Pb,"Lead Cylinder",0,0,0);

  new G4PVPlacement(0,                             		// no rotation
                        lead_cylL,				// at (0,0,-31.75*mm)
                        LeadCyl_log,          	   		// Logical Volume
                        "LeadCyl L",            		// Name
                        LENALab_log,	           		// Mother Volume
                        false,                     		// no boolean operations
                        0);                        		// Copy number

  new G4PVPlacement(0,                             		// no rotation
                        lead_cylR,				// at (0,0,25.4*mm)
                        LeadCyl_log,          	   		// Logical Volume
                        "LeadCyl R",            		// Name
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
  LeadCyl_log->SetVisAttributes(LeadVisAtt);



}
void DetectorConstruction::BuildAero(){

  G4String name, symbol;             // a=mass of a mole;
  G4double a, z, density;            // z=mean number of protons;  
  //  G4int iz, n;                       // iz=nb of protons  in an isotope; 
  // n=nb of nucleons in an isotope;
  G4int ncomponents, natoms;
  //  G4double abundance, fractionmass;
  G4double fractionmass;
  //G4double temperature, pressure;


  // define Elements
  a = 1.01*g/mole;
  G4Element* elH  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a);

  a = 12.01*g/mole;
  G4Element* elC  = new G4Element(name="Carbon"  ,symbol="C" , z= 6., a);

  a = 16.00*g/mole;
  G4Element* elO  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);

  a = 28.09*g/mole;
  G4Element* elSi = new G4Element(name="Silicon", symbol="Si", z=14., a);


  // define a material from elements.   case 1: chemical molecule
  density = 1.000*g/cm3;
  G4Material* H2O = new G4Material(name="Water", density, ncomponents=2);
  H2O->AddElement(elH, natoms=2);
  H2O->AddElement(elO, natoms=1);

  density = 1.032*g/cm3;
  G4Material* Sci = new G4Material(name="Scintillator", density, ncomponents=2);
  Sci->AddElement(elC, natoms=9);

  // define a material from elements and/or others materials (mixture of mixtures)
  density = 0.200*g/cm3;
  G4Material* Aerog = new G4Material(name="Aerogel", density, ncomponents=3);
  Aerog->AddMaterial(SiO2, fractionmass=62.5*perCent);
  Aerog->AddMaterial(H2O , fractionmass=37.4*perCent);
  Aerog->AddElement (elC , fractionmass= 0.1*perCent);
  Sci->AddElement(elH, natoms=10);

  density = 2.200*g/cm3;
  G4Material* SiO2 = new G4Material(name="quartz", density, ncomponents=2);
  SiO2->AddElement(elSi, natoms=1);
  SiO2->AddElement(elO , natoms=2);


}
void DetectorConstruction::BuildMagnet(){

  double halfMagnetLength=13.75*cm;
  double halfMagnetGap=10.0*cm;
  double translateMagnet=halfMagnetLength+halfMagnetGap;
  double innerRadius=2.0*cm;
  double outerRadius=10.0*cm;

  G4Tubs* magSolid = new G4Tubs("magSolid",innerRadius,outerRadius,halfMagnetLength,0,360);

  G4LogicalVolume* magSolid_log =
    new G4LogicalVolume(magSolid,steel1010A, "Magnet log");

  new G4PVPlacement(0,
		    G4ThreeVector(0,0,translateMagnet),
		    magSolid_log,
		    "Source",
		    LENALab_log,
		    false,
		    0);

  new G4PVPlacement(0,
		    G4ThreeVector(0,0,-translateMagnet),
		    magSolid_log,
		    "Source",
		    LENALab_log,
		    false,
		    0);

  G4VisAttributes* MagVisAtt = new G4VisAttributes(G4Colour(0.8,0.0,0.8));
  MagVisAtt->SetForceSolid(1);
  magSolid_log->SetVisAttributes(MagVisAtt);

  //Bore hole in magnet

}

void DetectorConstruction::BuildSource(){

  double pipeLength=30.0*cm;
  double halfGap=30.0*cm;
  double halfPipeLength=pipeLength/2.0;
  double shiftPipe=0;
  double pipeRadius=1.0*cm;
  G4double a; //atomic mass
  G4double z; //atomic number
  G4double density, fractionmass;
  G4int ncomponents,natoms;

  Al = 
     new G4Material("Aluminum", z= 13., a= 26.98*g/mole, density= 2.7*g/cm3);

  G4Tubs* sourcePipeSolid = new G4Tubs("sourcePipeSolid",0,pipeRadius,halfPipeLength,0,360);

  G4LogicalVolume* sourcePipeSolid_log =
    new G4LogicalVolume(sourcePipeSolid,Al, "sourcePipe log");

  /*  new G4PVPlacement(0,
		    G4ThreeVector(0,0,shiftPipe),
		    sourcePipeSolid_log,
		    "Source",
		    LENALab_log,
		    false,
		    0);


  new G4PVPlacement(0,
		    G4ThreeVector(0,0,-shiftPipe),
		    sourcePipeSolid_log,
		    "Source",
		    LENALab_log,
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

  G4double innerSourceTubeRadius=4.5*cm;
  G4double outerSourceTubeRadius=5.0*cm;
  G4double sourceTubeLength=5*cm;
  //Hollow source container
  G4Tubs* sourceTube_1 = new G4Tubs("sourcePipeSolid",innerSourceTubeRadius,outerSourceTubeRadius,sourceTubeLength,3.14159,6.3);


  G4LogicalVolume* sourceTube_log_1 =
    new G4LogicalVolume(sourceTube_1,steel1010A, "sourcePipe log");

    new G4PVPlacement(0,
		    G4ThreeVector(0,0,0),
		    sourceTube_log_1,
		    "Source",
		    LENALab_log,
		    false,
		    0);

  // Aluminium Tube Color is red
  G4VisAttributes* TubeVisAtt = new G4VisAttributes(G4Colour(0.8,0.0,0.0));
  TubeVisAtt->SetForceSolid(1);
  sourceTube_log_1->SetVisAttributes(TubeVisAtt);

  //=============================================

  //=============================================
  // Source Endcaps
  //=============================================

  G4double innerSourceEndcapRadius=0*cm;
  G4double outerSourceEndcapRadius=5*cm;
  G4double sourceEndcapLength=1.0*cm;

  G4Tubs* sourceTubeEndcap_1 = new G4Tubs("sourcePipeSolid",innerSourceEndcapRadius,outerSourceEndcapRadius,sourceEndcapLength,3.14159,6.3);
  G4Tubs* sourceTubeEndcap_2 = new G4Tubs("sourcePipeSolid",innerSourceEndcapRadius,outerSourceEndcapRadius,sourceEndcapLength,3.14159,6.3);

  G4LogicalVolume* sourceTubeEndcap_log_1 =
    new G4LogicalVolume(sourceTubeEndcap_1,steel1010A, "sourcePipe log");

  G4LogicalVolume* sourceTubeEndcap_log_2 =
    new G4LogicalVolume(sourceTubeEndcap_2,steel1010A, "sourcePipe log");


  new G4PVPlacement(0,
		    G4ThreeVector(0,0,-6*cm),
		    sourceTubeEndcap_log_1,
		    "Source",
		    LENALab_log,
		    false,
		    0);


  new G4PVPlacement(0,
		    G4ThreeVector(0,0,6*cm),
		    sourceTubeEndcap_log_2,
		    "Source",
		    LENALab_log,
		    false,
		    0);

  sourceTubeEndcap_log_1->SetVisAttributes(TubeVisAtt);

  sourceTubeEndcap_log_2->SetVisAttributes(TubeVisAtt);
  
  //===========================================================





  G4String name, symbol;             // a=mass of a mole;

  // n=nb of nucleons in an isotope;


  // define Elements
  a = 1.01*g/mole;
  G4Element* elH  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a);
  a = 12.01*g/mole;
  G4Element* elC  = new G4Element(name="Carbon"  ,symbol="C" , z= 6., a);
  a = 16.00*g/mole;
  G4Element* elO  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);
  a = 28.09*g/mole;
  G4Element* elSi = new G4Element(name="Silicon", symbol="Si", z=14., a);

  // define a material from elements.   case 1: chemical molecule
  density = 1.000*g/cm3;
  G4Material* H2O = new G4Material(name="Water", density, ncomponents=2);
  H2O->AddElement(elH, natoms=2);
  H2O->AddElement(elO, natoms=1);

  density = 1.032*g/cm3;
  G4Material* Sci = new G4Material(name="Scintillator", density, ncomponents=2);
  Sci->AddElement(elC, natoms=9);

  // define a material from elements and/or others materials (mixture of mixtures)
  density = 0.200*g/cm3;
  G4Material* Aerog = new G4Material(name="Aerogel", density, ncomponents=3);
  Aerog->AddMaterial(SiO2, fractionmass=62.5*perCent);
  Aerog->AddMaterial(H2O , fractionmass=37.4*perCent);
  Aerog->AddElement (elC , fractionmass= 0.1*perCent);
  Sci->AddElement(elH, natoms=10);

  density = 2.200*g/cm3;
  G4Material* SiO2 = new G4Material(name="quartz", density, ncomponents=2);
  SiO2->AddElement(elSi, natoms=1);
  SiO2->AddElement(elO , natoms=2);

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

  //size of aerogel 'cube' which is not a disk

  G4Tubs* aeroDisk = new G4Tubs("aeroDiskSolid",0.0*cm,45.0*mm,4.0*cm,0,6.3);

  G4LogicalVolume* aeroDisk_log =
    new G4LogicalVolume(aeroDisk,vacuum, "Aerog log");

  double halfAerogelLength=1.0*cm;
  double halfAerogelGap=1.00*cm;
  double aerogelShift=1.5*cm+halfAerogelGap;

  new G4PVPlacement(0,
		    G4ThreeVector(0,0,aerogelShift),
		    aeroDisk_log,
		    "Source",
		    LENALab_log,
		    false,
		    0);

  new G4PVPlacement(0,
		    G4ThreeVector(0,0,-aerogelShift),
		    aeroDisk_log,
		    "Source",
		    LENALab_log,
		    false,
		    0);
  //Aero Color
  G4VisAttributes* AeroColor = new G4VisAttributes(G4Colour(0.0,0.8,0.8));
  AeroColor->SetForceSolid(1);
  aeroDisk_log->SetVisAttributes(AeroColor);

  //Scintillator
  double halfScintWidth=0.25*cm;
  double halfScintGap=0.1*cm;
  double scintShift=halfScintWidth+halfScintGap;

  G4double outerScintRadius = 31.0*mm;
  G4double innerScintRadius = 0.0*mm;

  G4NistManager* manager = G4NistManager::Instance();
  anthracene  = manager->FindOrBuildMaterial("G4_ANTHRACENE");

  G4Tubs* scintDisk = new G4Tubs("scintDiskSolid",innerScintRadius,outerScintRadius,halfScintWidth,0,6.3);

  G4LogicalVolume* scintDisk_log =
    new G4LogicalVolume(scintDisk,vacuum, "scint log");

  new G4PVPlacement(0,
		    G4ThreeVector(0,0,scintShift),
		    scintDisk_log,
		    "Source",
		    LENALab_log,
		    false,
		    0);
  
  new G4PVPlacement(0,
		    G4ThreeVector(0,0,-scintShift),
		    scintDisk_log,
		    "Source",
		    LENALab_log,
		    false,
		    0);

  // Scintillator Color is Green
  G4VisAttributes* ScintVisAtt = new G4VisAttributes(G4Colour(0.0,0.8,0.0));
  ScintVisAtt->SetForceSolid(1);
  scintDisk_log->SetVisAttributes(ScintVisAtt);

// Kapton Dupont de Nemur (density: 1.396-1.430, get middle )

  a = 14.01*g/mole;
  G4Element* elN  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);
  G4int nel;
  density = 1.413*g/cm3;
  G4Material* Kapton = new G4Material(name="Kapton", density, nel=4);
  Kapton->AddElement(elO,5);
  Kapton->AddElement(elC,22);
  Kapton->AddElement(elN,2);
  Kapton->AddElement(elH,10);

  double halfKaptonWidth=0.25*cm;
  double innerFoilRadius=0*mm;
  double outerFoilRadius=31*mm;
  G4Tubs* kaptonFoil = new G4Tubs("kaptonFoilSolid",innerFoilRadius,outerFoilRadius,halfKaptonWidth,0,6.3);

  G4LogicalVolume* kaptonFoil_log =
    new G4LogicalVolume(kaptonFoil,vacuum, "kapton log");

  new G4PVPlacement(0,
		    G4ThreeVector(0,0,0),
		    kaptonFoil_log,
		    "Source",
		    LENALab_log,
		    false,
		    0);
  

  // Kapton Color is Green
  G4VisAttributes* KapVisAtt = new G4VisAttributes(G4Colour(0.8,0.0,0.8));
  KapVisAtt->SetForceSolid(1);
  kaptonFoil_log->SetVisAttributes(KapVisAtt);

}
