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
//#include "GermSD.hh"
//#include "NaIannSD.hh"
//#include "ScintSD.hh"
//#include "APEXSD.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "G4Navigator.hh"
#include "DetectorConstruction.hh"


using namespace CLHEP;

SphericalSourceHolderConstruction::SphericalSourceHolderConstruction()
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

 
  DefineMaterials();
  SDflag=false;

  elecEkinMin = 0*eV;

  //needed in multiple funcitons
  LeadHthick = 25.4/4.*mm;
  HlenLeadTopx = 457.2/2.+2*LeadHthick;
  xposLeadSide  = HlenLeadTopx - LeadHthick;
  
  aNavigator= new G4Navigator();
 
}

SphericalSourceHolderConstruction::~SphericalSourceHolderConstruction()
{
  delete aNavigator;


 }



void SphericalSourceHolderConstruction::DefineMaterials()
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
    ("Steel1010A", density=7.87*g/cm3, ncomponents=1);
  steel1010A->AddElement(elFe,1.0);
  /*  steel1010A->AddElement(elFe,0.9918);
  steel1010A->AddElement(elMn,0.0030);
  steel1010A->AddElement(elS,0.0005);
  steel1010A->AddElement(elP,0.0004);
  steel1010A->AddElement(elC,0.0008);*/
  

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

  //Make meteorite air, for source runs
  density = 1.290*mg/cm3;
  meteorite = new G4Material("meteorite",density,ncomponents=2);
  meteorite->AddElement(elN, fractionmass=0.7);
  meteorite->AddElement(elO, fractionmass=0.3);
  
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

void SphericalSourceHolderConstruction::BuildSource(G4LogicalVolume* LENALabVolume){

  G4cout << "BUILDING SOURCE FOR THE SPHERICAL CONSTRUCTION" << G4endl;

  //===================================
  // Define vacuum
  //===================================

  G4double density=1.e-25*g/cm3;
  G4double atomicNumber=1.;
  G4double massOfMole=1.008*g/mole;
  //  G4double density=1.e-25*g/cm3;
  G4double temperature=77.*kelvin;
  G4double pressure=0.00133*pascal;
  vacuum=new G4Material("Vacuum",atomicNumber,massOfMole,density,kStateGas,temperature,pressure);

  //===================================================
  // Surround source with vacuum solid (two halves)
  //===================================================

  G4double vacuumVolumeRadius=4.701*cm;
  G4double vacuumVolumeLength=5.051*cm;

  G4Sphere* vacuumSpherething = new G4Sphere("vacuumSpherething",0,vacuumVolumeRadius,0,2*pi,0,pi);
  G4LogicalVolume* vacuumSphere_log = new G4LogicalVolume(vacuumSpherething,vacuum,"vacuumVolumeExt");

  new G4PVPlacement(0,
		   G4ThreeVector(0,0,0),
		   vacuumSphere_log,
		   "Source",
		   LENALabVolume,
		   false,
		   0);

  G4VisAttributes* VacVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,0.8));
  VacVisAtt->SetForceSolid(1);
  VacVisAtt->SetVisibility(1);
  vacuumSphere_log->SetVisAttributes(VacVisAtt);

  G4double a; //atomic mass
  G4double z; //atomic number
  G4double fractionmass;
  G4int ncomponents,natoms;

  G4RotationMatrix noRotate;

  //==================================
  // Source Tube Parts
  //==================================


 G4Sphere* aluminumShell = new G4Sphere("aluminumShell",4.3*cm,4.6*cm,0,2*pi,0,pi);
 G4LogicalVolume* aluminumShell_log = new G4LogicalVolume(aluminumShell,Al,"vacuumVolumeExt");

 new G4PVPlacement(0,
		   G4ThreeVector(0,0,0),
		   aluminumShell_log,
		   "Source",
		   vacuumSphere_log,
		   false,
		   0);

  G4VisAttributes* AlShellVisAtt = new G4VisAttributes(G4Colour(0.8,0.0,0.0));
  AlShellVisAtt->SetForceSolid(1);
  AlShellVisAtt->SetVisibility(1);
  aluminumShell_log->SetVisAttributes(AlShellVisAtt);



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


  density = 1.032*g/cm3;
  G4Material* Sci = new G4Material(name="Scintillator", density, ncomponents=2);
      Sci->AddElement(elC, natoms=9);

  // define a material from elements and/or others materials (mixture of mixtures)
  density = .010*g/cm3;
  G4Material* Aerog = new G4Material(name="Aerogel", density, ncomponents=3);
  Aerog->AddMaterial(SiO2, fractionmass=62.5*perCent);
  Aerog->AddMaterial(Water , fractionmass=37.4*perCent);
  Aerog->AddElement (elC , fractionmass= 0.1*perCent);
  Sci->AddElement(elH, natoms=10);

  G4Sphere* aeroHemisphere1 = new G4Sphere("aeroHemisphere1",0,4.22*cm,0,2*pi,0,pi/2);
  G4Sphere* aeroHemisphere = new G4Sphere("aeroHemisphere",0,4.22*cm,0,2*pi,pi/2,pi);
  G4LogicalVolume* aeroSphere1_log = new G4LogicalVolume(aeroHemisphere1,Aerog,"aeroVolumeExt1");
  G4LogicalVolume* aeroSphere_log = new G4LogicalVolume(aeroHemisphere,Aerog,"aeroVolumeExt");

  new G4PVPlacement(0,
		    G4ThreeVector(0*cm,0*cm,-2*mm),
		    aeroSphere_log,
		    "Source",
		    vacuumSphere_log,
		    false,
		    0);

  new G4PVPlacement(0,
		    G4ThreeVector(0*cm,0*cm,2*mm),
		    aeroSphere1_log,
		    "Source",
		    vacuumSphere_log,
		    false,
		    0);


  G4VisAttributes* AeroVisAtt = new G4VisAttributes(G4Colour(0.0,0.8,0.8));
  G4VisAttributes* Aero1VisAtt = new G4VisAttributes(G4Colour(0.0,0.8,0.8));
  AeroVisAtt->SetForceSolid(1);
  AeroVisAtt->SetVisibility(1);
  aeroSphere_log->SetVisAttributes(AeroVisAtt);
  Aero1VisAtt->SetForceSolid(1);
  Aero1VisAtt->SetVisibility(1);
  aeroSphere1_log->SetVisAttributes(Aero1VisAtt);

  //Scintillator
  double halfScintWidth=0.15*mm;
  double halfScintGap=0.5*mm;
  double scintShift=halfScintWidth+halfScintGap;

  G4double outerScintRadius = 4.21*cm;
  G4double innerScintRadius = 0.0*mm;

  G4NistManager* manager = G4NistManager::Instance();
  anthracene  = manager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

  G4Tubs* scintDisk = new G4Tubs("scintDiskSolid",innerScintRadius,outerScintRadius,halfScintWidth,0,2*pi);

  G4LogicalVolume* scintDisk_log =
    new G4LogicalVolume(scintDisk,anthracene, "scint log");

  /*  new G4PVPlacement(0,
		    G4ThreeVector(0,0,scintShift),
		    scintDisk_log,
		    "Source",
		    vacuumSphere_log,
		    false,
		    0);
  
  new G4PVPlacement(0,
		    G4ThreeVector(0,0,-scintShift),
		    scintDisk_log,
		    "Source",
		    vacuumSphere_log,
		    false,
		    0);*/

  // Scintillator Color is Green
  G4VisAttributes* ScintVisAtt = new G4VisAttributes(G4Colour(0.0,0.8,0.0));
  ScintVisAtt->SetForceSolid(1);
  scintDisk_log->SetVisAttributes(ScintVisAtt);

// Kapton Dupont de Nemur (density: 1.396-1.430, get middle )
  Kapton  = manager->FindOrBuildMaterial("G4_KAPTON");
  /*  a = 14.01*g/mole;
  G4Element* elN  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);
  G4int nel;
  density = 1.413*g/cm3;
  G4Material* Kapton = new G4Material(name="Kapton", density, nel=4);
  Kapton->AddElement(elO,5);
  Kapton->AddElement(elC,22);
  Kapton->AddElement(elN,2);
  Kapton->AddElement(elH,10);*/

  double halfKaptonWidth=0.02*cm;
  double innerFoilRadius=0*mm;
  double outerFoilRadius=31*mm;
  G4Tubs* kaptonFoil = new G4Tubs("kaptonFoilSolid",innerFoilRadius,outerFoilRadius,halfKaptonWidth,0,2*pi);

  G4LogicalVolume* kaptonFoil_log =
    new G4LogicalVolume(kaptonFoil,Kapton, "kapton log");

  /*  new G4PVPlacement(0,
		    G4ThreeVector(0,0,0),
		    kaptonFoil_log,
		    "Source",
		    vacuumSphere_log,
		    false,
		    0);*/
  

  // Kapton Color is Green
  G4VisAttributes* KapVisAtt = new G4VisAttributes(G4Colour(0.8,0.0,0.8));
  KapVisAtt->SetForceSolid(1);
  kaptonFoil_log->SetVisAttributes(KapVisAtt);

}
