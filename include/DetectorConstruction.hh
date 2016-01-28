//--------------------------------------------------------------------
// DetectorConstruction.hh
//
// Description: The detector definitions, materials etc.
// Changes: 7/18/05 None yet
//--------------------------------------------------------------------

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4LogicalVolume.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Navigator.hh"
#include <string>
#include "CylindricalSourceHolderConstruction.hh"
#include "SphericalSourceHolderConstruction.hh"
#include "AsymmetricalCylindricalSourceHolderConstruction.hh"
#include "Electromagnet.hh"
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
#include "DetectorMessenger.hh"

class G4Box;
class G4Tubs;
class G4UnionSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class DetectorMessenger;

class GermSD;
class NaIannSD;
class ScintSD;
class APEXSD;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:

  DetectorConstruction();
  ~DetectorConstruction();

public:

  G4VPhysicalVolume* Construct();

  // The functions to build various parts of the crystal
  void BuildHPGe();
  void BuildWalls();
//  void BuildSource();
  void BuildPbShield();
  void BuildAPEX();
  void BuildAPEXcradle();
  void BuildAero();
 // void BuildMagnet();
//  void SetSourceHolder(G4String);
  
  // Update the geometry when dimensions have been changed
  void UpdateGeometry();
  G4bool SDflag;

  G4Navigator *aNavigator;
  G4String GetVolumeMaterial(G4double x,G4double y,G4double z);

  G4ThreeVector GetWorldXYZ();

private:
  void DefineMaterials();
  //CylindricalSourceHolderConstruction* cylindricalHolder;
  //AsymmetricalCylindricalSourceHolderConstruction* asymmetricalCylindricalHolder;
  //SphericalSourceHolderConstruction* sphericalHolder;
//  Electromagnet* electromagnet;
  
  G4LogicalVolume* LENALab_log;
  G4LogicalVolume* ActGeCrys_log;
  G4LogicalVolume* GeCrys_log;
  G4LogicalVolume *NaISeg_log;
  G4LogicalVolume *NaIPlug_log;
  G4LogicalVolume *SciBoxTop_log;
  G4LogicalVolume *SciBoxSide_log;
  G4LogicalVolume *SciBoxBeam_log;
  G4LogicalVolume *SciBoxDet_log;
  G4LogicalVolume *APEXNaI_log;

  G4UserLimits *userLimits;

  G4VPhysicalVolume* LENALab_phys;

  // Materials
  std::string  printMaterialTable;
  G4Material* C;
  G4Material* Al;
  G4Material* Pb;
  G4Material* Cu;
  G4Material* Ta;
  G4Material* Ge;
  G4Material* Li;
  G4Material* B;
  G4Material* Air;
  G4Material* Water;
  G4Material* NaITl;
  G4Material* Sci;
  G4Material* ssteel;
  G4Material* steel1010A;
  G4Material* mylar;
  G4Material* brass;
  G4Material* ceramic;
  G4Material* mumetal;
  G4Material* concrete;
  G4Material* rubber;
  G4Material* anthracene;
  G4Material* styrofoam;
  G4Material* Al2O3;
  G4Material* CaO;
  G4Material* FeO;
  G4Material* MgO;
  G4Material* P2O5;
  G4Material* SiO2;
  G4Material* Na2O;
  G4Material* NiO;
  G4Material* meteorite;
  G4Material* vacuum;
  G4Material* Quartz;
  G4Material* delrin;

  // The changeable dimensions
  G4double GeHLength;
  G4double GeRadius;
  G4double GeOffset; // The distance of Ge from the detector face (default 6mm)
  G4double GeEndRad; // The radius of the curved end
  G4double GeHoleRad;
  G4double GeHoleHLength;
  G4double GeFingerRad;
  G4double GeVertDisp;
  G4double GeHorisDisp;
  G4double GeDeadLayerHThick;
  G4bool Visuals;

  G4double SourceDetectorDist;
  G4double DetectorAngle;

  G4int ContactPinType;

  G4double LENALabx;                // x length of world
  G4double LENALaby;                // y length of world
  G4double LENALabz;                // z length of world
  G4double Source_Detectorz;         // The dist between source and det face

  G4double VolMat;
  G4double elecEkinMin;

  //needed in multiple functions
  G4double LeadHthick;
  G4double HlenLeadTopx;
  G4double xposLeadSide;

  // The Sensitive Detectors
  GermSD* GeSD;
  NaIannSD* NaISD;
  ScintSD* SciSD;
  APEXSD* ApexSD;

  DetectorMessenger* detectorMessenger; // pointer to messenger
  G4String sourceHolderName;

  //#####################################################################
// The Collimator
void BuildCollimator(){

  //---------------- The Collimator --------------------------

   G4NistManager* NISTman = G4NistManager::Instance();
   G4double density;
   G4double a;
   G4double z;

   delrin = NISTman->FindOrBuildMaterial("G4_POLYOXYMETHYLENE");
   Al = 
     new G4Material("Aluminum", z= 13., a= 26.98*g/mole, density= 2.7*g/cm3);
   Pb =
     new G4Material("Lead", z=82.,a=207.2*g/mole, density=11340.*kg/m3);

  G4double pi  =3.141592653589793238463;

  // The Aluminum Collimator Pipe
  G4double CollimatorPipe_outer = 4.5*25.4/2.*mm;
  G4double CollimatorPipe_inner = 4.25*25.4/2.*mm;
  G4double CollimatorPipe_Hlen  = 570.0/2.*mm;
  G4RotationMatrix noRotate;

  G4Tubs* CollimatorPipe_solid
    = new G4Tubs("Collimator Pipe",CollimatorPipe_inner,CollimatorPipe_outer,
		 CollimatorPipe_Hlen,0,2*pi);
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
		 DelrinContainer_Hlen,0,2*pi);
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
  G4ThreeVector lead_cylL(0,0,-25.9*mm);
  G4ThreeVector lead_cylR(0,0,25.9*mm);


  G4Tubs* lead_cyl
    = new G4Tubs("Collimator Cylinder",0.*mm,LeadRad,LeadThick,0,2*pi);

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
//####################################################################################
  
};


#endif
