//--------------------------------------------------------------------
// DetectorConstruction.hh
//
// Description: The detector definitions, materials etc.
// Changes: 7/18/05 None yet
//--------------------------------------------------------------------

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Navigator.hh"

class G4Box;
class G4Tubs;
class G4UnionSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class DetectorMessenger;

class GermSD;
class NaIannSD;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:

  DetectorConstruction();
  ~DetectorConstruction();

public:

  G4VPhysicalVolume* Construct();

  // The functions to build various parts of the setup
  void BuildHPGe();
  void BuildNaIAnnulus();
  void BuildGammaCube();
  void BuildAPEX();
  void BuildAPEXcradle();
  void BuildTargetHolder();
  void BuildCollimator();
  void BuildRadioactiveSource();

  // Functions to set dimensions of user definable parts
  void SetGeHLength(G4double);
  void SetGeRadius(G4double);
  void SetGeOffset(G4double);
  void SetGeEndRad(G4double);
  void SetGeHoleHLength(G4double);
  void SetGeHoleRad(G4double);
  void SetGeFingerRad(G4double);
  void SetGeVertDisp(G4double);
  void SetGeHorisDisp(G4double);
  void SetGeDeadLayerHThick(G4double);
  void SetVisuals(G4bool);
  void SetSourceDetectorDist(G4double);
  void SetDetectorAngle(G4double);
  void SetContactPinType(G4int);

  // Update the geometry when dimensions have been changed
  void UpdateGeometry();
  G4bool SDflag;

  G4Navigator *aNavigator;
  G4String GetVolumeMaterial(G4double x,G4double y,G4double z);

private:
  void DefineMaterials();

  G4VPhysicalVolume* LENALab_phys;
  G4LogicalVolume* LENALab_log;
  G4LogicalVolume* ActGeCrys_log;
  G4LogicalVolume* GeCrys_log;
  G4LogicalVolume* NaISeg_log;
  G4LogicalVolume* scintillator_log;

  // Materials
  G4Material* Al;
  G4Material* Pb;
  G4Material* Cu;
  G4Material* Ta;
  G4Material* Ge;
  G4Material* Li;
  G4Material* B;
  G4Material* Air;
  G4Material* Water;
  G4Material* Quartz;
  G4Material* ceramic;
  G4Material* NaITl;
  G4Material* Sci;
  G4Material* ssteel;
  G4Material* mylar;
  G4Material* brass;
  G4Material* vacuum;
  G4Material* delrin;
  G4Material* plastic;

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

  // The Sensitive Detectors
  GermSD* GeSD;
  NaIannSD* NaISD;

  DetectorMessenger* detectorMessenger; // pointer to messenger

};

#endif
