//--------------------------------------------------------------------
// DetectorConstruction.hh
//
// Description: The detector definitions, materials etc.
// Changes: 7/18/05 None yet
//--------------------------------------------------------------------

#ifndef SphericalSourceHolderConstruction_h
#define SphericalSourceHolderConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Navigator.hh"
#include <string>

class G4Box;
class G4Tubs;
class G4UnionSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;

class GermSD;
class NaIannSD;
class ScintSD;
class APEXSD;

class SphericalSourceHolderConstruction 
{
public:

  SphericalSourceHolderConstruction();
  ~SphericalSourceHolderConstruction();

public:

  void BuildSource(G4LogicalVolume*);
  G4bool SDflag;
  G4Navigator *aNavigator;

private:
  void DefineMaterials();

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
  G4Material* Kapton;
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

  G4String sourceHolder;


};


#endif
