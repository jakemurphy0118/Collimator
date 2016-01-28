//--------------------------------------------------------------------
// DetectorConstruction.hh
//
// Description: The detector definitions, materials etc.
// Changes: 7/18/05 None yet
//--------------------------------------------------------------------

#ifndef CylindricalSourceHolderConstruction_h
#define CylindricalSourceHolderConstruction_h 1

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

class CylindricalSourceHolderConstruction 
{
public:

  CylindricalSourceHolderConstruction();
  ~CylindricalSourceHolderConstruction();

  void BuildSource(G4LogicalVolume*);
  G4LogicalVolume* BuildVacuum(G4LogicalVolume*);
  G4LogicalVolume* BuildKapton(G4LogicalVolume*);
  G4LogicalVolume* BuildScintillator(G4LogicalVolume*);
  G4LogicalVolume* BuildAerogel(G4LogicalVolume*);
  void BuildShell(G4LogicalVolume*);

private:
  void DefineMaterials();

  G4double a;
  G4double z;
  G4double density;
  G4double fractionmass;
  G4int ncomponents;
  G4int natoms;
  G4String name;
  G4String symbol;

  // Materials
  G4Material* Al;
  G4Material* Air;
  G4Material* Water;
  G4Material* NaITl;
  G4Material* Sci;
  G4Material* Kapton;
  G4Material* anthracene;
  G4Material* SiO2;
  G4Material* vacuum;
  G4Material* Aerog;
  G4Material* delrin;
  G4Material* rubber;
  G4String sourceHolder;


};


#endif
