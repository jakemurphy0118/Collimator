#ifndef ELECTROMAGNET_H
#define ELECTROMAGNET_H 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Navigator.hh"
#include <string>

class Electromagnet
{

public:

  Electromagnet();
  ~Electromagnet();
  void BuildCALIOPEMagnetLog(G4LogicalVolume*);

private:
  void DefineMaterials();

};

#endif
