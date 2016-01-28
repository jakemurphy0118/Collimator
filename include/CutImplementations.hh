#ifndef CutImplementations_h
#define CutImplementations_h 1

#include "globals.hh"
#include <vector>

class CutImplementations
{
  
public:
  
  CutImplementations();
  ~CutImplementations();
  
  void DoAnalysisCuts(G4double, G4double, G4double, G4double, G4double, G4double,std::vector<G4double>,std::vector<G4double>,std::vector<G4double>,std::vector<G4int>,G4int);
  void RyanReconstruction(std::vector<G4double>,std::vector<G4double>,std::vector<G4double>,std::vector<G4int>, G4double, G4String);
  void YamazakiReconstruction(std::vector<G4double>,std::vector<G4double>,std::vector<G4double>,std::vector<G4int>,G4double, G4String);
  
private:
  G4int IDNum;

};

#endif
