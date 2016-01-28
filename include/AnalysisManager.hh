//----------------------------------------------------------------------
// AnalysisManager.hh
//
// Description: Sorts out the analysis, Include file
//
// Version 0.1, 8/3/05
// Changes: Nothing!
//----------------------------------------------------------------------
#ifndef AnalysisManager_h
#define AnalysisManager_h 1

//#ifdef G4ANALYSIS_USE

#include "G4ThreeVector.hh"
#include "globals.hh"
#include "g4root.hh"

class G4Track;
class AnalysisMessenger;

class AnalysisManager {
public:

  virtual ~AnalysisManager();
  static AnalysisManager* getInstance();
  static void dispose();

public:

  // All the commands 
  void SetFileName(G4String filename) {fileName = filename;};
  void SetasciiFile(G4bool asciibool){ asciiFileBool = asciibool;  };
  void SetBinWidth(G4double binwidth){binWidth = binwidth;};
  void SetDetailBool(G4bool detailbool){detailBool = detailbool;};
  void SetResolutionBool(G4bool resolutionbool){resolutionBool = resolutionbool;};
  void SetGeResolution(G4double geresolution){GeResolution = geresolution*keV;}
  void SetNaIResolution(G4double nairesolution){NaIResolution = nairesolution;};
  void SetNaIPosResolution(G4double naiposresolution){NaIPosResolution = naiposresolution;};
  void SetSciResolution(G4double sciresolution){SciResolution = sciresolution;};

  // Get what the commands are set to
  G4bool GetasciiFileState(){ return asciiFileBool ;};
  G4double GetBinWidth(){ return binWidth;};
  G4bool GetDetailBool(){return detailBool;};
  G4bool GetResolutionBool(){return resolutionBool;};
  G4double GetGeResolution(){return GeResolution;};
  G4double GetNaIResolution(){return NaIResolution;};
  G4double GetNaIPosResolution(){return NaIPosResolution;};
  G4double GetSciResolution(){return SciResolution;};
  G4String GetFileName(){return fileName;};

private:

  AnalysisManager();
  static AnalysisManager* instance;

  G4String fileName;
  G4bool asciiFileBool;
  G4double binWidth;
  G4bool detailBool;
  G4bool resolutionBool;
  G4double GeResolution;
  G4double NaIResolution;
  G4double NaIPosResolution;
  G4double SciResolution;

  AnalysisMessenger* analysisMessenger;
};
//#endif // G4ANALYSIS_USE

#endif

