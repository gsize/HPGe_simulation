//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef HPGeDetectorSD_h
#define HPGeDetectorSD_h 1

#include "G4VSensitiveDetector.hh"
#include "Hit.hh"

class G4Step;
class G4HCofThisEvent;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HPGeDetectorSD : public G4VSensitiveDetector
{
  public:
      HPGeDetectorSD(G4String);
     ~HPGeDetectorSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*, G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);

  private:
      HPGeDetectorHitsCollection* trackerCollection;
  G4double e_g;
  G4double l_g;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

