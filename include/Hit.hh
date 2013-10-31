//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Hit_h
#define Hit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HPGeDetectorHit : public G4VHit
{
  public:

      HPGeDetectorHit();
     ~HPGeDetectorHit();
      HPGeDetectorHit(const HPGeDetectorHit&);
      const HPGeDetectorHit& operator=(const HPGeDetectorHit&);
      G4int operator==(const HPGeDetectorHit&) const;

      inline void* operator new(size_t);
      inline void  operator delete(void*);

      void Draw();
      void Print();

  public:

      void SetTrackID  (G4int track)      { trackID = track; };
      void SetEdep     (G4double de)      { edep = de; };
      void SetPos      (G4ThreeVector xyz){ pos = xyz; };

      G4int GetTrackID()    { return trackID; };
      G4double GetEdep()    { return edep; };
      G4ThreeVector GetPos(){ return pos; };

  private:

      G4int         trackID;
      G4double      edep;
      G4ThreeVector pos;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<HPGeDetectorHit> HPGeDetectorHitsCollection;

extern G4Allocator<HPGeDetectorHit> HPGeDetectorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* HPGeDetectorHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) HPGeDetectorHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void HPGeDetectorHit::operator delete(void *aHit)
{
  HPGeDetectorHitAllocator.FreeSingle((HPGeDetectorHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
