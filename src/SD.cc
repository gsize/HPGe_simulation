//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4ParticleChange.hh"
#include "G4ParticleChangeForMSC.hh"

#include "G4VProcess.hh"

#include "Analysis.hh"

//extern TTree mtree;
//extern test_block test;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HPGeDetectorSD::HPGeDetectorSD(G4String name)
    :G4VSensitiveDetector(name)
    ,e_g(0.0)
    ,l_g(0.0)
{
    G4String HCname;
    collectionName.insert(HCname="trackerCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HPGeDetectorSD::~HPGeDetectorSD() { }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HPGeDetectorSD::Initialize(G4HCofThisEvent* HCE)
{
    trackerCollection = new HPGeDetectorHitsCollection
    (SensitiveDetectorName,collectionName[0]);
    static G4int HCID = -1;
    if(HCID<0)
    {
        HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    }
    HCE->AddHitsCollection( HCID, trackerCollection );
    //Intitalize the total energy and range
    e_g  = 0.0;
    l_g  = 0.0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool HPGeDetectorSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
//Fetch the energy deposit
    G4double edep = aStep->GetTotalEnergyDeposit();
//Fetch the step length (range)
    G4double ltra = aStep->GetStepLength();
//Fetch the kinetic energy
    G4double kne   = aStep->GetTrack()->GetKineticEnergy();
//Fetch the volume (sub detector) name
    G4String pname = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();
//Fetch the G4 track ID
    G4int    trackid = aStep->GetTrack()->GetTrackID();
//Fetch the corresponding particle name
    G4String p_name = aStep->GetTrack()->GetDefinition()->GetParticleName();
//Fetch the G4 parent particle name
//G4int   parentid = aStep->GetTrack()->GetParentID();
//Fetch the step position
//G4ThreeVector pos = aStep->GetPreStepPoint()->GetPosition();
//Fetch the momentum direction
    G4ThreeVector mom = aStep->GetPreStepPoint()->GetMomentumDirection();

    if(edep==0.) return false;

    e_g+=edep;
    l_g+=ltra;

    HPGeDetectorHit* newHit = new HPGeDetectorHit();
    newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());
    newHit->SetEdep     (edep);
    newHit->SetPos      (aStep->GetPostStepPoint()->GetPosition());
    trackerCollection->insert( newHit );

//newHit->Print();
//newHit->Draw();

    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HPGeDetectorSD::EndOfEvent(G4HCofThisEvent*)
{
    if (verboseLevel>0)
    {
        G4int NbHits = trackerCollection->entries();
        G4cout << "\n-------->Hits Collection: in this event they are " << NbHits
               << " hits in the tracker: " << G4endl;
        for (G4int i=0; i<NbHits; i++) (*trackerCollection)[i]->Print();
    }

        // get analysis manager
        G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

     if(e_g>6.*keV)
    {
        //fill hist
        analysisManager->FillH1(2, e_g);
        // fill ntuple
        //analysisManager->FillNtupleDColumn(2, e_g);
        //analysisManager->FillNtupleDColumn(1, l_g);

        //analysisManager->AddNtupleRow();
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

