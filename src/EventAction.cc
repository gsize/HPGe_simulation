//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id:   EventAction.cc,v 1.11 2006-06-29 17:48:05 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"

#include "HistoManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

	EventAction::  EventAction(HistoManager *histo)
:G4UserEventAction(),
	fHistManager(histo),
	fHPGeEdepHCID(-1)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~  EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void   EventAction::BeginOfEventAction(const G4Event* /*event*/)
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "Analysis.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4RunManager.hh"

void   EventAction::EndOfEventAction(const G4Event* event)
{
//	G4int event_id = event->GetEventID();
// Get the energy of primary particle 
	G4double primary_energy = event->GetPrimaryVertex()->GetPrimary()->GetTotalEnergy();
	fHistManager->FillSourceData(primary_energy /MeV);

	// Get hist collections IDs
	G4double energyShield = 1. *keV;
	if ( fHPGeEdepHCID < 0  ) {
		fHPGeEdepHCID 
			= G4SDManager::GetSDMpointer()->GetCollectionID("HPGe/Edep");
		}
		// Get sum values from hits collections
		//
		G4double HPGeEdep = GetSum(GetHitsCollection(fHPGeEdepHCID, event));

		// fill histograms
		//
		if(HPGeEdep > energyShield) {
		fHistManager->FillSDData(HPGeEdep/ MeV);
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4THitsMap<G4double>* 
EventAction::GetHitsCollection(G4int hcID,
		const G4Event* event) const
{
	G4THitsMap<G4double>* hitsCollection 
		= static_cast<G4THitsMap<G4double>*>(
				event->GetHCofThisEvent()->GetHC(hcID));

	if ( ! hitsCollection ) {
		G4ExceptionDescription msg;
		msg << "Cannot access hitsCollection ID " << hcID; 
		G4Exception("B4dEventAction::GetHitsCollection()",
				"MyCode0003", FatalException, msg);
	}         

	return hitsCollection;
}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double EventAction::GetSum(G4THitsMap<G4double>* hitsMap) const
{
	G4double sumValue = 0;
	std::map<G4int, G4double*>::iterator it;
	for ( it = hitsMap->GetMap()->begin(); it != hitsMap->GetMap()->end(); it++) {
		sumValue += *(it->second);
	}
	return sumValue;  
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventAction::PrintEventStatistics(
		G4double absoEdep) const
{
	// Print event statistics
	//
	G4cout
		<< "   Absorber: total energy: " 
		<< std::setw(7) << G4BestUnit(absoEdep, "Energy")
		/*     << "       total track length: " 
			   << std::setw(7) << G4BestUnit(absoTrackLength, "Length")
			   << G4endl
			   << "        Gap: total energy: " 
			   << std::setw(7) << G4BestUnit(gapEdep, "Energy")
			   << "       total track length: " 
			   << std::setw(7) << G4BestUnit(gapTrackLength, "Length")
			   */
		<< G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
