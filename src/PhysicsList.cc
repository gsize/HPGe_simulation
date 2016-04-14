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
/// \file electromagnetic/TestEm0/src/PhysicsList.cc
/// \brief Implementation of the PhysicsList class
//
// 
// $Id: PhysicsList.cc 66241 2012-12-13 18:34:42Z gunter $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"

#include "PhysListEmStandard.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"

#include "G4LossTableManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList() 
	: G4VModularPhysicsList(),fCutForGamma(0),fCutForElectron(0),fCutForPositron(0),
	fCurrentDefaultCut(0),fEmPhysicsList(0),fEmName("local"),fMessenger(0)
{    
	G4LossTableManager::Instance();

	fCurrentDefaultCut   = 0.01*mm;
	fCutForGamma         = fCurrentDefaultCut;
	fCutForElectron      = fCurrentDefaultCut;
	fCutForPositron      = fCurrentDefaultCut;

	fMessenger = new PhysicsListMessenger(this);

	SetVerboseLevel(1);

	// EM physics
	fEmName = G4String("local");
	fEmPhysicsList = new PhysListEmStandard(fEmName);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{
	delete fMessenger;
	delete fEmPhysicsList;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Bosons
#include "G4ChargedGeantino.hh"
#include "G4Geantino.hh"
#include "G4Gamma.hh"
#include "G4GenericIon.hh"
#include "G4BosonConstructor.hh"  
#include "G4LeptonConstructor.hh" 
#include "G4MesonConstructor.hh" 
#include "G4BaryonConstructor.hh" 
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
	// pseudo-particles
	G4Geantino::GeantinoDefinition();
	G4ChargedGeantino::ChargedGeantinoDefinition();

	// gamma
	G4Gamma::GammaDefinition();

	//bosons 
	
	G4BosonConstructor  pBosonConstructor;
	pBosonConstructor.ConstructParticle();

	// leptons
	G4LeptonConstructor pLeptonConstructor;
	pLeptonConstructor.ConstructParticle();

	// mesons
	G4MesonConstructor pMesonConstructor;
	pMesonConstructor.ConstructParticle();

	// barions
	G4BaryonConstructor pBaryonConstructor;
	pBaryonConstructor.ConstructParticle();

	// ions
	G4IonConstructor pIonConstructor;
	pIonConstructor.ConstructParticle();

	G4ShortLivedConstructor pShortLivedConstructor;
	pShortLivedConstructor.ConstructParticle();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4EmProcessOptions.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"

void PhysicsList::ConstructProcess()
{
	// Transportation
	//
	AddTransportation();

	// Electromagnetic physics list
	//
	fEmPhysicsList->ConstructProcess();

	// Decay
	RegisterPhysics(new G4DecayPhysics());
	//
	//     // Radioactive decay
	RegisterPhysics(new G4RadioactiveDecayPhysics());
/*
	// Em options
	//
	// Main options and setting parameters are shown here.
	// Several of them have default values.
	//
	G4EmProcessOptions emOptions;

	//physics tables
	//
	//emOptions.SetMinEnergy(100*eV);        //default    
	//emOptions.SetMaxEnergy(100*TeV);        //default  
	//emOptions.SetDEDXBinning(12*20);        //default=12*7  
	//emOptions.SetLambdaBinning(12*20);        //default=12*7

	emOptions.SetBuildCSDARange(true);     
	emOptions.SetMaxEnergyForCSDARange(10*GeV);
	//emOptions.SetDEDXBinningForCSDARange(12*20);

	//emOptions.SetSplineFlag(true);        //default

	emOptions.SetVerbose(0);  
*/
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "G4PhysListFactory.hh"
void PhysicsList::AddPackage(const G4String& name)
{
	G4PhysListFactory factory;
	G4VModularPhysicsList* phys =factory.GetReferencePhysList(name);
	G4int i=0;
	const G4VPhysicsConstructor* elem= phys->GetPhysics(i);
	G4VPhysicsConstructor* tmp = const_cast<G4VPhysicsConstructor*> (elem);
	while (elem !=0)
	{
		RegisterPhysics(tmp);
		elem= phys->GetPhysics(++i) ;
		tmp = const_cast<G4VPhysicsConstructor*> (elem);
	}
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddPhysicsList(const G4String& name)
{
	if (verboseLevel>0) {
		G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
	}

	if (name == fEmName) return;

	if (name == "local") {

		fEmName = name;
		delete fEmPhysicsList;
		fEmPhysicsList = new PhysListEmStandard(name);

	} else if (name == "emstandard_opt0"){
		fEmName = name;
		delete fEmPhysicsList;
		fEmPhysicsList = new G4EmStandardPhysics();

	} else if (name == "emstandard_opt1"){
		fEmName = name;
		delete fEmPhysicsList;
		fEmPhysicsList = new G4EmStandardPhysics_option1();

	} else if (name == "emstandard_opt2"){
		fEmName = name;
		delete fEmPhysicsList;
		fEmPhysicsList = new G4EmStandardPhysics_option2();

	} else if (name == "emstandard_opt3"){
		fEmName = name;
		delete fEmPhysicsList;
		fEmPhysicsList = new G4EmStandardPhysics_option3();

	} else if (name == "emstandard_opt4"){
		fEmName = name;
		delete fEmPhysicsList;
		fEmPhysicsList = new G4EmStandardPhysics_option4();

	} else if (name == "empenelope"){
		fEmName = name;
		delete fEmPhysicsList;
		fEmPhysicsList = new G4EmPenelopePhysics();

	} else if (name == "emlivermore"){
		fEmName = name;
		delete fEmPhysicsList;
		fEmPhysicsList = new G4EmLivermorePhysics();

	} else {

		G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
			<< " is not defined"
			<< G4endl;
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

void PhysicsList::SetCuts()
{ 
/*	// fixe lower limit for cut
	G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(100*eV, 1*GeV);

	// set cut values for gamma at first and for e- second and next for e+,
	// because some processes for e+/e- need cut values for gamma
	SetCutValue(fCutForGamma, "gamma");
	SetCutValue(fCutForElectron, "e-");
	SetCutValue(fCutForPositron, "e+");
	DumpCutValuesTable();
	*/
	SetCutsWithDefault();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForGamma(G4double cut)
{
	fCutForGamma = cut;
	SetParticleCuts(fCutForGamma, G4Gamma::Gamma());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForElectron(G4double cut)
{
	fCutForElectron = cut;
	SetParticleCuts(fCutForElectron, G4Electron::Electron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForPositron(G4double cut)
{
	fCutForPositron = cut;
	SetParticleCuts(fCutForPositron, G4Positron::Positron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
