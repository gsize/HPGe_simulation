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
/// \file analysis/AnaEx01/src/HistoManager.cc
/// \brief Implementation of the HistoManager class
//
//
// $Id$
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "HistoManager.hh"
//#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
	:fileName("HPGe_data")
	 ,factoryOn(false)
{
	// histograms
	for (G4int k=0; k<MaxHisto; k++) {
		fHistId[k] = 0;
		fHistPt[k] = 0;    
	}
	// ntuple
	for (G4int k=0; k<MaxNtCol; k++) {
		fNtColId[k] = 0;
	}  

	Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{ 
	delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
	// Create or get analysis manager
	// The choice of analysis technology is done via selection of a namespace
	// in HistoManager.hh
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	analysisManager->SetVerboseLevel(1);
	//analysisManager->SetActivation(true);   //enable inactivation of histograms
	G4cout << "Using"<< analysisManager->GetFileType() << G4endl;
	analysisManager->SetFileName(fileName);

	// Create directories 
	analysisManager->SetHistoDirectoryName("histo");
	analysisManager->SetNtupleDirectoryName("ntuple");

	// create selected histograms
	//
	analysisManager->SetFirstHistoId(1);
	analysisManager->SetFirstNtupleId(1);
	fHistId[0] = analysisManager->CreateH1("source","Gamma source (MeV)",
			8192, 0., 2.0*MeV);
	fHistPt[0] = analysisManager->GetH1(fHistId[0]);

	fHistId[1] = analysisManager->CreateH1("HPGe","Gamma HPGe (MeV)",
			8192, 0., 2.0*MeV);
	fHistPt[1] = analysisManager->GetH1(fHistId[1]);

	// Creating ntuple ID=1
	//
	analysisManager->CreateNtuple("101", "source");
	fNtColId[0] = analysisManager->CreateNtupleDColumn("Edep_source");
	analysisManager->FinishNtuple();
	// Creating ntuple ID=2
	analysisManager->CreateNtuple("102", "HPGe");
	fNtColId[1] = analysisManager->CreateNtupleDColumn("Edep_HPGe");
	analysisManager->FinishNtuple();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void HistoManager::OpenFile()
{
	// Create or get analysis manager
	// The choice of analysis technology is done via selection of a namespace
	// in HistoManager.hh
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	fileName = 	analysisManager->GetFileName();
	// Open an output file
	//
	//	if ( analysisManager->IsActive() ) {
	G4bool fileOpen = analysisManager->OpenFile();
	if (!fileOpen) {
		G4cout << "\n---> HistoManager::OpenFile(): cannot open " << fileName 
			<< G4endl;
		return;
	}
	//	}

	factoryOn = true;       
	G4cout << "\n----> Histogram Tree is opened in " << fileName << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Save()
{
	if (factoryOn) {
		G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
		//	if ( analysisManager->IsActive() ) {    
		analysisManager->Write();
		analysisManager->CloseFile();  
		G4cout << "\n----> Histogram Tree is saved in " << fileName << G4endl;
		//	}
		factoryOn = false;
	}                    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Normalize(G4int ih, G4double fac)
{
	if (ih >= MaxHisto) {
		G4cout << "---> warning from HistoManager::Normalize() : histo " << ih
			<< "  fac= " << fac << G4endl;
		return;
	}

	//  if (fHistPt[ih]) fHistPt[ih]->scale(fac);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::FillHisto(G4int ih, G4double e, G4double weight)
{
	if (ih > MaxHisto) {
		G4cout << "---> warning from HistoManager::FillHisto() : histo " << ih
			<< "does note xist; xbin= " << e << " w= " << weight << G4endl;
		return;
	}

	//	G4cout<<ih<<" histo:"<<fHistId[ih]<<" engry:"<< e <<G4endl;
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	if (analysisManager) analysisManager->FillH1(ih, e, weight);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::FillNtuple(int ih ,int NtColID,G4double energy)
{                
	//	G4cout<<ih<<" Ntuple:"<<fNtColId[NtColID]<<G4endl;
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	analysisManager->FillNtupleDColumn(ih,fNtColId[NtColID], energy);
	analysisManager->AddNtupleRow(ih);  
}  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::FillSourceData(G4double energy)
{
	this -> FillHisto(1,energy);
	this -> FillNtuple(1,0,energy);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::FillSDData(G4double energy)
{
	this -> FillHisto(2,energy);
	this -> FillNtuple(2,1,energy);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::PrintStatistic()
{
	/*
	   if(factoryOn) {
	   G4cout << "\n ----> print histograms statistic \n" << G4endl;

	   G4cout 
	   << " EAbs : mean = " << G4BestUnit(fHistPt[1]->mean(), "Energy") 
	   << " rms = " << G4BestUnit(fHistPt[1]->rms(),  "Energy") 
	   << G4endl;
	   G4cout                
	   << " EGap : mean = " << G4BestUnit(fHistPt[2]->mean(), "Energy") 
	   << " rms = " << G4BestUnit(fHistPt[2]->rms(),  "Energy") 
	   << G4endl;
	   G4cout 
	   << " LAbs : mean = " << G4BestUnit(fHistPt[3]->mean(), "Length") 
	   << " rms = " << G4BestUnit(fHistPt[3]->rms(),  "Length") 
	   << G4endl;
	   G4cout 
	   << " LGap : mean = " << G4BestUnit(fHistPt[4]->mean(), "Length") 
	   << " rms = " << G4BestUnit(fHistPt[4]->rms(),  "Length") 
	   << G4endl;
	   }
	   */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


