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
// $Id: exampleN02.cc,v 1.16 2009-10-30 14:59:59 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4ScoringManager.hh"

//#include "FTFP_BERT.hh"
//#include "QGSP_BERT_HP.hh"
#include "PhysicsList.hh"

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

#include "G4UImanager.hh"

// MPI session
//#include "G4MPImanager.hh"
//#include "G4MPIsession.hh"


#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
    // random engine
   // CLHEP::Ranlux64Engine randomEngine;
   // CLHEP::HepRandom::setTheEngine(&randomEngine);
  //G4Random::setTheEngine(new CLHEP::RanecuEngine);
  
    // User Verbose output class
    //
    //G4VSteppingVerbose* verbosity = new  SteppingVerbose;
    //G4VSteppingVerbose::SetInstance(verbosity);
     
  // Construct the default run manager
  //
#ifdef G4MULTITHREADED
  G4int nThreads = 3;
  G4MTRunManager* runManager = new G4MTRunManager;
  
    // Number of threads can be defined via 3rd argument
  if (argc==3) {
    nThreads = G4UIcommand::ConvertToInt(argv[2]);
  }
  runManager->SetNumberOfThreads(nThreads);
  G4cout << "##### Application started for " << runManager->GetNumberOfThreads() 
         << " threads" << " #####" << G4endl;
#else
  G4RunManager* runManager = new G4RunManager;
#endif

// Activate UI-command base scorer
 G4ScoringManager * scManager = G4ScoringManager::GetScoringManager();
 scManager->SetVerboseLevel(1);
 
//====================================================================
// Un-comment this line for user defined score writer
//    scManager->SetScoreWriter(new RE03UserScoreWriter());
//====================================================================

    // User Initialization classes (mandatory)
    //
     DetectorConstruction* detector = new  DetectorConstruction;
    runManager->SetUserInitialization(detector);
    
      G4VModularPhysicsList* physicsList = new PhysicsList();/*QGSP_BERT_HP;//FTFP_BERT;*/
  runManager->SetUserInitialization(physicsList);

    ActionInitialization* actionInitialization
     = new ActionInitialization();
  runManager->SetUserInitialization(actionInitialization);

    // Initialize G4 kernel
    //
    //runManager->Initialize();

#ifdef G4VIS_USE
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();
#endif

// --------------------------------------------------------------------
    // ready for go
    // MPIsession treats both interactive and batch modes.
    // Just start your session as below.
    // --------------------------------------------------------------------
//    session-> SessionStart();
    
      // Get the pointer to the User Interface manager
      //
      G4UImanager * UImanager = G4UImanager::GetUIpointer();

      if (argc!=1)   // batch mode
        {
          G4String command = "/control/execute ";
          G4String fileName = argv[1];
          UImanager->ApplyCommand(command+fileName);
        }
      else           // interactive mode : define UI session
        {
    #ifdef G4UI_USE
          G4UIExecutive * ui = new G4UIExecutive(argc,argv);
    #ifdef G4VIS_USE
          UImanager->ApplyCommand("/control/execute vis.mac");
    #endif
          ui->SessionStart();
          delete ui;
    #endif
        }
    
#ifdef G4VIS_USE
    delete visManager;
#endif

    // Free the store: user actions, physics_list and detector_description are
    //                 owned and deleted by the run manager, so they should not
    //                 be deleted in the main() program !

    delete runManager;
    //delete verbosity;

    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

