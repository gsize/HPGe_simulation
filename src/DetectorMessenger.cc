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
/// \file hadronic/Hadr03/src/DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
// $Id: DetectorMessenger.cc 70755 2013-06-05 12:17:48Z ihrivnac $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:G4UImessenger(), 
 fDetector(Det), fDetDir(0), fOutDeadLayerThicknessCmd(0)
	,fFlagPbShieldCmd(0)
{ 
  
  G4bool broadcast = false;
  fDetDir = new G4UIdirectory("/HPGe_simulation/det/",broadcast);
  fDetDir->SetGuidance("detector construction commands");
        
  fOutDeadLayerThicknessCmd = new G4UIcmdWithADoubleAndUnit("/HPGe_simulation/det/setOutDeadLayerThickness",this);
  fOutDeadLayerThicknessCmd->SetGuidance("Set out dead layer thickness of HPGe .");
  fOutDeadLayerThicknessCmd->SetParameterName("outDeadLayerThickness",false);
  fOutDeadLayerThicknessCmd->SetDefaultUnit("mm");
  fOutDeadLayerThicknessCmd->SetUnitCategory("Length");
  fOutDeadLayerThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  fFlagPbShieldCmd = new G4UIcmdWithABool("/HPGe_simulation/det/setPbShield",this);
  fFlagPbShieldCmd->SetGuidance("add Pb Shield .");
  fFlagPbShieldCmd->SetParameterName("flagPbShield",false);
  fFlagPbShieldCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fOutDeadLayerThicknessCmd;
  delete fDetDir;
  delete fFlagPbShieldCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
   
  if( command == fOutDeadLayerThicknessCmd )
   { fDetector->SetOutDeadLayerThickness(fOutDeadLayerThicknessCmd->GetNewDoubleValue(newValue));}
   if(command == fFlagPbShieldCmd) 
   {
	   fDetector->SetPbShield(fFlagPbShieldCmd->GetNewBoolValue(newValue));
   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
