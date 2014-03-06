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
// $Id:   DetectorConstruction.hh,v 1.10 2008-09-22 16:41:20 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef   DetectorConstruction_h
#define   DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

class G4GlobalMagFieldMessenger;

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4VPVParameterisation;
class G4UserLimits;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class   DetectorConstruction : public G4VUserDetectorConstruction
{
	public:

		DetectorConstruction();
		~  DetectorConstruction();

	public:

		G4VPhysicalVolume* Construct();

		virtual void ConstructSDandField();

	private:
		void DefineMaterials();
		void ConstructWorld();
		void ConstructPbShield();
		void ConstructHPGeDetector();
	private:
		// data members
		//
		G4Box*             solidWorld;    // pointer to the solid envelope
		G4LogicalVolume*   logicWorld;    // pointer to the logical envelope
		G4VPhysicalVolume* physiWorld;    // pointer to the physical envelope
		G4UserLimits* stepLimit;             // pointer to user step limits
		G4Material* Shield_Fe;
		G4Material* Shield_Cu;
		G4Material* Shield_Sn;
		G4Material* Shield_Pb;
		G4Material* Shield_Air;
		G4Material* HPGe_detector_Ge;

		static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger; 
		// magnetic field messenger
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
