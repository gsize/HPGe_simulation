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

class DetectorMessenger;
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
		void SetPbShield(G4bool value);

		void SetOutDeadLayerThickness(double value);
		void SetCoverThick(G4double value);
		void SetShellRadius(G4double value);
		void SetShellLength(G4double value);
		void SetShellThick(G4double value);
		void SetEndGap(G4double value);
		void SetCUPLength(G4double value);
		void SetCUPThick(G4double value);
		void SetCUPTopThick(G4double value);
		void SetCUPBottomThick(G4double value);
		void SetMylarThick (G4double value);
		void SetCrystalRadius (G4double value);
		void SetCrystalHalfLength (G4double value);
		void SetCrystalEndRadius (G4double value);
		void SetHoleDepth(G4double value);
		void SetHoleRadius(G4double value);
		void SetOuterDeadLayerThick(G4double value);
		void SetInnerDeadLayerThick(G4double value);
		void UpdateGeometry();

	private:
		void DefineMaterials();
		G4VPhysicalVolume* ConstructWorld();
		G4VPhysicalVolume* ConstructPbShield(G4LogicalVolume* );
		void ConstructHPGeDetector(G4LogicalVolume* );

	private:
		// data members
		//
		G4VPhysicalVolume* physiWorld;    // pointer to the physical envelope
		G4UserLimits* stepLimit;             // pointer to user step limits
		DetectorMessenger* detectorMessenger;	
		G4Material* Shield_Fe;
		G4Material* Shield_Cu;
		G4Material* Shield_Sn;
		G4Material* Shield_Pb;
		G4Material* Shield_Air;
		G4Material* GeCrystal;

		//cover
		G4double coverThick;

		//HPGe parameter
		//detector shell
		G4double shellRadius;
		G4double shellLength;
		G4double shellThick ;
		G4double endGap;
		G4double detectorMove;
		//CUP
		G4double CUPLength;
		G4double CUPThick ;
		G4double CUPTopThick ;
		G4double CUPBottomThick ;
		G4double mylarThick;

		//Ge crystal
		G4double crystalRadius;
		G4double crystalHalfLength;
		G4double crystalEndRadius;
		G4double holeDepth;
		G4double holeRadius;
		G4double outerDeadLayerThick;
		G4double innerDeadLayerThick ;

		G4bool fCheckOverlaps;
		G4bool flagPbShield;

		static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger; 
		// magnetic field messenger
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
