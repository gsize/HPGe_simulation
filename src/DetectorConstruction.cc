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
// $Id:   DetectorConstruction.cc,v 1.22 2010-01-22 11:57:03 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
//#include "SD.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4SDManager.hh"
//#include "G4SDChargedFilter.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"

#include "G4NistManager.hh"
//#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal 
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = 0; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

	DetectorConstruction::  DetectorConstruction()
	:solidWorld(0)
	,  logicWorld(0)
	,  physiWorld(0)
	 ,   stepLimit(0)
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~  DetectorConstruction()
{
	delete stepLimit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

G4VPhysicalVolume*   DetectorConstruction::Construct()
{
	// Clean old geometry, if any
	G4GeometryManager::GetInstance()->OpenGeometry();
	G4PhysicalVolumeStore::GetInstance()->Clean();
	G4LogicalVolumeStore::GetInstance()->Clean();
	G4SolidStore::GetInstance()->Clean();

	//--------- Material definition ---------
	DefineMaterials();
	ConstructWorld();
	ConstructPbShield();
	ConstructHPGeDetector();

	return physiWorld;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* DetectorConstruction::ConstructWorld()
{
	G4NistManager* nist = G4NistManager::Instance();
	G4Material* Air = nist->FindOrBuildMaterial("G4_AIR");

	//--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------
	//------------------------------
	// World
	//------------------------------
	G4double fWorldLength = 640.*mm;
	G4double HalfWorldLength = 0.5*fWorldLength;

	G4GeometryManager::GetInstance()->SetWorldMaximumExtent(fWorldLength);
	G4cout << "Computed tolerance = "
		<< G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm
		<< " mm" << G4endl;

	solidWorld= new G4Box("world",HalfWorldLength,HalfWorldLength,HalfWorldLength);
	logicWorld= new G4LogicalVolume( solidWorld, Air, "World", 0, 0, 0);

	//  Must place the World Physical volume unrotated at (0,0,0).
	//
	physiWorld = new G4PVPlacement(0,               // no rotation
			G4ThreeVector(), // at (0,0,0)
			logicWorld,      // its logical volume
			"World",         // its name
			0,               // its mother  volume
			false,           // no boolean operations
			0);              // copy number
	return physiWorld;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume*   DetectorConstruction::ConstructPbShield()
{
	G4NistManager* nist = G4NistManager::Instance();

	G4double Shield_Tubs_rmin = 0.*mm;
	G4double Shield_Length = 630.*mm;
	G4double Shield_rMax = 0.5*510.*mm;
	G4double  Shield_Tubs_sphi =   0.*deg;
	G4double  Shield_Tubs_dphi = 360.*deg;

	G4double Shield_Fe_Len = 9.5*mm;
	G4double Shield_Pb_Len = 101.*mm;
	G4double Shield_Cu_Len = 1.6 *mm;
	G4double Shield_Sn_Len = 0.5*mm;

	//Shield Fe
	G4double Shield_Fe_Tubs_rMax = Shield_rMax;
	G4double Shield_Fe_Tubs_dz = 0.5*Shield_Length;
	Shield_Fe = nist->FindOrBuildMaterial("G4_Fe");

	G4VSolid * Shield_Fe_tubs
		= new G4Tubs("Shield_Fe_Tubs_tubs",Shield_Tubs_rmin,Shield_Fe_Tubs_rMax,Shield_Fe_Tubs_dz,
				Shield_Tubs_sphi,Shield_Tubs_dphi);
	G4LogicalVolume * Shield_Fe_log
		= new G4LogicalVolume(Shield_Fe_tubs,Shield_Fe,"Shield_Fe_log",0,0,0);

	G4VPhysicalVolume * physiWorldPbShield =
		new G4PVPlacement(0,G4ThreeVector(),Shield_Fe_log,"Shield_Fe_phys",
				logicWorld,false,0);
	G4VisAttributes* Shield_Fe_logVisAtt
		= new G4VisAttributes(G4Colour(1.0,0.0,1.0));
	Shield_Fe_log->SetVisAttributes(Shield_Fe_logVisAtt);

	//Shield Pb
	G4double Shield_Pb_Tubs_rMax = Shield_Fe_Tubs_rMax - Shield_Fe_Len;
	G4double Shield_Pb_Tubs_dz = Shield_Fe_Tubs_dz - Shield_Fe_Len;
	Shield_Pb = nist->FindOrBuildMaterial("G4_Pb");

	G4VSolid * Shield_Pb_tubs
		= new G4Tubs("PbTubs_tubs",Shield_Tubs_rmin,Shield_Pb_Tubs_rMax,Shield_Pb_Tubs_dz,
				Shield_Tubs_sphi,Shield_Tubs_dphi);
	G4LogicalVolume * Shield_Pb_log
		= new G4LogicalVolume(Shield_Pb_tubs,Shield_Pb,"Shield_Pb_log",0,0,0);
	// G4VPhysicalVolume * tracker_phys =
	new G4PVPlacement(0,G4ThreeVector(),Shield_Pb_log,"Shield Pb phys",
			Shield_Fe_log,false,0);
	G4VisAttributes* Shield_Pb_logVisAtt
		= new G4VisAttributes(G4Colour(1.0,0.0,0.50));
	Shield_Pb_log->SetVisAttributes(Shield_Pb_logVisAtt);

	//Shield Sn
	G4double Shield_Sn_Tubs_rMax = Shield_Pb_Tubs_rMax - Shield_Pb_Len;
	G4double Shield_Sn_Tubs_dz = Shield_Pb_Tubs_dz - Shield_Pb_Len;
	Shield_Sn = nist->FindOrBuildMaterial("G4_Sn");

	G4VSolid * Shield_Sn_tubs
		= new G4Tubs("SnTubs_tubs",Shield_Tubs_rmin,Shield_Sn_Tubs_rMax,Shield_Sn_Tubs_dz,
				Shield_Tubs_sphi,Shield_Tubs_dphi);
	G4LogicalVolume * Shield_Sn_log
		= new G4LogicalVolume(Shield_Sn_tubs,Shield_Sn,"Shield Sn",0,0,0);
	// G4VPhysicalVolume * tracker_phys =
	new G4PVPlacement(0,G4ThreeVector(),Shield_Sn_log,"Shield Sn phys",
			Shield_Pb_log,false,0);
	//  G4VisAttributes* Shield_Sn_logVisAtt
	//    = new G4VisAttributes(G4Colour(1.0,0.0,0.30));
	//  Shield_Sn_log->SetVisAttributes(Shield_Sn_logVisAtt);


	//Shield Cu
	G4double Shield_Cu_Tubs_rMax = Shield_Sn_Tubs_rMax - Shield_Sn_Len;
	G4double Shield_Cu_Tubs_dz = Shield_Sn_Tubs_dz - Shield_Sn_Len;
	Shield_Cu = nist->FindOrBuildMaterial("G4_Cu");

	G4VSolid * Shield_Cu_tubs
		= new G4Tubs("CuTubs_tubs",Shield_Tubs_rmin,Shield_Cu_Tubs_rMax,Shield_Cu_Tubs_dz,
				Shield_Tubs_sphi,Shield_Tubs_dphi);
	G4LogicalVolume * Shield_Cu_log
		= new G4LogicalVolume(Shield_Cu_tubs,Shield_Cu,"Shield Cu",0,0,0);
	// G4VPhysicalVolume * tracker_phys =
	new G4PVPlacement(0,G4ThreeVector(),Shield_Cu_log,"Shield Cu phys",
			Shield_Sn_log,false,0);
	//  G4VisAttributes* Shield_Cu_logVisAtt
	//    = new G4VisAttributes(G4Colour(1.0,0.50,0.60));
	//  Shield_Cu_log->SetVisAttributes(Shield_Cu_logVisAtt);

	//Shield Air
	G4double Shield_Air_Tubs_rMax = Shield_Cu_Tubs_rMax - Shield_Cu_Len;
	G4double Shield_Air_Tubs_dz = Shield_Cu_Tubs_dz - Shield_Cu_Len;
	Shield_Air = nist->FindOrBuildMaterial("G4_AIR");

	G4VSolid * Shield_Air_tubs
		= new G4Tubs("AirTubs_tubs",Shield_Tubs_rmin,Shield_Air_Tubs_rMax,Shield_Air_Tubs_dz,
				Shield_Tubs_sphi,Shield_Tubs_dphi);
	G4LogicalVolume * Shield_Air_log
		= new G4LogicalVolume(Shield_Air_tubs,Shield_Air,"Shield Air",0,0,0);
	// G4VPhysicalVolume * tracker_phys =
	new G4PVPlacement(0,G4ThreeVector(),Shield_Air_log,"Shield Air phys",
			Shield_Cu_log,false,0);
	G4VisAttributes* Shield_Air_logVisAtt
		= new G4VisAttributes(G4Colour(1.0,0.3,1.0));
	Shield_Air_log->SetVisAttributes(Shield_Air_logVisAtt);

	return physiWorldPbShield;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume*   DetectorConstruction::ConstructHPGeDetector()
{
	G4NistManager* nist = G4NistManager::Instance();

	G4double Tubs_rmin = 0.*mm;
	G4double  Tubs_sphi =   0.*deg;
	G4double  Tubs_dphi = 360.*deg;

	//shield Al
	G4double detector_move_len = 135 *mm;// 120. *mm;
	G4double Shell_Al_Tubs_rMax = 0.5 * 76. *mm;
	G4double Shell_Al_Tubs_dz = 0.5*110.*mm;
	G4double Shell_Al_Tubs_Thickness = 1.*mm;
	G4Material* Shell_Al = nist->FindOrBuildMaterial("G4_Al");

	G4VSolid * Shell_Al_tubs
		= new G4Tubs("Al_Tubs_tubs",Tubs_rmin,Shell_Al_Tubs_rMax,Shell_Al_Tubs_dz,
				Tubs_sphi,Tubs_dphi);
	G4LogicalVolume * Shell_Al_log
		= new G4LogicalVolume(Shell_Al_tubs,Shell_Al,"Shield Al",0,0,0);
	G4VPhysicalVolume * physiWorldHPGe =
		new G4PVPlacement(0,G4ThreeVector(0.,0.,detector_move_len),Shell_Al_log,"Detector Shell Al phys",
				logicWorld,false,0);
	G4VisAttributes* Shell_Al_logVisAtt
		= new G4VisAttributes(G4Colour(0.40,0.0,0.70));
	Shell_Al_log->SetVisAttributes(Shell_Al_logVisAtt);

	//shield Galactic
	G4double Shell_Galactic_Tubs_rMax = Shell_Al_Tubs_rMax - Shell_Al_Tubs_Thickness;
	G4double Shell_Galactic_Tubs_dz = Shell_Al_Tubs_dz - Shell_Al_Tubs_Thickness;
	G4Material* Shell_Galactic = nist->FindOrBuildMaterial("G4_Galactic");

	G4VSolid * Shell_Galactic_tubs
		= new G4Tubs("Galactic_Tubs_tubs",Tubs_rmin,Shell_Galactic_Tubs_rMax,Shell_Galactic_Tubs_dz,
				Tubs_sphi,Tubs_dphi);
	G4LogicalVolume * Shell_Galactic_log
		= new G4LogicalVolume(Shell_Galactic_tubs,Shell_Galactic,"Shell Galactic logic",0,0,0);
	// G4VPhysicalVolume * tracker_phys =
	new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),Shell_Galactic_log,"Shell Galactic phys",
			Shell_Al_log,false,0);
	//  G4VisAttributes* Shell_Galactic_logVisAtt
	//    = new G4VisAttributes(G4Colour(0.3,0.7,1.0));
	//  Shell_Galactic_log->SetVisAttributes(Shell_Galactic_logVisAtt);

	//HPGe dead_layer
	G4double HPGe_dead_layer_rMax = 0.5*63.0*mm;
	G4double HPGe_dead_layer_Tubs_dz = 0.5*44.10*mm;
	HPGe_detector_Ge = nist->FindOrBuildMaterial("G4_Ge");
	G4double HPGe_move = -( Shell_Galactic_Tubs_dz- HPGe_dead_layer_Tubs_dz - 4.*mm);//-27.95 *mm;
	G4double HPGe_dead_layer_Thickness = 0.70*mm;
	G4VSolid * HPGe_dead_layer_tubs
		= new G4Tubs("HPGe_dead_layer_tubs",Tubs_rmin,HPGe_dead_layer_rMax,HPGe_dead_layer_Tubs_dz,
				Tubs_sphi,Tubs_dphi);
	G4LogicalVolume * HPGe_dead_layer_log
		= new G4LogicalVolume(HPGe_dead_layer_tubs,HPGe_detector_Ge,"HPGe_dead_layer",0,0,0);
	// G4VPhysicalVolume * tracker_phys =
	new G4PVPlacement(0,G4ThreeVector(0. ,0. ,HPGe_move),HPGe_dead_layer_log,"HPGe_dead_layer_phys",
			Shell_Galactic_log,false,0);
	G4VisAttributes* HPGe_dead_layer_logVisAtt
		= new G4VisAttributes(G4Colour(8.0,6.0,1.20));
	HPGe_dead_layer_log->SetVisAttributes(HPGe_dead_layer_logVisAtt);

	//HPGe detector
	G4double HPGe_detector_rMax = HPGe_dead_layer_rMax - HPGe_dead_layer_Thickness;
	G4double HPGe_detector_Tubs_dz = HPGe_dead_layer_Tubs_dz - HPGe_dead_layer_Thickness;
	//HPGe_detector_Ge = nist->FindOrBuildMaterial("G4_Ge");
	//G4double HPGe_move = -30.45 *mm;

	G4VSolid * HPGe_detector_tubs
		= new G4Tubs("HPGe_detector_tubs",Tubs_rmin,HPGe_detector_rMax,HPGe_detector_Tubs_dz,
				Tubs_sphi,Tubs_dphi);
	G4LogicalVolume * HPGe_detector_log
		= new G4LogicalVolume(HPGe_detector_tubs,HPGe_detector_Ge,"HPGe_detector",0,0,0);
	// G4VPhysicalVolume * tracker_phys =
	new G4PVPlacement(0,G4ThreeVector(0. ,0. ,0.),HPGe_detector_log,"HPGe_detector_phys",
			HPGe_dead_layer_log,false,0);
	G4VisAttributes* HPGe_detector_logVisAtt
		= new G4VisAttributes(G4Colour(8.0,5.0,1.20));
	HPGe_detector_logVisAtt->G4VisAttributes::SetForceSolid(true);
	delete  HPGe_detector_log->GetVisAttributes();
	HPGe_detector_log->SetVisAttributes(HPGe_detector_logVisAtt);

	//HPGe detector inner dead_layer
	G4double HPGe_inner_dead_layer_thickness = 0.3*um;
	G4double HPGe_inner_dead_layer_rMax = 0.5 * 10.8*mm + HPGe_inner_dead_layer_thickness;
	G4double HPGe_inner_dead_layer_Tubs_dz = 0.5*27.80*mm + HPGe_inner_dead_layer_thickness;
	G4double HPGe_inner_dead_layer_move = HPGe_dead_layer_Tubs_dz - HPGe_inner_dead_layer_Tubs_dz;

	G4VSolid * HPGe_inner_dead_layer_tubs
		= new G4Tubs("HPGe_inner_dead_layer_tubs",Tubs_rmin,HPGe_inner_dead_layer_rMax,HPGe_inner_dead_layer_Tubs_dz,
				Tubs_sphi,Tubs_dphi);
	G4LogicalVolume * HPGe_inner_dead_layer_log
		= new G4LogicalVolume(HPGe_inner_dead_layer_tubs,HPGe_detector_Ge,"HPGe_inner_dead_layer_log",0,0,0);
	// G4VPhysicalVolume * tracker_phys =
	new G4PVPlacement(0,G4ThreeVector(0. ,0. ,HPGe_inner_dead_layer_move),HPGe_inner_dead_layer_log,"HPGe_inner_dead_layer_phys",
			HPGe_dead_layer_log,false,0);
	G4VisAttributes* HPGe_inner_dead_layer_logVisAtt
		= new G4VisAttributes(G4Colour(3.0,6.0,1.20));
	HPGe_inner_dead_layer_log->SetVisAttributes(HPGe_inner_dead_layer_logVisAtt);

	//HPGe detector inner Cu
	G4double HPGe_inner_tubs_Cu_rMax = HPGe_inner_dead_layer_rMax -  HPGe_inner_dead_layer_thickness;
	G4double HPGe_inner_tubs_Cu_dz = HPGe_inner_dead_layer_Tubs_dz -  HPGe_inner_dead_layer_thickness;
	G4double HPGe_inner_tubs_Cu_move = HPGe_dead_layer_Tubs_dz - HPGe_inner_dead_layer_Tubs_dz;
	G4Material* Mat_Cu = nist->FindOrBuildMaterial("G4_Cu");

	G4VSolid * HPGe_inner_tubs_Cu
		= new G4Tubs("HPGe_inner_tubs_Cu_tubs",Tubs_rmin,HPGe_inner_tubs_Cu_rMax,HPGe_inner_tubs_Cu_dz,
				Tubs_sphi,Tubs_dphi);
	G4LogicalVolume * HPGe_inner_tubs_Cu_log
		= new G4LogicalVolume(HPGe_inner_tubs_Cu,Mat_Cu,"HPGe_inner_tubs_Cu_tubs_log",0,0,0);
	// G4VPhysicalVolume * tracker_phys =
	new G4PVPlacement(0,G4ThreeVector(0. ,0. ,HPGe_inner_tubs_Cu_move),HPGe_inner_tubs_Cu_log,"HPGe_inner_tubs_Cu_tubs_phys",
			HPGe_dead_layer_log,false,0);
	G4VisAttributes* HPGe_inner_tubs_Cu_logVisAtt
		= new G4VisAttributes(G4Colour(4.0,5.0,3.20));
	HPGe_inner_tubs_Cu_logVisAtt->G4VisAttributes::SetForceSolid(true);
	delete  HPGe_inner_tubs_Cu_log->GetVisAttributes();
	HPGe_inner_tubs_Cu_log->SetVisAttributes(HPGe_inner_tubs_Cu_logVisAtt);

	//--------- example of User Limits -------------------------------

	// below is an example of how to set tracking constraints in a given
	// logical volume(see also in N02PhysicsList how to setup the processes
	// G4StepLimiter or G4UserSpecialCuts).

	// Sets a max Step length in the tracker region, with G4StepLimiter
	//
	//G4double maxStep = HPGe_detector_Tubs_dz;
	//stepLimit = new G4UserLimits(maxStep);
	//HPGe_detector_log->SetUserLimits(stepLimit);
	// Set additional contraints on the track, with G4UserSpecialCuts
	//
	// G4double maxLength = 2*fTrackerLength, maxTime = 0.1*ns, minEkin = 10*MeV;
	// logicTracker->SetUserLimits(new G4UserLimits(maxStep,maxLength,maxTime,
	//                                               minEkin));
	return physiWorldHPGe;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void   DetectorConstruction::DefineMaterials()
{
	// Lead material defined using NIST Manager
	G4NistManager* nistManager = G4NistManager::Instance();
	G4bool fromIsotopes = false;

	// material
	nistManager->FindOrBuildMaterial("G4_Fe", fromIsotopes);
	nistManager->FindOrBuildMaterial("G4_Cu", fromIsotopes);
	nistManager->FindOrBuildMaterial("G4_Pb", fromIsotopes);
	nistManager->FindOrBuildMaterial("G4_Sn", fromIsotopes);
	nistManager->FindOrBuildMaterial("G4_Ge", fromIsotopes);
	nistManager->FindOrBuildMaterial("G4_AIR");

	// Print materials
	G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
	G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
	// 
	// Scorers
	//

	// declare HPGe as a MultiFunctionalDetector scorer
	//  
	G4MultiFunctionalDetector* HPGeDetector 
		= new G4MultiFunctionalDetector("HPGe");

	G4VPrimitiveScorer* primitive;
	primitive = new G4PSEnergyDeposit("Edep");
	HPGeDetector->RegisterPrimitive(primitive);
	/*
	   primitive = new G4PSTrackLength("TrackLength");
	   G4SDChargedFilter* charged = new G4SDChargedFilter("chargedFilter");
	   primitive ->SetFilter(charged);
	   absDetector->RegisterPrimitive(primitive);  
	   */
	SetSensitiveDetector("HPGe_detector",HPGeDetector);

	// 
	// Magnetic field
	//
	// Create global magnetic field messenger.
	// Uniform magnetic field is then created automatically if
	// the field value is not zero.
	G4ThreeVector fieldValue = G4ThreeVector();
	fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
	fMagFieldMessenger->SetVerboseLevel(1);

	// Register the field messenger for deleting
	G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
