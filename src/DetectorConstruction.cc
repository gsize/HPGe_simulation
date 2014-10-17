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

#include "DetectorMessenger.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal 
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = 0; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

	DetectorConstruction::  DetectorConstruction()
	:G4VUserDetectorConstruction()
	,  physiWorld(0)
	,  stepLimit(0)
	,fCheckOverlaps(true)
	 ,flagPbShield(true)
{
	outDeadLayerThickness = 0.7 *mm;
	shellAlThickness = 1.0 *mm;
	detectorMessenger = new DetectorMessenger(this);	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~  DetectorConstruction()
{
	delete stepLimit;
	delete detectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

G4VPhysicalVolume* DetectorConstruction::Construct()
{
	// Clean old geometry, if any
	G4GeometryManager::GetInstance()->OpenGeometry();
	G4PhysicalVolumeStore::GetInstance()->Clean();
	G4LogicalVolumeStore::GetInstance()->Clean();
	G4SolidStore::GetInstance()->Clean();

	//--------- Material definition ---------
	DefineMaterials();
	G4LogicalVolume* pLV = 0;
	physiWorld = ConstructWorld();

	pLV = physiWorld->GetLogicalVolume();
	if(flagPbShield == true)
	{
		G4VPhysicalVolume* PVShieldDauther = ConstructPbShield(pLV);
		pLV = PVShieldDauther->GetLogicalVolume();
	}
	ConstructHPGeDetector(pLV);
	return physiWorld;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume*  DetectorConstruction::ConstructWorld()
{
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

	G4VSolid* solidWorld= new G4Box("world",HalfWorldLength,HalfWorldLength,HalfWorldLength);
	G4LogicalVolume* logicWorld= new G4LogicalVolume( solidWorld, Shield_Air, "LogicWorld", 0, 0, 0);

	//  Must place the World Physical volume unrotated at (0,0,0).
	//
	G4VPhysicalVolume* phyWorld = new G4PVPlacement(0,               // no rotation
			G4ThreeVector(), // at (0,0,0)
			logicWorld,      // its logical volume
			"PhysiWorld",         // its name
			0,               // its mother  volume
			false,           // no boolean operations
			0,              // copy number
			fCheckOverlaps);
	return phyWorld;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructPbShield(G4LogicalVolume* motherLogicalVolume)
{
	G4double Shield_Tubs_rmin = 0.*mm;
	G4double Shield_Length = 630.*mm;
	G4double Shield_rMax = 0.5*510.*mm;
	G4double  Shield_Tubs_sphi =   0.*deg;
	G4double  Shield_Tubs_dphi = 360.*deg;

	G4double Shield_Fe_Thickness = 9.5*mm;
	G4double Shield_Pb_Thickness = 101.*mm;
	G4double Shield_Cu_Thickness= 1.6 *mm;
	G4double Shield_Sn_Thickness= 0.5*mm;

	//Shield Fe
	G4double Shield_Fe_Tubs_rMax = Shield_rMax;
	G4double Shield_Fe_Tubs_dz = 0.5*Shield_Length;

	G4VSolid * Shield_Fe_tubs
		= new G4Tubs("Shield_Fe_tubs",Shield_Tubs_rmin,Shield_Fe_Tubs_rMax,Shield_Fe_Tubs_dz,
				Shield_Tubs_sphi,Shield_Tubs_dphi);
	G4LogicalVolume * Shield_Fe_log
		= new G4LogicalVolume(Shield_Fe_tubs,Shield_Fe,"Shield_Fe_log",0,0,0);

	//	G4VPhysicalVolume * physiWorldPbShield =
	new G4PVPlacement(0,G4ThreeVector(),Shield_Fe_log,"Shield_Fe_phys",
			motherLogicalVolume,false,0,fCheckOverlaps);

	G4VisAttributes* Shield_Fe_logVisAtt
		= new G4VisAttributes(G4Colour(1.0,0.0,1.0));
	Shield_Fe_log->SetVisAttributes(Shield_Fe_logVisAtt);

	//Shield Pb
	G4double Shield_Pb_Tubs_rMax = Shield_Fe_Tubs_rMax - Shield_Fe_Thickness;
	G4double Shield_Pb_Tubs_dz = Shield_Fe_Tubs_dz - Shield_Fe_Thickness;

	G4VSolid * Shield_Pb_tubs
		= new G4Tubs("Shield_Pb_tubs",Shield_Tubs_rmin,Shield_Pb_Tubs_rMax,Shield_Pb_Tubs_dz,
				Shield_Tubs_sphi,Shield_Tubs_dphi);
	G4LogicalVolume * Shield_Pb_log
		= new G4LogicalVolume(Shield_Pb_tubs,Shield_Pb,"Shield_Pb_log",0,0,0);
	// G4VPhysicalVolume * tracker_phys =
	new G4PVPlacement(0,G4ThreeVector(),Shield_Pb_log,"Shield_Pb_phys",
			Shield_Fe_log,false,0,fCheckOverlaps);
	G4VisAttributes* Shield_Pb_logVisAtt
		= new G4VisAttributes(G4Colour(1.0,0.0,0.50));
	Shield_Pb_log->SetVisAttributes(Shield_Pb_logVisAtt);


	//Shield Sn
	G4double Shield_Sn_Tubs_rMax = Shield_Pb_Tubs_rMax - Shield_Pb_Thickness;
	G4double Shield_Sn_Tubs_dz = Shield_Pb_Tubs_dz - Shield_Pb_Thickness;

	G4VSolid * Shield_Sn_tubs
		= new G4Tubs("Shield_Sn_tubs",Shield_Tubs_rmin,Shield_Sn_Tubs_rMax,Shield_Sn_Tubs_dz,
				Shield_Tubs_sphi,Shield_Tubs_dphi);
	G4LogicalVolume * Shield_Sn_log
		= new G4LogicalVolume(Shield_Sn_tubs,Shield_Sn,"Shield_Sn_log",0,0,0);
	// G4VPhysicalVolume * tracker_phys =
	new G4PVPlacement(0,G4ThreeVector(),Shield_Sn_log,"Shield_Sn_phys",
			Shield_Pb_log,false,0,fCheckOverlaps);
	G4VisAttributes* Shield_Sn_logVisAtt
		= new G4VisAttributes(G4Colour(1.0,0.0,0.30));
	Shield_Sn_log->SetVisAttributes(Shield_Sn_logVisAtt);


	//Shield Cu
	G4double Shield_Cu_Tubs_rMax = Shield_Sn_Tubs_rMax - Shield_Sn_Thickness;
	G4double Shield_Cu_Tubs_dz = Shield_Sn_Tubs_dz - Shield_Sn_Thickness;

	G4VSolid * Shield_Cu_tubs
		= new G4Tubs("Shield_Cu_tubs",Shield_Tubs_rmin,Shield_Cu_Tubs_rMax,Shield_Cu_Tubs_dz,
				Shield_Tubs_sphi,Shield_Tubs_dphi);
	G4LogicalVolume * Shield_Cu_log
		= new G4LogicalVolume(Shield_Cu_tubs,Shield_Cu,"Shield_Cu_log",0,0,0);
	// G4VPhysicalVolume * tracker_phys =
	new G4PVPlacement(0,G4ThreeVector(),Shield_Cu_log,"Shield_Cu_phys",
			Shield_Sn_log,false,0,fCheckOverlaps);
	G4VisAttributes* Shield_Cu_logVisAtt
		= new G4VisAttributes(G4Colour(1.0,0.50,0.60));
	Shield_Cu_log->SetVisAttributes(Shield_Cu_logVisAtt);


	//Shield Air
	G4double Shield_Air_Tubs_rMax = Shield_Cu_Tubs_rMax - Shield_Cu_Thickness;
	G4double Shield_Air_Tubs_dz = Shield_Cu_Tubs_dz - Shield_Cu_Thickness;

	G4VSolid * Shield_Air_tubs
		= new G4Tubs("Shield_Air_tubs",Shield_Tubs_rmin,Shield_Air_Tubs_rMax,Shield_Air_Tubs_dz,
				Shield_Tubs_sphi,Shield_Tubs_dphi);
	G4LogicalVolume * Shield_Air_log
		= new G4LogicalVolume(Shield_Air_tubs,Shield_Air,"Shield_Air",0,0,0);
	G4VPhysicalVolume *PVShieldAir =
		new G4PVPlacement(0,G4ThreeVector(),Shield_Air_log,"Shield_Air_phys",
				Shield_Cu_log,false,0,fCheckOverlaps);
	G4VisAttributes* Shield_Air_logVisAtt
		= new G4VisAttributes(G4Colour(1.0,0.3,1.0));
	Shield_Air_log->SetVisAttributes(Shield_Air_logVisAtt);

	G4bool shieldInvisible =1;

	if (shieldInvisible)
	{
		//Shield_Fe_log->SetVisAttributes(G4VisAttributes::Invisible);
		Shield_Pb_log->SetVisAttributes(G4VisAttributes::Invisible);
		Shield_Sn_log->SetVisAttributes(G4VisAttributes::Invisible);
		Shield_Cu_log->SetVisAttributes(G4VisAttributes::Invisible);
		//Shield_Air_log->SetVisAttributes(G4VisAttributes::Invisible);
	}

	return PVShieldAir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void   DetectorConstruction::ConstructHPGeDetector(G4LogicalVolume* matherLogicalVolume)
{
	G4NistManager* nist = G4NistManager::Instance();

	G4double Tubs_rmin = 0.*mm;
	G4double  Tubs_sphi =   0.*deg;
	G4double  Tubs_dphi = 360.*deg;

	G4double detectorRadius = 0.5 * 76. *mm;
	G4double detectorLength= 120. *mm;
	G4double GeCrystalRadius = 0.5 * 63. *mm;
	G4double GeCrystalLength = 44.1 *mm;
	//shield Al
	G4double detector_move_len = -140. *mm;
	G4double Shell_Al_Tubs_rMax =  detectorRadius;
	G4double Shell_Al_Tubs_dz = 0.5 * detectorLength;
	G4Material* Shell_Al = nist->FindOrBuildMaterial("G4_Al");

	G4VSolid * Shell_Al_tubs
		= new G4Tubs("Shell_Al_tubs",Tubs_rmin,Shell_Al_Tubs_rMax,Shell_Al_Tubs_dz,
				Tubs_sphi,Tubs_dphi);
	G4LogicalVolume * Shell_Al_log
		= new G4LogicalVolume(Shell_Al_tubs,Shell_Al,"Shell_Al_log",0,0,0);
	//G4VPhysicalVolume * physiWorldHPGe =
	new G4PVPlacement(0,G4ThreeVector(0.,0.,detector_move_len),Shell_Al_log,"Detector Shell Al phys",
			matherLogicalVolume,false,0,fCheckOverlaps);
	G4VisAttributes* Shell_Al_logVisAtt
		= new G4VisAttributes(G4Colour(0.40,0.0,0.70));
	Shell_Al_log->SetVisAttributes(Shell_Al_logVisAtt);

	//shell Galactic
	G4double Shell_Galactic_Tubs_rMax = Shell_Al_Tubs_rMax - shellAlThickness;
	G4double Shell_Galactic_Tubs_dz = Shell_Al_Tubs_dz - shellAlThickness;
	G4double GalacticTopThickness = 4. *mm;
	G4Material* Shell_Galactic = nist->FindOrBuildMaterial("G4_Galactic");

	G4VSolid * Shell_Galactic_tubs
		= new G4Tubs("Galactic_tubs",Tubs_rmin,Shell_Galactic_Tubs_rMax,Shell_Galactic_Tubs_dz,
				Tubs_sphi,Tubs_dphi);
	G4LogicalVolume * Shell_Galactic_log
		= new G4LogicalVolume(Shell_Galactic_tubs,Shell_Galactic,"Shell_Galactic_log",0,0,0);
	// G4VPhysicalVolume * tracker_phys =
	new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),Shell_Galactic_log,"Shell Galactic phys",
			Shell_Al_log,false,0,fCheckOverlaps);
	//  G4VisAttributes* Shell_Galactic_logVisAtt
	//    = new G4VisAttributes(G4Colour(0.3,0.7,1.0));
	//  Shell_Galactic_log->SetVisAttributes(Shell_Galactic_logVisAtt);

	//CUP Al
	G4double CUPThickness =0.8 *mm; 
	G4double CUPRadius = GeCrystalRadius + CUPThickness;
	G4double CUPLength =105. *mm;
	G4double CUPLength_dz =0.5 *CUPLength;
	G4double CUP_move_len = Shell_Galactic_Tubs_dz - (CUPLength_dz + GalacticTopThickness);
	G4VSolid * CUP_Al
		= new G4Tubs("CUP_Al",Tubs_rmin,CUPRadius,CUPLength_dz,
				Tubs_sphi,Tubs_dphi);
	G4LogicalVolume * CUP_Al_log
		= new G4LogicalVolume(CUP_Al,Shell_Al,"CUP_Al_log",0,0,0);
	//G4VPhysicalVolume * physiWorldHPGe =
	new G4PVPlacement(0,G4ThreeVector(0.,0.,CUP_move_len),CUP_Al_log,"Detector CUP Al phys",
			Shell_Galactic_log,false,0,fCheckOverlaps);
	G4VisAttributes* CUP_Al_logVisAtt
		= new G4VisAttributes(G4Colour(0.20,0.0,0.70));
	CUP_Al_logVisAtt ->G4VisAttributes::SetForceSolid(true);
	CUP_Al_log->SetVisAttributes(CUP_Al_logVisAtt);

	//CUP Galactic
	G4double CUP_Galactic_Radius = CUPRadius - CUPThickness ;
	G4double CUP_Galactic_Length = CUPLength  - 3.*mm -0.03*mm;
	G4double CUP_Galactic_Length_dz =0.5 * CUP_Galactic_Length;
	G4double CUP_Galactic_move_len = CUPLength_dz - (CUP_Galactic_Length_dz + 0.03*mm);
	G4VSolid * CUP_Galactic
		= new G4Tubs("CUP_Galactic",Tubs_rmin,CUP_Galactic_Radius ,CUP_Galactic_Length_dz,
				Tubs_sphi,Tubs_dphi);
	G4LogicalVolume * CUP_Galactic_log
		= new G4LogicalVolume(CUP_Galactic,Shell_Galactic,"CUP_Al_log",0,0,0);
	//G4VPhysicalVolume * physiWorldHPGe =
	new G4PVPlacement(0,G4ThreeVector(0.,0.,CUP_Galactic_move_len),CUP_Galactic_log,"Detector CUP_Galactic phys",
			CUP_Al_log/*logicWorld*/,false,0,fCheckOverlaps);
	G4VisAttributes* CUP_Galactic_logVisAtt
		= new G4VisAttributes(G4Colour(0.40,0.0,0.70));
	CUP_Galactic_log->SetVisAttributes(CUP_Galactic_logVisAtt);

	//HPGe dead_layer outer
	G4double HPGe_dead_layer_rMax = GeCrystalRadius;
	G4double HPGe_dead_layer_outer_dz = 0.5 * GeCrystalLength;
	G4double HPGe_move = CUP_Galactic_Length_dz - HPGe_dead_layer_outer_dz;

	G4VSolid * HPGe_dead_layer_outer
		= new G4Tubs("HPGe_dead_layer_outer",Tubs_rmin,HPGe_dead_layer_rMax,HPGe_dead_layer_outer_dz,
				Tubs_sphi,Tubs_dphi);
	G4LogicalVolume * HPGe_dead_layer_outer_log
		= new G4LogicalVolume(HPGe_dead_layer_outer,GeCrystal,"HPGe_dead_layer_log",0,0,0);
	// G4VPhysicalVolume * tracker_phys =
	new G4PVPlacement(0,G4ThreeVector(0. ,0. ,HPGe_move),HPGe_dead_layer_outer_log,"HPGe_dead_layer_phys",
			CUP_Galactic_log,false,0,fCheckOverlaps);
	G4VisAttributes* HPGe_dead_layer_logVisAtt
		= new G4VisAttributes(G4Colour(8.0,6.0,1.20));
	HPGe_dead_layer_logVisAtt->G4VisAttributes::SetForceSolid(true);
	HPGe_dead_layer_outer_log->SetVisAttributes(HPGe_dead_layer_logVisAtt);

	//Active HPGe Crystal
	G4double HPGe_detector_rMax = HPGe_dead_layer_rMax - outDeadLayerThickness;
	G4double HPGe_detector_Tubs_dz = HPGe_dead_layer_outer_dz - outDeadLayerThickness;
	//HPGe_detector_Ge = nist->FindOrBuildMaterial("G4_Ge");
	//G4double HPGe_move = -30.45 *mm;

	G4VSolid *ActiveHPGeCrystal_tubs
		= new G4Tubs("ActiveHPGeCrystal_tubs",Tubs_rmin,HPGe_detector_rMax,HPGe_detector_Tubs_dz,
				Tubs_sphi,Tubs_dphi);
	G4LogicalVolume *ActiveHPGeCrystal_log
		= new G4LogicalVolume(ActiveHPGeCrystal_tubs,GeCrystal,"HPGeDetector",0,0,0);
	// G4VPhysicalVolume * tracker_phys =
	new G4PVPlacement(0,G4ThreeVector(0. ,0. ,0.),ActiveHPGeCrystal_log,"HPGe_detector_phys",
			HPGe_dead_layer_outer_log,false,0,fCheckOverlaps);
	G4VisAttributes* HPGe_detector_logVisAtt
		= new G4VisAttributes(G4Colour(8.0,5.0,1.20));
	HPGe_detector_logVisAtt->G4VisAttributes::SetForceSolid(true);
	delete  ActiveHPGeCrystal_log->GetVisAttributes();
	ActiveHPGeCrystal_log->SetVisAttributes(HPGe_detector_logVisAtt);

	//HPGe detector inner dead_layer
	G4double GeDeadLayerInnerThickness = 0.3*um;
	G4double GeDeadLayerInnerRadius= 0.5* 10.8 *mm;
	G4double GeDeadLayerInnerDepth= 27.80 *mm;
	G4double HPGe_inner_dead_layer_rMax = GeDeadLayerInnerRadius+ GeDeadLayerInnerThickness;
	G4double HPGe_inner_dead_layer_Tubs_dz = 0.5* ( GeDeadLayerInnerDepth + GeDeadLayerInnerThickness );
	G4double HPGe_inner_dead_layer_move = -HPGe_detector_Tubs_dz + HPGe_inner_dead_layer_Tubs_dz;

	G4VSolid * HPGe_inner_dead_layer_tubs
		= new G4Tubs("HPGe_inner_dead_layer_tubs",Tubs_rmin,HPGe_inner_dead_layer_rMax,HPGe_inner_dead_layer_Tubs_dz,
				Tubs_sphi,Tubs_dphi);
	G4LogicalVolume * HPGe_inner_dead_layer_log
		= new G4LogicalVolume(HPGe_inner_dead_layer_tubs,GeCrystal,"HPGe_inner_dead_layer_log",0,0,0);
	// G4VPhysicalVolume * tracker_phys =
	new G4PVPlacement(0,G4ThreeVector(0. ,0. ,HPGe_inner_dead_layer_move),HPGe_inner_dead_layer_log,"HPGe_inner_dead_layer_phys",
			ActiveHPGeCrystal_log,false,0,fCheckOverlaps);
	G4VisAttributes* HPGe_inner_dead_layer_logVisAtt
		= new G4VisAttributes(G4Colour(3.0,6.0,1.20));
	HPGe_inner_dead_layer_log->SetVisAttributes(HPGe_inner_dead_layer_logVisAtt);

	//HPGe detector inner hole 
	G4double HPGe_inner_hole_rMax = GeDeadLayerInnerRadius;
	G4double HPGe_inner_hole_dz = 0.5 * GeDeadLayerInnerDepth;
	G4double HPGe_inner_hole_move = -HPGe_inner_dead_layer_Tubs_dz + HPGe_inner_hole_dz ;
	//G4Material* Mat_Cu = nist->FindOrBuildMaterial("G4_Cu");

	G4VSolid * HPGe_inner_hole
		= new G4Tubs("HPGe_inner_hole",Tubs_rmin,HPGe_inner_hole_rMax,HPGe_inner_hole_dz,
				Tubs_sphi,Tubs_dphi);
	G4LogicalVolume * HPGe_inner_hole_log
		= new G4LogicalVolume(HPGe_inner_hole,Shell_Galactic,"HPGe_inner_hole_log",0,0,0);
	// G4VPhysicalVolume * tracker_phys =
	new G4PVPlacement(0,G4ThreeVector(0. ,0. ,HPGe_inner_hole_move),HPGe_inner_hole_log,"HPGe_inner_hole_phys",
			HPGe_inner_dead_layer_log,false,0,fCheckOverlaps);
	G4VisAttributes* HPGe_inner_hole_logVisAtt
		= new G4VisAttributes(G4Colour(2.0,5.0,1.00));
	//HPGe_inner_hole_logVisAtt->G4VisAttributes::SetForceSolid(true);
	delete  HPGe_inner_hole_log->GetVisAttributes();
	HPGe_inner_hole_log->SetVisAttributes(HPGe_inner_hole_logVisAtt);

	G4bool detectorInvisible =1;

	if (detectorInvisible)
	{
		//Shell_Al_log -> SetVisAttributes(G4VisAttributes::Invisible);
		Shell_Galactic_log -> SetVisAttributes(G4VisAttributes::Invisible);
		CUP_Al_log -> SetVisAttributes(G4VisAttributes::Invisible);
		CUP_Galactic_log -> SetVisAttributes(G4VisAttributes::Invisible);
		HPGe_dead_layer_outer_log -> SetVisAttributes(G4VisAttributes::Invisible);
		//ActiveHPGeCrystal_log->SetVisAttributes(G4VisAttributes::Invisible);
		HPGe_inner_dead_layer_log -> SetVisAttributes(G4VisAttributes::Invisible);
		//HPGe_inner_tubs_Cu_log -> SetVisAttributes(G4VisAttributes::Invisible);
	}
	//--------- example of User Limits -------------------------------

	// below is an example of how to set tracking constraints in a given
	// logical volume(see also in N02PhysicsList how to setup the processes
	// G4StepLimiter or G4UserSpecialCuts).

	// Sets a max Step length in the tracker region, with G4StepLimiter
	//
	//G4double maxStep = HPGe_detector_Tubs_dz;
	//stepLimit = new G4UserLimits(maxStep);
	//HPGe_dead_layer_outer_log->SetUserLimits(stepLimit);

	// Set additional contraints on the track, with G4UserSpecialCuts
	//
	// G4double maxLength = 2*fTrackerLength, maxTime = 0.1*ns, minEkin = 10*MeV;
	// logicTracker->SetUserLimits(new G4UserLimits(maxStep,maxLength,maxTime,
	//                                               minEkin));
	//	return physiWorldHPGe;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void  DetectorConstruction::DefineMaterials()
{
	// Lead material defined using NIST Manager
	G4NistManager* nistManager = G4NistManager::Instance();
	G4bool fromIsotopes = false;

	// material
	Shield_Fe = nistManager->FindOrBuildMaterial("G4_Fe", fromIsotopes);
	Shield_Cu = nistManager->FindOrBuildMaterial("G4_Cu", fromIsotopes);
	Shield_Pb = nistManager->FindOrBuildMaterial("G4_Pb", fromIsotopes);
	Shield_Sn = nistManager->FindOrBuildMaterial("G4_Sn", fromIsotopes);
	GeCrystal = nistManager->FindOrBuildMaterial("G4_Ge", fromIsotopes);
	Shield_Air = nistManager->FindOrBuildMaterial("G4_AIR");

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
	SetSensitiveDetector("HPGeDetector",HPGeDetector);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetOutDeadLayerThickness(G4double value)
{
	outDeadLayerThickness = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetPbShield(G4bool value)
{
flagPbShield = value;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
