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
#include "G4Orb.hh"
#include "G4Torus.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
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
	//Pb Shield
	Shield_Length = 630.*mm;
	Shield_rMax = 0.5*510.*mm;
	Shield_Fe_Thickness = 9.5*mm;
	Shield_Pb_Thickness = 101.*mm;
	Shield_Cu_Thickness= 1.6 *mm;
	Shield_Sn_Thickness= 0.5*mm;
	//cover
	coverThick = 3.0 *mm;
	//detector shell
	shellRadius = 0.5 * 76. *mm;
	shellLength= 120. *mm;
	shellThick =1. *mm;
	endGap =4.0 *mm;
	detectorMove = -0.5 * shellLength;
	//CUP
	CUPLength =105.*mm;
	CUPThick =0.8 *mm;
	CUPTopThick =0.03 *mm;
	CUPBottomThick =3. *mm;
	mylarThick =0.03 *mm;
	//Ge crystal
	crystalRadius = 0.5 * 63. *mm;
	crystalHalfLength =0.5* 44.1 *mm;
	crystalEndRadius = 8. *mm;
	holeDepth = 27.8 *mm;
	holeRadius =0.5* 10.8 *mm;
	outerDeadLayerThick = 0.7 *mm;
	innerDeadLayerThick = 0.3 *um;
	//sample	
	shapeRad = 0.564 *cm;
	shapeHalfDepth = 0.5 *cm;
	sampleMove = 90.*mm;

	DefineMaterials();

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
	G4LogicalVolume* pLV = 0;
	physiWorld = ConstructWorld();

	pLV = physiWorld->GetLogicalVolume();
	if(flagPbShield == true)
	{
		G4VPhysicalVolume* PVShieldDauther = ConstructPbShield(pLV);
		pLV = PVShieldDauther->GetLogicalVolume();
	}
	ConstructHPGeDetector(pLV);
	ConstructSample(pLV);
	return physiWorld;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume*  DetectorConstruction::ConstructWorld()
{
	//--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------
	//------------------------------
	// World
	//------------------------------
	G4double fWorldLength = 1.2 * Shield_Length;
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
	G4double  Shield_Tubs_sphi =   0.*deg;
	G4double  Shield_Tubs_dphi = 360.*deg;


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
//sample Geo
void  DetectorConstruction::ConstructSample(G4LogicalVolume* motherLogicalVolume)
{
	G4NistManager* nist = G4NistManager::Instance();
	sampleMaterial= nist->FindOrBuildMaterial("G4_Galactic");

	G4VSolid * sampleShape
		= new G4Tubs("sampleShape", 0.*cm,shapeRad ,shapeHalfDepth ,
				0. *deg,360. *deg);
	G4LogicalVolume * logSample 
		= new G4LogicalVolume(sampleShape,sampleMaterial,"sample_log",0,0,0);

	//	G4VPhysicalVolume * physiWorldPbShield =
	new G4PVPlacement(0,G4ThreeVector(0.,0.,sampleMove +shapeHalfDepth ),logSample,"sample_phys",
			motherLogicalVolume,false,0,fCheckOverlaps);

	G4VisAttributes* visAttSample
		= new G4VisAttributes(G4Colour(0.3,0.2,1.0));
	visAttSample->G4VisAttributes::SetForceSolid(true);
	logSample->SetVisAttributes(visAttSample);

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void   DetectorConstruction::ConstructHPGeDetector(G4LogicalVolume* matherLogicalVolume)
{
	G4NistManager* nist = G4NistManager::Instance();
	G4Material* shellAl = nist->FindOrBuildMaterial("G4_Al");
	G4Material* vacuum = nist->FindOrBuildMaterial("G4_Galactic");
	G4Material* coverMat= nist->FindOrBuildMaterial("G4_PLEXIGLASS");
	G4Material* mylar= nist->FindOrBuildMaterial("G4_MYLAR");

	G4double  sphi =   0.*deg;
	G4double  dphi = 360.*deg;

	//Cover
	G4VSolid *cover = new G4Tubs("cover",
			0. *mm,
			shellRadius ,
			0.5* coverThick ,	
			sphi,
			dphi);
	G4LogicalVolume * logCover
		= new G4LogicalVolume(cover, coverMat,"logCover",0,0,0);
	new G4PVPlacement(0,G4ThreeVector(0. ,0. ,detectorMove + 0.5 *shellLength+ 0.5* coverThick),logCover,"physiCover",
			matherLogicalVolume,false,0,fCheckOverlaps);

	//detector
	G4VSolid *HPGe = new G4Tubs("HPGe",
			0.*mm,
			shellRadius,
			0.5*shellLength,
			sphi,
			dphi);
	G4LogicalVolume * logHPGe
		= new G4LogicalVolume(HPGe,vacuum,"logHPGe",0,0,0);
	new G4PVPlacement(0,G4ThreeVector(0. ,0. ,detectorMove),logHPGe,"physiHPGe",
			matherLogicalVolume,false,0,fCheckOverlaps);

	//Crystal
	G4VSolid *crystal1 = new G4Tubs("cyl1",
			0. *mm,
			crystalRadius - crystalEndRadius ,
			crystalHalfLength ,	
			sphi,
			dphi);
	G4VSolid *crystal2 = new G4Tubs("cyl2",
			0. *mm,
			crystalRadius ,
			crystalHalfLength - 0.5 * crystalEndRadius ,	
			sphi,
			dphi);
	G4VSolid *tor1 = new G4Torus("tor1",
			0. *mm,
			crystalEndRadius,
			crystalRadius - crystalEndRadius,
			sphi,
			dphi);
	G4VSolid *crystal3 = new G4UnionSolid("cry3",
			crystal1,crystal2,0,G4ThreeVector(0., 0., -0.5*(crystalEndRadius)));
	G4VSolid *crystal4 = new G4UnionSolid("cry4",
			crystal3,tor1,0,G4ThreeVector(0., 0., crystalHalfLength - crystalEndRadius));

	//Making the Active Crystal shap
	G4double activeRadius = crystalRadius - outerDeadLayerThick;
	G4double activeHalfLength = crystalHalfLength - 0.5*outerDeadLayerThick;
	G4double activeEndRadius = crystalEndRadius -outerDeadLayerThick;
	G4VSolid *activeCrystal1 = new G4Tubs("acyl1",
			0. *mm,
			activeRadius - activeEndRadius ,
			activeHalfLength ,	
			sphi,
			dphi);
	G4VSolid *activeCrystal2 = new G4Tubs("acyl2",
			0. *mm,
			activeRadius ,
			activeHalfLength - 0.5 * activeEndRadius ,	
			sphi,
			dphi);
	G4VSolid *activeTor1 = new G4Torus("activeTor1",
			0. *mm,
			activeEndRadius,
			activeRadius - activeEndRadius,
			sphi,
			dphi);
	G4VSolid *activeCrystal3 = new G4UnionSolid("cry3",
			activeCrystal1,activeCrystal2,0,G4ThreeVector(0., 0., -0.5*(activeEndRadius)));
	G4VSolid *activeCrystal4 = new G4UnionSolid("cry4",
			activeCrystal3,activeTor1,0,G4ThreeVector(0., 0., activeHalfLength - activeEndRadius));

	//making outer dead layer
	G4VSolid *outerDeadLayer = new G4SubtractionSolid("outerDeadLayer",
			crystal4 ,activeCrystal4 ,0,G4ThreeVector(0., 0., -0.5*outerDeadLayerThick));

	//making the hole
	G4VSolid *hole1 = new G4Tubs("hole1",
			0.*mm,
			holeRadius,
			0.5*(holeDepth- holeRadius),
			sphi,
			dphi);
	G4VSolid *hole2 = new G4Orb("hole2",holeRadius);
	G4VSolid *hole = new G4UnionSolid("hole",
			hole1,hole2,0 ,  G4ThreeVector(0., 0.,0.5*(holeDepth- holeRadius) ) );

	//Making the inner dead layer
	G4VSolid *innerDead1 = new G4Tubs("innerDead1",
			0.*mm,
			holeRadius+ innerDeadLayerThick,
			0.5*(holeDepth- holeRadius),
			sphi,
			dphi);
	G4VSolid *innerDead2 = new G4Orb("innerDead2",holeRadius +innerDeadLayerThick);
	G4VSolid *innerDead3 = new G4UnionSolid("innerDead3",
			innerDead1,innerDead2,0 ,  G4ThreeVector(0., 0.,0.5*(holeDepth- holeRadius) ) );
	G4VSolid *innerDeadLayer = new G4SubtractionSolid("innerDeadLayer",
			innerDead3,hole,0,G4ThreeVector(0., 0.,0.));

	//Making final detector shape
	G4VSolid * activeCrystal = new G4SubtractionSolid("activeCrystal",
			activeCrystal4 ,innerDead3,0,G4ThreeVector(0., 0.,-activeHalfLength+ 0.5*(holeDepth - holeRadius) ) );

	G4LogicalVolume * logOuterDeadLayer
		= new G4LogicalVolume(outerDeadLayer,GeCrystal,"logOuterDeadLayer",0,0,0);
	G4LogicalVolume * logInnerDeadLayer
		= new G4LogicalVolume(innerDeadLayer,GeCrystal,"logInnerDeadLayer",0,0,0);
	G4LogicalVolume * logActiveCrystal
		= new G4LogicalVolume(activeCrystal,GeCrystal,"logActiveCrystal",0,0,0);

	//mylar
	G4VSolid *mylarLayer = new G4Tubs("mylarLayer",
			0.*mm,
			CUPThick+ crystalRadius ,
			0.5*mylarThick,
			sphi,
			dphi);
	G4LogicalVolume * logMylar
		= new G4LogicalVolume(mylarLayer,mylar ,"logMylar",0,0,0);
	//CUP
	G4VSolid *CUP1 = new G4Tubs("CUP1",
			0.*mm,
			CUPThick+ crystalRadius ,
			0.5*CUPLength,
			sphi,
			dphi);
	G4VSolid *CUP2 = new G4Tubs("CUP2",
			0.*mm,
			crystalRadius ,
			0.5*(CUPLength- CUPTopThick- CUPBottomThick ),
			sphi,
			dphi);
	G4VSolid * CUP = new G4SubtractionSolid("CUP",
			CUP1 ,CUP2 ,0,G4ThreeVector(0., 0.,0.5 *(CUPBottomThick -CUPTopThick)) );
	G4LogicalVolume * logCUP
		= new G4LogicalVolume(CUP,shellAl,"logCUP",0,0,0);
	//detector shell
	G4VSolid *shell1 = new G4Tubs("shell1",
			0.*mm,
			shellRadius,
			0.5*shellLength,
			sphi,
			dphi);
	G4VSolid *shell2 = new G4Tubs("shell2",
			0.*mm,
			shellRadius-shellThick,
			0.5*shellLength - shellThick,
			sphi,
			dphi);
	G4VSolid * shell = new G4SubtractionSolid("shell",
			shell1 ,shell2 ,0,G4ThreeVector(0., 0.,0.0) );
	G4LogicalVolume * logShell
		= new G4LogicalVolume(shell,shellAl,"logShell",0,0,0);

	new G4PVPlacement(0,G4ThreeVector(),logShell,"physiShell",
			logHPGe,false,0,fCheckOverlaps);

	new G4PVPlacement(0,
			G4ThreeVector(0., 0.,0.5*(shellLength + mylarThick)-shellThick-endGap),
			logMylar,"physiMylarLayer",
			logHPGe,false,0,fCheckOverlaps);
	new G4PVPlacement(0,
			G4ThreeVector(0., 0.,0.5*(shellLength - CUPLength)-shellThick-endGap),
			logCUP,"physiCUP",
			logHPGe,false,0,fCheckOverlaps);
	new G4PVPlacement(0,
			G4ThreeVector(0., 0.,0.5*shellLength - crystalHalfLength -shellThick-endGap -CUPTopThick ),
			logOuterDeadLayer,"physiOuterDeadLayer",
			logHPGe,false,0,fCheckOverlaps);
	new G4PVPlacement(0,
			G4ThreeVector(0., 0.,0.5*shellLength - activeHalfLength -shellThick-endGap -CUPTopThick -outerDeadLayerThick),
			logActiveCrystal,"physiActiveCrystal",
			logHPGe,false,0,fCheckOverlaps);
	new G4PVPlacement(0,
			G4ThreeVector(0., 0.,0.5*(shellLength - holeDepth+holeRadius)-shellThick-endGap -CUPTopThick -(2.*crystalHalfLength -holeDepth+holeRadius)),
			logInnerDeadLayer,"physiInnerDeadLayer",
			logHPGe,false,0,fCheckOverlaps);

	//Detector Visualization Attributes
	G4VisAttributes* HPGeVisAtt
		= new G4VisAttributes(G4Colour(1.0,1.0,0.00));
	logHPGe->SetVisAttributes(HPGeVisAtt);

	G4VisAttributes* shellVisAtt
		= new G4VisAttributes(G4Colour(1.0,1.0,0.00));
	logShell->SetVisAttributes(shellVisAtt);

	G4VisAttributes* CUPVisAtt
		= new G4VisAttributes(G4Colour(0.2,1.0,0.00));
	logCUP->SetVisAttributes(CUPVisAtt);

	G4VisAttributes* outerDeadLayerVisAtt
		= new G4VisAttributes(G4Colour(0.9,1.0,0.0,0.80));
	outerDeadLayerVisAtt->G4VisAttributes::SetForceSolid(true);
	logOuterDeadLayer->SetVisAttributes(outerDeadLayerVisAtt);

	G4VisAttributes* activeCrystalVisAtt
		= new G4VisAttributes(G4Colour(0.0,1.0,0.3,0.20));
	activeCrystalVisAtt->G4VisAttributes::SetForceSolid(true);
	logActiveCrystal->SetVisAttributes(activeCrystalVisAtt);

	G4VisAttributes* innerDeadLayerVisAtt
		= new G4VisAttributes(G4Colour(0.9,1.0,0.0,0.80));
	innerDeadLayerVisAtt->G4VisAttributes::SetForceSolid(true);
	logInnerDeadLayer->SetVisAttributes(innerDeadLayerVisAtt);

	G4bool detectorInvisible =0;
	if (detectorInvisible)
	{
		logOuterDeadLayer-> SetVisAttributes(G4VisAttributes::Invisible);
		logActiveCrystal->SetVisAttributes(G4VisAttributes::Invisible);
		logInnerDeadLayer-> SetVisAttributes(G4VisAttributes::Invisible);
	}
	G4bool shellInvisible= 1;
	if(shellInvisible)
	{
		logHPGe-> SetVisAttributes(G4VisAttributes::Invisible);
		logShell-> SetVisAttributes(G4VisAttributes::Invisible);
		logCUP-> SetVisAttributes(G4VisAttributes::Invisible);
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
	SetSensitiveDetector("logActiveCrystal",HPGeDetector);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetPbShield(G4bool value)
{
	flagPbShield = value;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetOutDeadLayerThickness(G4double value)
{
	outerDeadLayerThick= value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetCoverThick(G4double value)
{
	coverThick = value;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetShellRadius(G4double value)
{
	shellRadius= value;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetShellLength(G4double value)
{
	shellLength= value;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetShellThick(G4double value)
{
	shellThick = value;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetEndGap(G4double value)
{
	endGap = value;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetCUPLength(G4double value)
{
	CUPLength = value;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetCUPThick(G4double value)
{
	CUPThick = value;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetCUPTopThick(G4double value)
{
	CUPTopThick = value;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetCUPBottomThick(G4double value)
{
	CUPBottomThick = value;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetMylarThick (G4double value)
{
	mylarThick = value;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetCrystalRadius (G4double value)
{
	crystalRadius = value;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetCrystalHalfLength (G4double value)
{
	crystalHalfLength = value;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetCrystalEndRadius (G4double value)
{
	crystalEndRadius = value;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetHoleDepth(G4double value)
{
	holeDepth = value;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetHoleRadius(G4double value)
{
	holeRadius= value;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetOuterDeadLayerThick(G4double value)
{
	outerDeadLayerThick = value;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetInnerDeadLayerThick(G4double value)
{
	innerDeadLayerThick = value;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"
void DetectorConstruction::UpdateGeometry()
{
	G4RunManager::GetRunManager()->DefineWorldVolume(this->Construct());
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
