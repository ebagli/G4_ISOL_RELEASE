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
// $Id: DetectorConstruction.cc,v 1.15 2009-01-22 17:41:43 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-04-patch-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"

#include "globals.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4CutTubs.hh"
#include "G4Torus.hh"
#include "G4ThreeVector.hh"

#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4SubtractionSolid.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4LogicalVolume.hh"

#include "G4AssemblyVolume.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"

#include "G4SolidStore.hh"
#include "G4TransportationManager.hh"
#include "G4UnitsTable.hh"
#include "G4NistManager.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

#include "G4SDManager.hh"

#include "SensitiveDetector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
	//This function illustrates the possible ways to define materials
    
	G4String symbol;             //a=mass of a mole;
	G4double a, z, density;      //z=mean number of protons;
    
	G4int ncomponents, natoms;
	G4double fractionmass;
	G4double temperature, pressure;
    
	//
	// define Elements
	//
    
	G4Element* C  = new G4Element("Carbon",  symbol="C",  z= 6, a=  12.01*g/mole);
	G4Element* U  = new G4Element("Uranium",  symbol="U",  z= 92, a=  238*g/mole);

    	//Graphite

	G4Material* Graphite =
	new G4Material("Graphite", density= 1.9*g/cm3, ncomponents=1);
	Graphite->AddElement(C, fractionmass=1.);
 
  
	//Uranium carbide

	G4Material* UC4 = new G4Material("UC4", density= 4.*g/cm3, ncomponents=2);
	UC4->AddElement(C, natoms=4.);
	UC4->AddElement(U, natoms=1.);


	//Tantalum

	G4Material* Ta= new G4Material("Tantalum", z=73., a=180.95*g/mole, density=16.69*g/cm3);
    
    
    
	//
	//Vacuum
	//
  	// Get nist material manager
  	G4NistManager* nist = G4NistManager::Instance();
	G4Material* Air = nist->FindOrBuildMaterial("G4_AIR");  
    	density     = 1.e-5*g/cm3;
    	pressure    = 1.e-2*bar;
    	temperature = 2273.15*kelvin;         //from PhysicalConstants.h
    	G4Material* TVac = new G4Material("TVac", density, ncomponents=1,
                                      kStateGas,temperature,pressure);
   	 TVac->AddMaterial(Air, fractionmass=1.);
    
    
    
    
    	WorldMaterial = TVac;
    	TubeMaterial = Ta;
    	BoxMaterial = Graphite;
    	DiskMaterial = UC4;
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    DefineMaterials();
    
    
    // World
    G4double WorldSizeX = 50*cm; G4double WorldSizeYZ=50*cm;
    G4Box* solidWorld = new G4Box("World",
                                  WorldSizeX/2,
                                  WorldSizeYZ/2,
                                  WorldSizeYZ/2);
    
    G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld,		//its solid
                                                      WorldMaterial,	//its material
                                                      "World");		//its name
    
    G4PVPlacement* physiWorld = new G4PVPlacement(0,			//no rotation
                                                  G4ThreeVector(),	//at (0,0,0)
                                                  "World",		//its name
                                                  logicWorld,		//its logical volume
                                                  0,			//its mother  volume
                                                  false,			//no boolean operation
                                                  0);			//copy number
    
    // Disks
    G4double DiskRmin = 0.*mm;    G4double DiskRmax = 40*mm;    G4double DiskDz   = 0.8*mm;
    G4Tubs* sDisk = new G4Tubs("Disk1",				//name
                               DiskRmin/2,
                               DiskRmax/2,
                               DiskDz/2,
                               0.*deg,
                               360.*deg);	//dimensions
    
    G4LogicalVolume* lDisk = new G4LogicalVolume(sDisk,			//shape
                                                 DiskMaterial,		//material
                                                 "Disk1");		//name
    G4double ZposDisk[7] = {55.68*mm, 36.38*mm, 10.08*mm, -15.22*mm, -32.52*mm, -49.82*mm, -66.12*mm};
    
    for(G4int i0=0;i0<7;i0++){
        new G4PVPlacement(0,				//no rotation
                          G4ThreeVector(0.,0.,ZposDisk[i0]),		//at (0,0,0)
                          lDisk,			//logical volume
                          "Disk1",			//name
                          logicWorld,	       		//mother  volume
                          false,			//no boolean operation
                          i0);				//copy number
    }
    
    // Target Window
    G4double WindowRmin = 0*mm;  G4double WindowRmax = 45*mm;    G4double WindowDz= 0.2*mm;
    G4Tubs* sWindow = new G4Tubs("Window",				//name
                                 WindowRmin/2,
                                 WindowRmax/2,
                                 WindowDz/2,
                                 0.*deg,
                                 360.*deg);	//dimensions
    
    G4LogicalVolume* lWindow = new G4LogicalVolume(sWindow,		//shape
                                                   BoxMaterial,		//material
                                                   "Window");		//name
    
    new G4PVPlacement(0,				//no rotation
                      G4ThreeVector(0.,0.,-77.32),	//at (0,0,0)
                      lWindow,				//logical volume
                      "Window",				//
                      logicWorld,	       		//mother  volume
                      false,				//no boolean operation
                      0);				//copy number

    // Target Window spacer
    G4double Wsp1Rmin = 42*mm;  G4double Wsp1Rmax = 45*mm;    G4double Wsp1Dz= 29.4*mm;
    G4Tubs* sWsp1 = new G4Tubs("Wsp1",				//name
                                 Wsp1Rmin/2,
                                 Wsp1Rmax/2,
                                 Wsp1Dz/2,
                                 0.*deg,
                                 360.*deg);	//dimensions
    
    G4LogicalVolume* lWsp1 = new G4LogicalVolume(sWsp1,		//shape
                                                   BoxMaterial,		//material
                                                   "Wsp1");		//name
    
    new G4PVPlacement(0,				//no rotation
                      G4ThreeVector(0.,0.,-92.03),	//at (0,0,0)
                      lWsp1,				//logical volume
                      "Window spacer 1",				//
                      logicWorld,	       		//mother  volume
                      false,				//no boolean operation
                      0);				//copy number

    G4double Wsp2Rmin = 41*mm;  G4double Wsp2Rmax = 42*mm;    G4double Wsp2Dz= 12*mm;
    G4Tubs* sWsp2 = new G4Tubs("Wsp2",				//name
                                 Wsp2Rmin/2,
                                 Wsp2Rmax/2,
                                 Wsp2Dz/2,
                                 0.*deg,
                                 360.*deg);	//dimensions
    
    G4LogicalVolume* lWsp2 = new G4LogicalVolume(sWsp2,		//shape
                                                   BoxMaterial,		//material
                                                   "Wsp2");		//name
    
    new G4PVPlacement(0,				//no rotation
                      G4ThreeVector(0.,0.,-83.33),	//at (0,0,0)
                      lWsp2,				//logical volume
                      "Window spacer 2",				//
                      logicWorld,	       		//mother  volume
                      false,				//no boolean operation
                      0);				//copy number

    // Target dumps    
    
    G4double DumpRmin = 0*mm;  G4double DumpRmax = 45*mm;    G4double Dump12Dz= 0.8*mm; G4double Dump3Dz= 1*mm;
    
    G4Tubs* sDump12 = new G4Tubs("Dump12",				//name
                                 DumpRmin/2,
                                 DumpRmax/2,
                                 Dump12Dz/2,
                                 0.*deg,
                                 360.*deg);	//dimensions
    
    G4LogicalVolume* lDump12 = new G4LogicalVolume(sDump12,		//shape
                                                   BoxMaterial,		//material
                                                   "Dump12");		//name
    
    
    new G4PVPlacement(0,				//no rotation
                      G4ThreeVector(0.,0.,68.28),	//at (0,0,0)
                      lDump12,				//logical volume
                      "Dump1",				//name
                      logicWorld,	       		//mother  volume
                      false,				//no boolean operation
                      0);				//copy number
    
    new G4PVPlacement(0,				//no rotation
                      G4ThreeVector(0.,0.,75.08),	//at (0,0,0)
                      lDump12,				//logical volume
                      "Dump2",				//name
                      logicWorld,	       		//mother  volume
                      false,				//no boolean operation
                      0);				//copy number
    G4Tubs* sDump3 = new G4Tubs("Dump3",				//name
                                 DumpRmin/2,
                                 DumpRmax/2,
                                 Dump3Dz/2,
                                 0.*deg,
                                 360.*deg);	//dimensions
    
    G4LogicalVolume* lDump3 = new G4LogicalVolume(sDump3,		//shape
                                                   BoxMaterial,		//material
                                                   "Dump3");		//name
    
    
    new G4PVPlacement(0,				//no rotation
                      G4ThreeVector(0.,0.,83.88),	//at (0,0,0)
                      lDump3,				//logical volume
                      "Dump3",				//name
                      logicWorld,	       		//mother  volume
                      false,				//no boolean operation
                      0);				//copy number    

    // Target Dump spacers

    G4double Dsp1Rmin = 42*mm;  G4double Dsp1Rmax = 45*mm;    G4double Dsp1Dz= 6*mm;
    G4Tubs* sDsp1 = new G4Tubs("Dsp1",				//name
                                 Dsp1Rmin/2,
                                 Dsp1Rmax/2,
                                 Dsp1Dz/2,
                                 0.*deg,
                                 360.*deg);	//dimensions
    
    G4LogicalVolume* lDsp1 = new G4LogicalVolume(sDsp1,		//shape
                                                   BoxMaterial,		//material
                                                   "Dsp1");		//name
    
    new G4PVPlacement(0,				//no rotation
                      G4ThreeVector(0.,0.,71.68),	//at (0,0,0)
                      lDsp1,				//logical volume
                      "Dump spacer 1",			//
                      logicWorld,	       		//mother  volume
                      false,				//no boolean operation
                      0);				//copy number

    G4double Dsp2Rmin = 42*mm;  G4double Dsp2Rmax = 45*mm;    G4double Dsp2Dz= 8*mm;
    G4Tubs* sDsp2 = new G4Tubs("Dsp2",				//name
                                 Dsp2Rmin/2,
                                 Dsp2Rmax/2,
                                 Dsp2Dz/2,
                                 0.*deg,
                                 360.*deg);	//dimensions
    
    G4LogicalVolume* lDsp2 = new G4LogicalVolume(sDsp2,		//shape
                                                   BoxMaterial,		//material
                                                   "Dsp2");		//name
    
    new G4PVPlacement(0,				//no rotation
                      G4ThreeVector(0.,0.,79.48),	//at (0,0,0)
                      lDsp2,				//logical volume
                      "Dump spacer 2",			//
                      logicWorld,	       		//mother  volume
                      false,				//no boolean operation
                      0);				//copy number

    G4double Dsp3Rmin = 34*mm;  G4double Dsp3Rmax = 45*mm;    G4double Dsp3Dz= 10.1*mm;
    G4Tubs* sDsp3 = new G4Tubs("Dsp3",				//name
                                 Dsp3Rmin/2,
                                 Dsp3Rmax/2,
                                 Dsp3Dz/2,
                                 0.*deg,
                                 360.*deg);	//dimensions
    
    G4LogicalVolume* lDsp3 = new G4LogicalVolume(sDsp3,		//shape
                                                   BoxMaterial,		//material
                                                   "Dsp2");		//name
    
    new G4PVPlacement(0,				//no rotation
                      G4ThreeVector(0.,0.,89.43),	//at (0,0,0)
                      lDsp3,				//logical volume
                      "Dump spacer 3",			//
                      logicWorld,	       		//mother  volume
                      false,				//no boolean operation
                      0);				//copy number

    // TRANSFER LINE: SECTION 1 + SECTION 2

    // Section 1

    G4double Sqrt3 = std::sqrt(3);
    G4double TransferRmin = 8*mm; G4double TransferRmax= 8.8*mm;    G4double TranSect1Dz= 60/Sqrt3*mm;


    G4CutTubs* Section1 = new G4CutTubs("Section1",				//name
                                       TransferRmin/2,
                                       TransferRmax/2,
                                       TranSect1Dz/2,
                                       0.*deg,
                                       360.*deg,
                                       G4ThreeVector(0,0,-1),
                                       G4ThreeVector(0,-0.5,0.5*Sqrt3));


    G4double CutSecRmin  = 0*mm;    G4double CutSecRmax  =50*mm;    G4double CutSecDz    =10*mm;
    G4Tubs* CutSection1 = new G4Tubs("CutSection1",				//name
                               CutSecRmin/2,
                               CutSecRmax/2,
                               CutSecDz/2,
                               0.*deg,
                               360.*deg);	//dimensions
 
    G4RotationMatrix* yRot = new G4RotationMatrix();
    yRot->rotateY(90.*deg);

    G4VSolid* sTransferline1 = new G4SubtractionSolid ("Transferline1",
    							Section1,
    							CutSection1,
    							yRot,
    							G4ThreeVector(0,0,-TranSect1Dz/2));

    G4LogicalVolume* lTransferline1 = new G4LogicalVolume(sTransferline1,		//shape
                                                         TubeMaterial,			//material
                                                         "Transferline1");		//name
    G4RotationMatrix* tRot = new G4RotationMatrix();
    tRot->rotateY(-90.*deg); 
    tRot->rotateX(60.*deg); 
   
    
    new G4PVPlacement(tRot,				// rotation
                      G4ThreeVector(15/Sqrt3*mm,15*mm,0*mm), 			//at (0,0,0)
                      lTransferline1,			//logical volume
                      "Transferline1",			//name
                      logicWorld,	       		//mother  volume
                      false,				//no boolean operation
                      0);				//copy number


    // Section 2

    G4double TranSect2Dz= 90-30/Sqrt3*mm;

    G4CutTubs* sTransferline2 = new G4CutTubs("sTransferline2",				//name
                                       TransferRmin/2,
                                       TransferRmax/2,
                                       TranSect2Dz/2,
                                       0.*deg,
                                       360.*deg,
                                       G4ThreeVector(0,-0.5/Sqrt3,-0.5),
                                       G4ThreeVector(0,0,1)); 
 
    G4LogicalVolume* lTransferline2 = new G4LogicalVolume(sTransferline2,		//shape
                                                         TubeMaterial,			//material
                                                         "Transferline2");		//name  
    G4RotationMatrix* tRot2 = new G4RotationMatrix();
    tRot2->rotateY(-90.*deg); 


    new G4PVPlacement(tRot2,								// rotation
                      G4ThreeVector(45+15/Sqrt3*mm,30*mm,0*mm), 			//at (0,0,0)
                      lTransferline2,							//logical volume
                      "Transferline2",							//name
                      logicWorld,	       						//mother  volume
                      false,								//no boolean operation
                      0);								//copy number

    
    // TARGET HEATER
    
    
    // Tantalum heater hole

    G4double ThHoleR = 8*mm; G4double ThHoleDz= 30*mm; 

    G4Tubs* sThHole = new G4Tubs("ThHole",				//name
                                 0,
				 ThHoleR/2,
                                 ThHoleDz,
                                 0.*deg,
                                 360.*deg);	                	//dimensions
        
       
    // Tantalum Heater
    G4double TubeRmin  = 49.6*mm;    G4double TubeRmax  =50*mm;    G4double TubeDz    =169*mm;
    G4Tubs* sTaTube = new G4Tubs("TaTube",				//name
                               TubeRmin/2,
                               TubeRmax/2,
                               TubeDz/2,
                               0.*deg,
                               360.*deg);	//dimensions
    
    
    
    G4SubtractionSolid* sTaHeater = new G4SubtractionSolid("TaHeater",
                                                               sTaTube,
                                                               sThHole,
                                                               tRot,
                                                               G4ThreeVector(15/Sqrt3*mm,15*mm,21.5*mm));
    
    G4LogicalVolume* lTaHeater = new G4LogicalVolume(sTaHeater,			//shape
                                                         TubeMaterial,		//material
                                                         "TaHeater");		//name
    
    new G4PVPlacement(0,				//no rotation
                      G4ThreeVector(0,0,-21.5),		//at (0,0,0)
                      lTaHeater,			//logical volume
                      "Tantalum Heater",		//name
                      logicWorld,	       		//mother  volume
                      false,				//no boolean operation
                      0);				//copy number
    
    
    
    
    // GRAPHITE BOX
    
    // Graphite box hole

    G4double GbHoleR = 10*mm; G4double GbHoleDz= 30*mm; 

    G4Tubs* sGbHole = new G4Tubs("GbHole",				//name
                                 0,
				 GbHoleR/2,
                                 GbHoleDz,
                                 0.*deg,
                                 360.*deg);	                	//dimensions
            
    // Graphite Box

    G4double BoxRmin  = 45*mm;    G4double BoxRmax  =49.5*mm;    G4double BoxDz    =201.5*mm;
    G4Tubs* sBox = new G4Tubs("Box",				//name
                              BoxRmin/2,
                              BoxRmax/2,
                              BoxDz/2,
                              0.*deg,
                              360.*deg);	//dimensions
    
    
    G4SubtractionSolid* CBox = new G4SubtractionSolid("CBox",
                                                               sBox,
                                                               sGbHole,
                                                               tRot,
                                                               G4ThreeVector(15/Sqrt3*mm,15*mm,6.27*mm));
    
    G4LogicalVolume* lCBox = new G4LogicalVolume(CBox,			//shape
                                                         BoxMaterial,		//material
                                                         "CBox");		//name
    
    new G4PVPlacement(0,				//no rotation
                      G4ThreeVector(0,0,-6.27),			//at (0,0,0)
                      lCBox,			//logical volume
                      "CBox",			//name
                      logicWorld,	       		//mother  volume
                      false,				//no boolean operation
                      0);				//copy number
    
    
   
    
    // Detector
    
    G4double detX = 1.*mm;
    G4double detYZ = 30.*mm;
    G4Box* solidDetector = new G4Box("detector",
                                  detX/2,
                                  detYZ/2,
                                  detYZ/2);
    
    G4LogicalVolume* logicDetector = new G4LogicalVolume(solidDetector,		//its solid
                                                      WorldMaterial,	//its material
                                                      "detector");		//its name
    
    new G4PVPlacement(0,			//no rotation
                                                  G4ThreeVector(90.5*mm,30.*mm,0.),	//at (0,0,0)
                      logicDetector,		//its logical volume
                                                  "detector",		//its name
                                                  logicWorld,			//its mother  volume
                                                  false,			//no boolean operation
                                                  0);			//copy number

    


//
// VISUALIZATION ATTRIBUTES
//  
 


  	G4Colour white   (1.0, 1.0, 1.0);
 	G4Colour gray    (0.5, 0.5, 0.5);
  	G4Colour black   (0.0, 0.0, 0.0);
  	G4Colour red     (1.0, 0.0, 0.0);
 	G4Colour green   (0.0, 0.1, 0.0);
  	G4Colour blue    (0.0, 0.0, 1.0);
  	G4Colour cyan    (0.0, 0.1, 1.0);
  	G4Colour magenta (1.0, 0.0, 1.0);
 	G4Colour yellow  (1.0, 1.0, 0.0);
  	G4Colour orange  (1.0, 0.5, 0.0);
   


  	G4VisAttributes* worldVisAtt = new G4VisAttributes(white);
  	worldVisAtt->SetVisibility(true);
  	logicWorld->SetVisAttributes(worldVisAtt); 

 	G4VisAttributes* DiskVisAtt = new G4VisAttributes(yellow);
 	DiskVisAtt->SetVisibility(true);
 	lDisk->SetVisAttributes(DiskVisAtt);

 	G4VisAttributes* BoxVisAtt = new G4VisAttributes(black);
 	BoxVisAtt->SetVisibility(true);
 	lCBox->SetVisAttributes(BoxVisAtt);

	G4VisAttributes* Dump12VisAtt = new G4VisAttributes(black);
  	Dump12VisAtt->SetVisibility(true);
  	lDump12->SetVisAttributes(Dump12VisAtt);

  	G4VisAttributes* Dump3VisAtt = new G4VisAttributes(black);
  	Dump3VisAtt->SetVisibility(true);
  	lDump3->SetVisAttributes(Dump3VisAtt);

  	G4VisAttributes* Dsp1VisAtt = new G4VisAttributes(black);
  	Dsp1VisAtt->SetVisibility(true);
  	lDsp1->SetVisAttributes(Dsp1VisAtt);

  	G4VisAttributes* Dsp2VisAtt = new G4VisAttributes(black);
  	Dsp2VisAtt->SetVisibility(true);
  	lDsp2->SetVisAttributes(Dsp2VisAtt);

  	G4VisAttributes* Dsp3VisAtt = new G4VisAttributes(black);
  	Dsp3VisAtt->SetVisibility(true);
  	lDsp3->SetVisAttributes(Dsp3VisAtt);

  	G4VisAttributes* Wsp1VisAtt = new G4VisAttributes(black);
  	Wsp1VisAtt->SetVisibility(true);
  	lWsp1->SetVisAttributes(Wsp1VisAtt);

  	G4VisAttributes* Wsp2VisAtt = new G4VisAttributes(black);
  	Wsp2VisAtt->SetVisibility(true);
  	lWsp2->SetVisAttributes(Wsp2VisAtt);

  	G4VisAttributes* WindowVisAtt = new G4VisAttributes(black);
  	WindowVisAtt->SetVisibility(true);
  	lWindow->SetVisAttributes(WindowVisAtt);

  	G4VisAttributes* TaHeaterVisAtt = new G4VisAttributes(gray);
  	TaHeaterVisAtt->SetVisibility(true);
  	lTaHeater->SetVisAttributes(TaHeaterVisAtt);
 
  	G4VisAttributes* Transferline1VisAtt = new G4VisAttributes(gray);
  	Transferline1VisAtt->SetVisibility(true);
  	lTransferline1->SetVisAttributes(Transferline1VisAtt);
  
  	G4VisAttributes* Transferline2VisAtt = new G4VisAttributes(gray);
  	Transferline2VisAtt->SetVisibility(true);
  	lTransferline2->SetVisAttributes(Transferline2VisAtt);
    
    return physiWorld;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField(){
    G4String SDname;
    G4VSensitiveDetector* telescope = new SensitiveDetector(SDname="/telescope");
    G4SDManager::GetSDMpointer()->AddNewDetector(telescope);
    
    G4LogicalVolume* fDetLogic = G4LogicalVolumeStore::GetInstance()->GetVolume("detector");
    fDetLogic->SetSensitiveDetector(telescope);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

