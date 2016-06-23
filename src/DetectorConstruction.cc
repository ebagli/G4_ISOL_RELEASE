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

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Torus.hh"

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
    
    G4Element* H  = new G4Element("Hydrogen",symbol="H",  z= 1, a=   1.01*g/mole);
    G4Element* C  = new G4Element("Carbon",  symbol="C",  z= 6, a=  12.01*g/mole);
    G4Element* N  = new G4Element("Nitrogen",symbol="N",  z= 7, a=  14.01*g/mole);
    G4Element* O  = new G4Element("Oxygen",  symbol="O",  z= 8, a=  16.00*g/mole);
    G4Element* Na = new G4Element("Sodium",  symbol="Na", z=11, a=  22.99*g/mole);
    G4Element* Ar = new G4Element("Argon",   symbol="Ar", z=18, a=  39.95*g/mole);
    G4Element* I  = new G4Element("Iodine",  symbol="I" , z=53, a= 126.90*g/mole);
    
    //
    // define simple materials
    //
    
    new G4Material("H2Liq"    , z= 1, a= 1.01*g/mole, density= 70.8*mg/cm3);
    new G4Material("Beryllium", z= 4, a= 9.01*g/mole, density= 1.848*g/cm3);
    new G4Material("Aluminium", z=13, a=26.98*g/mole, density= 2.700*g/cm3);
    new G4Material("Silicon"  , z=14, a=28.09*g/mole, density= 2.330*g/cm3);
    
    G4Material* lAr =
    new G4Material("liquidArgon", density= 1.390*g/cm3, ncomponents=1);
    lAr->AddElement(Ar, natoms=1);
    
    new G4Material("Iron",     z=26, a= 55.85*g/mole, density= 7.870*g/cm3);
    new G4Material("Copper",   z=29, a= 63.55*g/mole, density= 8.960*g/cm3);
    new G4Material("Germanium",z=32, a= 72.61*g/mole, density= 5.323*g/cm3);
    new G4Material("Silver",   z=47, a=107.87*g/mole, density= 10.50*g/cm3);
    new G4Material("Tungsten", z=74, a=183.85*g/mole, density= 19.30*g/cm3);
    new G4Material("Gold",     z=79, a=196.97*g/mole, density= 19.32*g/cm3);
    new G4Material("Lead",     z=82, a=207.19*g/mole, density= 11.35*g/cm3);
    
    // define a material from elements.   case 1: chemical molecule
    G4Material* H2O = new G4Material("Water", density= 1.000*g/cm3, ncomponents=2);
    H2O->AddElement(H, natoms=2);
    H2O->AddElement(O, natoms=1);
    H2O->GetIonisation()->SetMeanExcitationEnergy(75*eV);
    
    G4Material* CH = new G4Material("Plastic", density= 1.04*g/cm3, ncomponents=2);
    CH->AddElement(C, natoms=1);
    CH->AddElement(H, natoms=1);
    
    G4Material* NaI = new G4Material("NaI", density= 3.67*g/cm3, ncomponents=2);
    NaI->AddElement(Na, natoms=1);
    NaI->AddElement(I , natoms=1);
    NaI->GetIonisation()->SetMeanExcitationEnergy(452*eV);
    
    // define a material from elements.   case 2: mixture by fractional mass
    G4Material* Air = new G4Material("Air", density= 1.290*mg/cm3, ncomponents=2);
    Air->AddElement(N, fractionmass=0.7);
    Air->AddElement(O, fractionmass=0.3);
    
    G4Material* Air20 =
    new G4Material("Air20", density= 1.205*mg/cm3, ncomponents=2,
                   kStateGas, 293.*kelvin, 1.*atmosphere);
    Air20->AddElement(N, fractionmass=0.7);
    Air20->AddElement(O, fractionmass=0.3);
    
    //graphite
    G4Material* Graphite =
    new G4Material("Graphite", density= 1.9*g/cm3, ncomponents=1);
    Graphite->AddElement(C, fractionmass=1.);
    
    //tantalum
    G4Material* Ta=
    new G4Material("Tantalum", z=73., a=180.95*g/mole, density=16.69*g/cm3);
    
    //target
    G4Material* Target=
    new G4Material("Target", z=74., a=182.95*g/mole, density=14.69*g/cm3);

    
    // example of vacuum
    density     = 1.e-5*g/cm3;
    pressure    = 1.e-2*bar;
    temperature = 2273.15*kelvin;         //from PhysicalConstants.h
    G4Material* TVac = new G4Material("TVac", density, ncomponents=1,kStateGas,temperature,pressure);
    TVac->AddMaterial(Air, fractionmass=1.);
    
    
    
    
    WorldMaterial = TVac;
    TubeMaterial = Ta;
    BoxMaterial = Graphite;
    DiskMaterial = Target;
    
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
    G4double DiskRmin = 0.*mm;    G4double DiskRmax = 40*mm;    G4double DiskDz   = 1.3*mm;
    G4Tubs* sDisk = new G4Tubs("Disk1",				//name
                               DiskRmin/2,
                               DiskRmax/2,
                               DiskDz/2,
                               0.*deg,
                               360.*deg);	//dimensions
    
    G4LogicalVolume* lDisk = new G4LogicalVolume(sDisk,			//shape
                                                 DiskMaterial,		//material
                                                 "Disk1");		//name
    G4double ZposDisk[7] = {48.73*mm,33.43*mm, 12.13*mm, -9.17*mm, -26.47*mm, -39.77*mm, -53.07*mm};
    
    for(G4int i0=0;i0<7;i0++){
        new G4PVPlacement(0,				//no rotation
                          G4ThreeVector(0.,0.,ZposDisk[i0]),		//at (0,0,0)
                          lDisk,			//logical volume
                          "Disk1",			//name
                          logicWorld,	       		//mother  volume
                          false,			//no boolean operation
                          i0);				//copy number
    }
    
    // Blocks
    G4double Block1Rmin = 0*mm;  G4double Block1Rmax = 45*mm;    G4double Block1Dz= 46.03*mm;
    G4Tubs* sBlock1 = new G4Tubs("Block1",				//name
                                 Block1Rmin/2,
                                 Block1Rmax/2,
                                 Block1Dz/2,
                                 0.*deg,
                                 360.*deg);	//dimensions
    
    G4LogicalVolume* lBlock1 = new G4LogicalVolume(sBlock1,			//shape
                                                   BoxMaterial,		//material
                                                   "Block1");		//name
    
    new G4PVPlacement(0,				//no rotation
                      G4ThreeVector(0.,0.,-81.235),		//at (0,0,0)
                      lBlock1,			//logical volume
                      "Block1",			//
                      logicWorld,	       		//mother  volume
                      false,			//no boolean operation
                      0);				//copy number
    
    
    G4double Block2Rmin = 0*mm;  G4double Block2Rmax = 45*mm;    G4double Block2Dz= 49.37*mm;
    
    G4Tubs* sBlock2 = new G4Tubs("Block2",				//name
                                 Block2Rmin/2,
                                 Block2Rmax/2,
                                 Block2Dz/2,
                                 0.*deg,
                                 360.*deg);	//dimensions
    
    G4LogicalVolume* lBlock2 = new G4LogicalVolume(sBlock2,			//shape
                                                   BoxMaterial,		//material
                                                   "Block2");		//name
    
    
    new G4PVPlacement(0,				//no rotation
                      G4ThreeVector(0.,0.,79.565),		//at (0,0,0)
                      lBlock2,			//logical volume
                      "Block2",			//name
                      logicWorld,	       		//mother  volume
                      false,			//no boolean operation
                      0);				//copy number
    
    
    
    
    // TRANSFER LINE: TRANSFER TORUS + TRANSFER TUBE
    
    G4RotationMatrix* aRot = new G4RotationMatrix();
    aRot->rotateY(180.*deg);
    
    G4RotationMatrix* bRot = new G4RotationMatrix();
    bRot->rotateY(90.*deg);
    
    G4double pRmin    = 8*mm; G4double pRmax    = 8.8*mm;    G4double pRtor    = 8*mm;
    
    G4Torus* sTransfer = new G4Torus("Transfer",				//name
                                     pRmin/2,
                                     pRmax/2,
                                     pRtor,
                                     30.*deg,
                                     60.*deg);	//dimensions
    
    G4LogicalVolume* lTransfer = new G4LogicalVolume(sTransfer,			//shape
                                                     TubeMaterial,		//material
                                                     "Transfer");		//name
    
    
    
    new G4PVPlacement(aRot,				// rotation
                      G4ThreeVector (18.*mm,16.*mm,0.),//4*mm,10*mm,0.),		//at (0,0,0)
                      lTransfer,			//logical volume
                      "Transfer",			//name
                      logicWorld,	       		//mother  volume
                      false,			//no boolean operation
                      0);				//copy number
    
    
    
    G4double TransferRmin = 8*mm; G4double TransferRmax= 8.8*mm;    G4double TransferDz= 71.3*mm;
    
    G4Tubs* sTransfertube = new G4Tubs("Transfertube",				//name
                                       TransferRmin/2,
                                       TransferRmax/2,
                                       TransferDz/2,
                                       0.*deg,
                                       360.*deg);	//dimensions
    
    G4LogicalVolume* lTransfertube = new G4LogicalVolume(sTransfertube,			//shape
                                                         TubeMaterial,		//material
                                                         "Transfertube");		//name
    
    
    
    new G4PVPlacement(bRot,				// rotation
                      G4ThreeVector(18.+35.65*mm,24.*mm,0.),//12*mm,27*mm, 0. ),//40 *mm,16.5*mm, 0.),		//at (0,0,0)
                      lTransfertube,			//logical volume
                      "Transfertube",			//name
                      logicWorld,	       		//mother  volume
                      false,			//no boolean operation
                      0);				//copy number
    
    
    
    
    // "SOLID" HOLE
    
    
    
    G4Torus* sHole = new G4Torus("Hole",				//name
                                 0,pRmin/2,
                                 pRtor,
                                 0.*deg,
                                 90.*deg);	//dimensions
    
    
    G4ThreeVector transtransfer(18.*mm,16.*mm,0.);
    
    
    
    
    // HOLE: TANTALUM TUBE - "SOLID" HOLE
    
    // Tantalium Tube
    G4double TubeRmin  = 49.5*mm;    G4double TubeRmax  =50*mm;    G4double TubeDz    =215.5*mm;
    G4Tubs* sTube = new G4Tubs("Tube",				//name
                               TubeRmin/2,
                               TubeRmax/2,
                               TubeDz/2,
                               0.*deg,
                               360.*deg);	//dimensions
    
    
    
    G4SubtractionSolid* ssubtraction1 = new G4SubtractionSolid("Tube-Transfer",
                                                               sTube,
                                                               sHole,
                                                               aRot,
                                                               transtransfer);
    
    G4LogicalVolume* lsubtraction1 = new G4LogicalVolume(ssubtraction1,			//shape
                                                         TubeMaterial,		//material
                                                         "Tube-Transfer");		//name
    
    new G4PVPlacement(0,				//no rotation
                      G4ThreeVector(),		//at (0,0,0)
                      lsubtraction1,			//logical volume
                      "Tube-Transfer",			//name
                      logicWorld,	       		//mother  volume
                      false,			//no boolean operation
                      0);				//copy number
    
    
    
    
    // HOLE: GRAPHITE BOX - "SOLID" HOLE
    
    // Graphite Box
    G4double BoxRmin  = 45*mm;    G4double BoxRmax  =49.5*mm;    G4double BoxDz    =208.5*mm;
    G4Tubs* sBox = new G4Tubs("Box",				//name
                              BoxRmin/2,
                              BoxRmax/2,
                              BoxDz/2,
                              0.*deg,
                              360.*deg);	//dimensions
    
    
    G4SubtractionSolid* ssubtraction2 = new G4SubtractionSolid("Box-Transfer",
                                                               sBox,
                                                               sHole,
                                                               aRot,
                                                               transtransfer);
    
    G4LogicalVolume* lsubtraction2 = new G4LogicalVolume(ssubtraction2,			//shape
                                                         BoxMaterial,		//material
                                                         "Box-Transfer");		//name
    
    new G4PVPlacement(0,				//no rotation
                      G4ThreeVector(),		//at (0,0,0)
                      lsubtraction2,			//logical volume
                      "Box-Transfer",			//name
                      logicWorld,	       		//mother  volume
                      false,			//no boolean operation
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
                                                  G4ThreeVector(90.*mm,30.*mm,0.),	//at (0,0,0)
                      logicDetector,		//its logical volume
                                                  "detector",		//its name
                                                  logicWorld,			//its mother  volume
                                                  false,			//no boolean operation
                                                  0);			//copy number

    
    
    
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

