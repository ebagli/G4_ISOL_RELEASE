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
#include "TargetSensitiveDetector.hh"

#include "EffusionOptrMultiParticleChangeCrossSection.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction():
fTemperature(1600.*CLHEP::kelvin)
{
    fMessenger = new DetectorConstructionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials(){
    
	// Elements
    G4Material* Graphite = new G4Material("Graphite",
                                          6., //z
                                          12.011*g/mole, //a
                                          1.9*g/cm3, //density
                                          kStateSolid, //state,
                                          fTemperature); //fTemperature
    
	//Uranium carbide
    G4Element* C = G4NistManager::Instance()->FindOrBuildElement("C");
    G4Element* U = G4NistManager::Instance()->FindOrBuildElement("U");

    G4Material* UC4 = new G4Material("UC4",
                                     4.*g/cm3, //density
                                     2, //nComponents
                                     kStateSolid, //state,
                                     fTemperature); //fTemperature
    UC4->AddElement(C,
                        4); //natoms
    UC4->AddElement(U,
                        1); //natoms


	//Tantalum
    G4Material* Tantalum = new G4Material("Tantalum",
                                          73., //z
                                          180.95*g/mole, //a
                                          16.69*g/cm3, //density
                                          kStateSolid, //state,
                                          fTemperature); //fTemperature
    
	//
	//Vacuum
	//
  	// Get nist material manager
	G4Material* Air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
    G4Material* TVac = new G4Material("TVac",
                                      1.e-5*g/cm3, //density
                                      1, //ncomponents
                                      kStateGas, //state
                                      fTemperature, //fTemperature
                                      1.e-2*bar); //pressure
    TVac->AddMaterial(Air,
                      1.); //fractionmass
    
    WorldMaterial = TVac;
    TubeMaterial = Tantalum;
    BoxMaterial = Graphite;
    DiskMaterial = UC4;
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    DefineMaterials();

    G4Colour color_white   (1.0, 1.0, 1.0);
    G4Colour color_gray    (0.5, 0.5, 0.5);
    G4Colour color_black   (0.0, 0.0, 0.0);
    G4Colour color_red     (1.0, 0.0, 0.0);
    G4Colour color_green   (0.0, 0.1, 0.0);
    G4Colour color_blue    (0.0, 0.0, 1.0);
    G4Colour color_cyan    (0.0, 0.1, 1.0);
    G4Colour color_magenta (1.0, 0.0, 1.0);
    G4Colour color_yellow  (1.0, 1.0, 0.0);
    G4Colour color_orange  (1.0, 0.5, 0.0);

    //*********************************************************//
    //
    // Definition of the enclosing box (world)
    //
    G4double WorldSizeX = 100*cm; G4double WorldSizeYZ=100*cm;
    G4Box* solidWorld = new G4Box("World",
                                  WorldSizeX/2,
                                  WorldSizeYZ/2,
                                  WorldSizeYZ/2);
    
    G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld,
                                                      WorldMaterial,
                                                      "World");
    
    G4PVPlacement* physiWorld = new G4PVPlacement(0,
                                                  G4ThreeVector(),
                                                  "World",
                                                  logicWorld,
                                                  0,
                                                  false,
                                                  0);

    G4VisAttributes* worldVisAtt = new G4VisAttributes(color_white);
    worldVisAtt->SetVisibility(false);
    logicWorld->SetVisAttributes(worldVisAtt);

    
    //*********************************************************//
    //
    // Definition of the Target Disks
    //
    
    G4double DiskDmin =  0.0 * mm;
    G4double DiskDmax = 40.0 * mm;
    G4double DiskDz   =  0.8 * mm;

    
    G4bool bTest = true;
    
    if(bTest == true){
        //DiskDmin =  0.0 * mm;
        //DiskDmax = 26.0 * mm;
        DiskDz   =  1.0 * mm;
    }
    

    G4double DiskPosZ[7] = {
         55.68 * mm,
         36.38 * mm,
         10.08 * mm,
        -15.22 * mm,
        -32.52 * mm,
        -49.82 * mm,
        -66.12 * mm};
    
    for(G4int i0=0;i0<7;i0++){

        G4String diskName = "Disk";
        diskName += std::to_string(i0);

        CreateTub(diskName,
                  DiskDmin,
                  DiskDmax,
                  DiskDz,
                  DiskPosZ[i0],
                  DiskMaterial,
                  logicWorld,
                  color_yellow);
    }

    //*********************************************************//

    
    //*********************************************************//
    //
    // Definition of the Target Window
    //

    //Window
    G4double WindowDmin =   0.0  * mm;
    G4double WindowDmax = DiskDmax + 5.0  * mm; // Original 45.0 * mm
    G4double WindowDz   =   0.2  * mm;
    G4double WindowPosz = DiskPosZ[6] - 10. * mm; // Original -77.32 * mm

    if(bTest == true){
        WindowDz   =  0.4 * mm;
    }

    CreateTub("Window",
              WindowDmin,
              WindowDmax,
              WindowDz,
              WindowPosz,
              BoxMaterial,
              logicWorld,
              color_red);

    /*
    //Spacer
    G4double WindowSpacer1Dmin =  DiskDmax + 2.0  * mm; // Original 42.0 * mm
    G4double WindowSpacer1Dmax =  DiskDmax + 5.0  * mm; // Original 45.0 * mm
    G4double WindowSpacer1Dz   =  29.4 * mm; // Original 29.02 * mm - with 29.4 mm it overlaps the Window
    G4double WindowSpacer1Posz =  WindowPosz - WindowSpacer1Dz * 0.5; // Original -92.03 * mm
    
    G4double WindowSpacer2Dmin =  DiskDmax + 1.0  * mm; // Original 41.0  * mm
    G4double WindowSpacer2Dmax =  DiskDmax + 2.0  * mm; // Original 42.0  * mm
    G4double WindowSpacer2Dz   =  12.0 * mm; // Original 11.62 * mm - with 12.0 mm it overlaps the Window
    G4double WindowSpacer2Posz =  WindowPosz + WindowSpacer2Dz * 0.5; // Original -83.33 * mm
    
    CreateTub("WindowSpacer1",
              WindowSpacer1Dmin,
              WindowSpacer1Dmax,
              WindowSpacer1Dz,
              WindowSpacer1Posz,
              BoxMaterial,
              logicWorld,
              color_blue);

    CreateTub("WindowSpacer2",
              WindowSpacer2Dmin,
              WindowSpacer2Dmax,
              WindowSpacer2Dz,
              WindowSpacer2Posz,
              BoxMaterial,
              logicWorld,
              color_blue);
    */
    

    //*********************************************************//
    //
    // Definition of the Target Dumps
    //
    
    // Dump
    G4double DumpD1min  =  0.0  * mm;
    G4double DumpD1max  = 45.0  * mm;
    G4double Dump1Dz    =  0.8  * mm;
    G4double Dump1Posz  = 68.28 * mm;

    G4double DumpD2min  =  0.0  * mm;
    G4double DumpD2max  = 45.0  * mm;
    G4double Dump2Dz    =  0.8  * mm;
    G4double Dump2Posz  = 75.08 * mm;

    G4double DumpD3min  =  0.0  * mm;
    G4double DumpD3max  = 45.0  * mm;
    G4double Dump3Dz    =  1.0  * mm;
    G4double Dump3Posz  = 83.88 * mm;
    
    CreateTub("Dump1",
              DumpD1min,
              DumpD1max,
              Dump1Dz,
              Dump1Posz,
              BoxMaterial,
              logicWorld,
              color_red);

    CreateTub("Dump2",
              DumpD2min,
              DumpD2max,
              Dump2Dz,
              Dump2Posz,
              BoxMaterial,
              logicWorld,
              color_red);

    CreateTub("Dump3",
              DumpD3min,
              DumpD3max,
              Dump3Dz,
              Dump3Posz,
              BoxMaterial,
              logicWorld,
              color_red);

    /*
    // Dump Spacer
    G4double DumpSpacer1Dmin =  42.0  * mm;
    G4double DumpSpacer1Dmax =  45.0  * mm;
    G4double DumpSpacer1Dz   =   6.0  * mm;
    G4double DumpSpacer1Posz =  71.68 * mm;
    
    
    G4double DumpSpacer2Dmin =  42.0  * mm;
    G4double DumpSpacer2Dmax =  45.0  * mm;
    G4double DumpSpacer2Dz   =   7.9  * mm; // with 8.0 mm it overlaps the Dump
    G4double DumpSpacer2Posz =  79.43 * mm; // with 79.48 mm it overlaps the Dump

    G4double DumpSpacer3Dmin =  34.0  * mm;
    G4double DumpSpacer3Dmax =  45.0  * mm;
    G4double DumpSpacer3Dz   =  10.1  * mm;
    G4double DumpSpacer3Posz =  89.43 * mm;

    CreateTub("DumpSpacer1",
              DumpSpacer1Dmin,
              DumpSpacer1Dmax,
              DumpSpacer1Dz,
              DumpSpacer1Posz,
              BoxMaterial,
              logicWorld,
              color_blue);
    
    CreateTub("DumpSpacer2",
              DumpSpacer2Dmin,
              DumpSpacer2Dmax,
              DumpSpacer2Dz,
              DumpSpacer2Posz,
              BoxMaterial,
              logicWorld,
              color_blue);
    
    CreateTub("DumpSpacer3",
              DumpSpacer3Dmin,
              DumpSpacer3Dmax,
              DumpSpacer3Dz,
              DumpSpacer3Posz,
              BoxMaterial,
              logicWorld,
              color_blue);
    */
    
    //*********************************************************//
    //
    // Definition of the Box
    //
    

    // Hole Parameters to be modified
    G4double holeAngle = 90 * deg;
    G4double HolePosz = 0.;
    G4double TransferSect1Dz = 10.0 * mm;

    // Box
    G4double BoxDmin  =   45.0 * mm;
    G4double BoxDmax  =   49.5 * mm;
    G4double BoxDz    =  201.5 * mm;
    G4Tubs*  sBox = new G4Tubs("Box.Solid",
                               BoxDmin * 0.5,
                               BoxDmax * 0.5,
                               BoxDz   * 0.5,
                               0.   * deg,
                               360. * deg);


    // Hole
    //G4double holeTan = std::tan(holeAngle);
    G4double holeCos = std::cos(holeAngle);
    G4double holeSin = std::sin(holeAngle);
    G4double holeCosHalf = std::cos(holeAngle*0.5);
    G4double holeSinHalf = std::sin(holeAngle*0.5);
    G4double holePos = (BoxDmax+BoxDmin) * 0.5 * 0.5; // initial value 15.

    G4double HoleDmin  =  0.0 * mm;
    G4double HoleDmax  =  8.8 * mm;
    G4double HoleDz    = 30.0 * mm;
    
    G4Tubs*  sHole = new G4Tubs("Hole.Solid",
                                HoleDmin,
                                HoleDmax * 0.5,
                                HoleDz,
                                0.   * deg,
                                360. * deg);
    
    
    // Box - Hole
    G4RotationMatrix* sHoleRot = new G4RotationMatrix();
    G4double BoxHolePosz = -6.27 * mm;
    sHoleRot->rotateY(-90. * deg);
    sHoleRot->rotateX(holeAngle);
    G4SubtractionSolid* sBoxHole = new G4SubtractionSolid("BoxHole.Solid",
                                                          sBox,
                                                          sHole,
                                                          sHoleRot,
                                                          G4ThreeVector(holePos*holeCos,
                                                                        holePos*holeSin,
                                                                        -BoxHolePosz+HolePosz));
    
    G4LogicalVolume* lBoxHole = new G4LogicalVolume(sBoxHole,
                                                    BoxMaterial,
                                                    "BoxHole.Logic");
    
    new G4PVPlacement(0,
                      G4ThreeVector(0,0,BoxHolePosz),
                      lBoxHole,
                      "BoxHole.Physical",
                      logicWorld,
                      false,
                      0);
    
    G4VisAttributes* BoxVisAtt = new G4VisAttributes(color_gray);
    BoxVisAtt->SetVisibility(true);
    lBoxHole->SetVisAttributes(BoxVisAtt);

    
    //*********************************************************//
    //
    // Definition of the Transfer Line
    //

    G4double TransferDmin = 8.0 * mm;
    G4double TransferDmax = 8.8 * mm;
    
    // Section 1
    G4CutTubs* sSection1 = new G4CutTubs("Section1.Solid",
                                         TransferDmin * 0.5,
                                         TransferDmax * 0.5,
                                         TransferSect1Dz * 0.5,
                                         0.   * deg,
                                         360. * deg,
                                         G4ThreeVector( 0.0,
                                                        0.0,
                                                       -1.0),
                                         G4ThreeVector( 0.0,
                                                       -holeSinHalf,
                                                        holeCosHalf));
  
    G4LogicalVolume* lTransfer1 = new G4LogicalVolume(sSection1,
                                                      TubeMaterial,
                                                      "Transfer1.Logic");
    
    G4RotationMatrix* Transfer1Rot = new G4RotationMatrix();
    Transfer1Rot->rotateY(-90. * deg);
    Transfer1Rot->rotateX(holeAngle);
   
    
    G4ThreeVector Transfer1Posz = G4ThreeVector((holePos+TransferSect1Dz*0.5)*holeCos,
                                                (holePos+TransferSect1Dz*0.5)*holeSin,
                                                 HolePosz);

    new G4PVPlacement(Transfer1Rot,
                      Transfer1Posz,
                      lTransfer1,
                      "Transfer1.Physical",
                      logicWorld,
                      false,
                      0);


    // Section 2
    G4double TransferSect2Dz = 90. * mm - (TransferSect1Dz + holePos) * holeCos;

    G4CutTubs* sTransfer2 = new G4CutTubs("Transfer2.Solid",
                                          TransferDmin * 0.5,
                                          TransferDmax * 0.5,
                                          TransferSect2Dz * 0.5,
                                          0.   * deg,
                                          360. * deg,
                                          G4ThreeVector( 0.0,
                                                        -holeSinHalf,
                                                        -holeCosHalf),
                                          G4ThreeVector( 0.0,
                                                         0.0,
                                                         1.0));
 
    G4LogicalVolume* lTransfer2 = new G4LogicalVolume(sTransfer2,
                                                      TubeMaterial,
                                                      "Transfer2.Logic");
    G4RotationMatrix* Transfer2Rot = new G4RotationMatrix();
    Transfer2Rot->rotateY(-90.*deg);

    G4ThreeVector Transfer2Posz = G4ThreeVector((TransferSect2Dz*0.5 + (holePos+TransferSect1Dz)*holeCos),
                                                ((holePos+TransferSect1Dz)*holeSin),
                                                  HolePosz);
    

    new G4PVPlacement(Transfer2Rot,
                      Transfer2Posz,
                      lTransfer2,
                      "Transfer2.Physical",
                      logicWorld,
                      false,
                      0);

    
    //*********************************************************//
    //
    // Definition of the Heater
    //
    
    
    // Tantalum Heater
    G4double TubeDmin  =  49.6 * mm;
    G4double TubeDmax  =  50.0 * mm;
    G4double TubeDz    = 169.0 * mm;
    G4double TubePosz  = -21.5 * mm;
    
    G4Tubs* sTaTube = new G4Tubs("TaTube.Solid",
                                 TubeDmin/2,
                                 TubeDmax/2,
                                 TubeDz/2,
                                 0.*deg,
                                 360.*deg);
    
    
    
    G4SubtractionSolid* sTaHeater = new G4SubtractionSolid("TaHeater.Solid",
                                                           sTaTube,
                                                           sHole,
                                                           sHoleRot,
                                                           G4ThreeVector(holePos*holeCos,
                                                                         holePos*holeSin,
                                                                         -TubePosz + HolePosz));

    G4LogicalVolume* lTaHeater = new G4LogicalVolume(sTaHeater,
                                                     TubeMaterial,
                                                     "TaHeater.Logic");
    
    G4ThreeVector TaHeaterPosz = G4ThreeVector(0.,
                                               0.,
                                               TubePosz + HolePosz);

    new G4PVPlacement(0,
                      TaHeaterPosz,
                      lTaHeater,
                      "TaHeater.Physical",
                      logicWorld,
                      false,
                      0);
    
    
   
    

    
   
    
    //*********************************************************//
    //
    // Definition of the Detector
    //
    
    G4double detX  =  1. * mm;
    G4double detYZ = 100. * mm;
    G4ThreeVector detPos = G4ThreeVector( 90.5 * mm,
                                          30.0 * mm,
                                          HolePosz);
    
    G4Box* sDetector = new G4Box("Detector.Solid",
                                 detX  * 0.5,
                                 detYZ * 0.5,
                                 detYZ * 0.5);
    
    G4LogicalVolume* lDetector = new G4LogicalVolume(sDetector,
                                                     WorldMaterial,
                                                     "Detector.Logic");
    
    new G4PVPlacement(0,
                      detPos,
                      lDetector,
                      "Detector.Physical",
                      logicWorld,
                      false,
                      0);
    
    G4VisAttributes* DetectorVisAtt = new G4VisAttributes(color_white);
    DetectorVisAtt->SetVisibility(true);
    lDetector->SetVisAttributes(DetectorVisAtt);

    return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField(){
    
    G4String SDname;
    G4VSensitiveDetector* telescope = new TargetSensitiveDetector(SDname="/telescope",1);
    G4SDManager::GetSDMpointer()->AddNewDetector(telescope);
    
    G4LogicalVolume* fDetLogic = G4LogicalVolumeStore::GetInstance()->GetVolume("Detector.Logic");
    if(fDetLogic!=NULL){
        fDetLogic->SetSensitiveDetector(telescope);
        G4cout << "--- Attaching sensitive detector " << telescope->GetName()
        << " to logical volume " << fDetLogic->GetName() << G4endl;
    }
    
    G4VSensitiveDetector* ucxdet = new TargetSensitiveDetector(SDname="/ucx",0);
    G4SDManager::GetSDMpointer()->AddNewDetector(ucxdet);

    for(G4int i0=0;i0<7;i0++){
        
        G4String diskName = "Disk";
        diskName += std::to_string(i0);
        diskName += ".Logic";
        G4LogicalVolume* detLogic = G4LogicalVolumeStore::GetInstance()->GetVolume(diskName);
        if(detLogic!=NULL){
            detLogic->SetSensitiveDetector(ucxdet);
            G4cout << "--- Attaching sensitive detector " << ucxdet->GetName()
            << " to logical volume " << diskName << G4endl;

        }
    }
    
    EffusionOptrMultiParticleChangeCrossSection* effusionXSchange = new EffusionOptrMultiParticleChangeCrossSection();
    effusionXSchange->AddParticle("GenericIon");
    // Modify Radioactive In-Flight Decay with Sticking Time
    for (auto lv : *G4LogicalVolumeStore::GetInstance()){
        G4String lvName = lv->GetName();
        effusionXSchange->AttachTo(lv);
        G4cout << "--- Attaching biasing operator " << effusionXSchange->GetName()
        << " to logical volume " << lvName << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::CreateTub(G4String name,
                                     G4double Dmin,
                                     G4double Dmax,
                                     G4double Dz,
                                     G4double PosZ,
                                     G4Material* material,
                                     G4LogicalVolume* motherVolume,
                                     G4Colour color){

    G4String sname = name + ".Solid";
    G4Tubs* solid = new G4Tubs(sname,
                               Dmin * 0.5,
                               Dmax * 0.5,
                               Dz   * 0.5,
                               0.   * deg,
                               360. * deg);
    
    G4String lname = name + ".Logic";
    G4LogicalVolume* logic = new G4LogicalVolume(solid,
                                                 material,
                                                 lname);
    G4String pname = name + ".Physical";

    new G4PVPlacement(0,
                      G4ThreeVector(0.,0.,PosZ),
                      logic,
                      pname,
                      motherVolume,
                      false,
                      0);

    G4VisAttributes* visAtt = new G4VisAttributes(color);
    visAtt->SetVisibility(true);
    logic->SetVisAttributes(visAtt);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

