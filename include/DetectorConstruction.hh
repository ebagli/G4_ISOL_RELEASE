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
// $Id: DetectorConstruction.hh,v 1.6 2006-06-29 16:54:31 gunter Exp $
// GEANT4 tag $Name: geant4-09-04-patch-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Colour.hh"
#include "DetectorConstructionMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction();
    ~DetectorConstruction();
    
    G4VPhysicalVolume* Construct();
    void DefineMaterials();
    void ConstructSDandField();

private:
    G4Material* WorldMaterial;
    G4Material* TubeMaterial;
    G4Material* BoxMaterial;
    G4Material* DiskMaterial;

private:
    void CreateTub(G4String,
                   G4double,
                   G4double,
                   G4double,
                   G4double,
                   G4Material*,
                   G4LogicalVolume*,
                   G4Colour);

private:
    G4double fTemperature;
    
public:
    void SetTemperature(G4double aDouble) {fTemperature=aDouble;}
    G4double GetTemperature() {return fTemperature;}

private:
    DetectorConstructionMessenger* fMessenger;
    
private:
    G4bool bPrimaries;
public:
    void SetPrimaries(G4bool aBool) {bPrimaries=aBool;}
    G4bool GetPrimaries() {return bPrimaries;}

    /* Variables To Be Changed Via Messenger */

private:
    G4String fTargetMaterialName;
public:
    void SetTargetMaterialName(G4String aString) {fTargetMaterialName=aString;}
    G4String GetTargetMaterialName() {return fTargetMaterialName;}

private:
    G4int fTargetDiskNumber;
public:
    void SetTargetDiskNumber(G4int aInt) {fTargetDiskNumber=aInt;}
    G4int GetTargetDiskNumber() {return fTargetDiskNumber;}

private:
    G4double fTargetDensity;
public:
    void SetTargetDensity(G4double aDouble) {fTargetDensity=aDouble;}
    G4double GetTargetDensity() {return fTargetDensity;}

private:
    G4double fTargetBoxInit;
public:
    void SetTargetBoxInit(G4double aDouble) {fTargetBoxInit=aDouble;}
    G4double GetTargetBoxInit() {return fTargetBoxInit;}

private:
    G4double fTargetBoxEnd;
public:
    void SetTargetBoxEnd(G4double aDouble) {fTargetBoxEnd=aDouble;}
    G4double GetTargetBoxEnd() {return fTargetBoxEnd;}

private:
    G4double fTargetDiskRadius;
public:
    void SetTargetDiskRadius(G4double aDouble) {fTargetDiskRadius=aDouble;}
    G4double GetTargetDiskRadius() {return fTargetDiskRadius;}

private:
    G4double fTargetDiskPosition[MAX_DISK_NUMBER];
public:
    void SetTargetDiskPosition(G4int aInt, G4double aDouble) {fTargetDiskPosition[aInt]=aDouble;}
    G4double GetTargetDiskPosition(G4int aInt) {return fTargetDiskPosition[aInt];}

private:
    G4double fTargetDiskThickness[MAX_DISK_NUMBER];
public:
    void SetTargetDiskThickness(G4int aInt, G4double aDouble) {fTargetDiskThickness[aInt]=aDouble;}
    G4double GetTargetDiskThickness(G4int aInt) {return fTargetDiskThickness[aInt];}
    





};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif


