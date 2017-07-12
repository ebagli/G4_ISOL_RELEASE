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
#include "G4ExtendedMaterial.hh"
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
    G4ExtendedMaterial* TubeMaterial;
    G4ExtendedMaterial* BoxMaterial;
    G4ExtendedMaterial* DiskMaterial;

private:
    void CreateTub(G4String,
                   G4double,
                   G4double,
                   G4double,
                   G4double,
                   G4ExtendedMaterial*,
                   G4LogicalVolume*,
                   G4Colour);

private:
    G4double fTemperature;
    G4double fAdsorptionTime;
    G4double fDiffusionCoefficient;
    G4double fPorousDiffusionCoefficient;
    
public:
    void SetTemperature(G4double aDouble) {fTemperature=aDouble;}
    G4double GetTemperature() {return fTemperature;}

    void SetAdsorptionTime(G4double aDouble) {fAdsorptionTime=aDouble;}
    G4double GetAdsorptionTime() {return fAdsorptionTime;}
    
    void SetDiffusionCoefficient(G4double aDouble) {fDiffusionCoefficient=aDouble;}
    G4double GetDiffusionCoefficient() {return fDiffusionCoefficient;}
    
    void SetPorousDiffusionCoefficient(G4double aDouble) {fPorousDiffusionCoefficient=aDouble;}
    G4double GetPorousDiffusionCoefficient() {return fPorousDiffusionCoefficient;}

private:
    DetectorConstructionMessenger* fMessenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif


