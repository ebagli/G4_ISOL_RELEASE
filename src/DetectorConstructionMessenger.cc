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

#include "DetectorConstructionMessenger.hh"
#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4RunManager.hh"

#include "G4ios.hh"

DetectorConstructionMessenger::
DetectorConstructionMessenger(
                              DetectorConstruction* mpga)
:fTarget(mpga){
    fMyDetDirectory = new G4UIdirectory("/det/");
    fMyDetDirectory->SetGuidance("Detector setup control commands.");
    
    
    fTempCmd = new G4UIcmdWithADoubleAndUnit("/det/setTemperature",this);
    fTempCmd->SetGuidance("Set detector temperature.");
    fTempCmd->SetParameterName("temp",
                               true);
    fTempCmd->SetDefaultValue(273.);
    fTempCmd->SetDefaultUnit("kelvin");

    fAdsTimeCmd = new G4UIcmdWithADoubleAndUnit("/det/setAdsorptionTime",this);
    fAdsTimeCmd->SetGuidance("Set adsorpion time.");
    fAdsTimeCmd->SetParameterName("adstime",
                               true);
    fAdsTimeCmd->SetDefaultValue(0.);
    fAdsTimeCmd->SetDefaultUnit("second");


    fDiffLengthCmd = new G4UIcmdWithADoubleAndUnit("/det/setDiffusionCoefficient",this);
    fDiffLengthCmd->SetGuidance("Set diffusion coefficient [cm2/s]");
    fDiffLengthCmd->SetParameterName("diffcoeff",
                                  true);
    fDiffLengthCmd->SetDefaultValue(0.);
    fDiffLengthCmd->SetDefaultUnit("cm2");


    fPDiffLengthCmd = new G4UIcmdWithADoubleAndUnit("/det/setPorousDiffusionCoefficient",this);
    fPDiffLengthCmd->SetGuidance("Set porous diffusion coefficient [cm2/s]");
    fPDiffLengthCmd->SetParameterName("porousdiffcoeff",
                                     true);
    fPDiffLengthCmd->SetDefaultValue(0.);
    fPDiffLengthCmd->SetDefaultUnit("cm2");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstructionMessenger::
~DetectorConstructionMessenger(){
    delete fTempCmd;
    delete fAdsTimeCmd;
    delete fDiffLengthCmd;
    delete fPDiffLengthCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstructionMessenger::SetNewValue(G4UIcommand *command,
                                                G4String newValue){
    if(command==fTempCmd ){
        fTarget->SetTemperature(fTempCmd->GetNewDoubleValue(newValue));
    }
    if(command==fAdsTimeCmd ){
        fTarget->SetAdsorptionTime(fAdsTimeCmd->GetNewDoubleValue(newValue));
    }
    if(command==fDiffLengthCmd ){
        fTarget->SetDiffusionCoefficient(fDiffLengthCmd->GetNewDoubleValue(newValue) / CLHEP::second);
    }
    if(command==fPDiffLengthCmd ){
        fTarget->SetPorousDiffusionCoefficient(fPDiffLengthCmd->GetNewDoubleValue(newValue) / CLHEP::second);
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4String DetectorConstructionMessenger::GetCurrentValue(G4UIcommand * command){
    G4String cv;
    
    if( command==fTempCmd ){
        cv = fTempCmd->ConvertToString(fTarget->GetTemperature(),"kelvin");
    }
    if( command==fAdsTimeCmd ){
        cv = fAdsTimeCmd->ConvertToString(fTarget->GetAdsorptionTime(),"second");
    }
    if( command==fDiffLengthCmd ){
        cv = fDiffLengthCmd->ConvertToString(fTarget->GetDiffusionCoefficient(),"m");
    }
    if( command==fPDiffLengthCmd ){
        cv = fPDiffLengthCmd->ConvertToString(fTarget->GetPorousDiffusionCoefficient(),"m");
    }

    return cv;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
