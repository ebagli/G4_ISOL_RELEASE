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

    fTargetMaterialCmd = new G4UIcmdWithAString("/det/setTargetMaterial",this);
    fTargetMaterialCmd->SetGuidance("Set target material.");
    fTargetMaterialCmd->SetParameterName("targetmat",
                                         true);
    fTargetMaterialCmd->SetDefaultValue("UC4");

    fTargetDiskNumberCmd = new G4UIcmdWithADouble("/det/setNumberOfDisks",this);
    fTargetDiskNumberCmd->SetGuidance("Set number of disks.");
    fTargetDiskNumberCmd->SetParameterName("numberofdisks",
                                        true);
    fTargetDiskNumberCmd->SetDefaultValue(7);

    fTargetDensityCmd = new G4UIcmdWithADoubleAndUnit("/det/setTargetDensity",this);
    fTargetDensityCmd->SetGuidance("Set target density.");
    fTargetDensityCmd->SetParameterName("targetdensity",
                                        true);
    fTargetDensityCmd->SetDefaultValue(4.);
    fTargetDensityCmd->SetDefaultUnit("g/cm3");

    fTargetBoxInitCmd = new G4UIcmdWithADoubleAndUnit("/det/setBoxInit",this);
    fTargetBoxInitCmd->SetGuidance("Set target box init.");
    fTargetBoxInitCmd->SetParameterName("targetboxinit",
                                        true);
    fTargetBoxInitCmd->SetDefaultValue(10.702);
    fTargetBoxInitCmd->SetDefaultUnit("cm");

    fTargetBoxEndCmd = new G4UIcmdWithADoubleAndUnit("/det/setBoxEnd",this);
    fTargetBoxEndCmd->SetGuidance("Set target box end.");
    fTargetBoxEndCmd->SetParameterName("targetboxinit",
                                        true);
    fTargetBoxEndCmd->SetDefaultValue(9.448);
    fTargetBoxEndCmd->SetDefaultUnit("cm");

    fTargetDiskRadiusCmd = new G4UIcmdWithADoubleAndUnit("/det/setTargetRadius",this);
    fTargetDiskRadiusCmd->SetGuidance("Set target radius.");
    fTargetDiskRadiusCmd->SetParameterName("targetradius",
                                       true);
    fTargetDiskRadiusCmd->SetDefaultValue(2.);
    fTargetDiskRadiusCmd->SetDefaultUnit("cm");

    G4double defaultDistances[MAX_DISK_NUMBER] = {-6.682,-5.052,-3.322,-1.592, +0.938,+3.568,+5.498,0.,0.,0.,
                                                  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

    for(int i = 0;i<MAX_DISK_NUMBER;i++){
        G4String command = "/det/setDiskPosition" + std::to_string(i+1);
        fTargetDiskPositionCmd[i] = new G4UIcmdWithADoubleAndUnit(command,this);
        fTargetDiskPositionCmd[i]->SetGuidance("Set disk position.");
        fTargetDiskPositionCmd[i]->SetParameterName("detSize",
                                                    true);
        fTargetDiskPositionCmd[i]->SetDefaultValue(defaultDistances[i]);
        fTargetDiskPositionCmd[i]->SetDefaultUnit("cm");

        command = "/det/setDiskThickness" + std::to_string(i+1);
        fTargetDiskThicknessCmd[i] = new G4UIcmdWithADoubleAndUnit(command,this);
        fTargetDiskThicknessCmd[i]->SetGuidance("Set disk position.");
        fTargetDiskThicknessCmd[i]->SetParameterName("detSize",
                                                    true);
        fTargetDiskThicknessCmd[i]->SetDefaultValue(0.08);
        fTargetDiskThicknessCmd[i]->SetDefaultUnit("cm");
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorConstructionMessenger::
~DetectorConstructionMessenger(){
    delete fTempCmd;
    delete fTargetMaterialCmd;
    delete fTargetDiskNumberCmd;

    delete fTargetDensityCmd;
    delete fTargetBoxInitCmd;
    delete fTargetBoxEndCmd;
    delete fTargetDiskRadiusCmd;
    
    for(int i = 0;i<MAX_DISK_NUMBER;i++){
        delete fTargetDiskPositionCmd[i];
        delete fTargetDiskThicknessCmd[i];
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorConstructionMessenger::SetNewValue(G4UIcommand *command,
                                                G4String newValue){
    if(command==fTempCmd ){
        fTarget->SetTemperature(fTempCmd->GetNewDoubleValue(newValue));
    }

    if(command==fTargetMaterialCmd ){
        fTarget->SetTargetMaterialName(newValue);
    }

    if(command==fTargetDiskNumberCmd ){
        fTarget->SetTargetDiskNumber(int(fTargetDiskNumberCmd->GetNewDoubleValue(newValue)));
    }

    if(command==fTargetDensityCmd ){
        fTarget->SetTargetDensity(fTargetDensityCmd->GetNewDoubleValue(newValue));
    }

    if(command==fTargetDensityCmd ){
        fTarget->SetTemperature(fTargetDensityCmd->GetNewDoubleValue(newValue));
    }

    if(command==fTargetBoxInitCmd ){
        fTarget->SetTargetBoxInit(fTargetBoxInitCmd->GetNewDoubleValue(newValue));
    }

    if(command==fTargetBoxEndCmd ){
        fTarget->SetTargetBoxEnd(fTargetBoxEndCmd->GetNewDoubleValue(newValue));
    }

    if(command==fTargetDiskRadiusCmd ){
        fTarget->SetTargetDiskRadius(fTargetDiskRadiusCmd->GetNewDoubleValue(newValue));
    }
    
    for(int i = 0;i<MAX_DISK_NUMBER;i++){
        if(command==fTargetDiskPositionCmd[i]){
            fTarget->SetTargetDiskPosition(i,fTargetDiskPositionCmd[i]->GetNewDoubleValue(newValue));
        }
        if(command==fTargetDiskThicknessCmd[i]){
            fTarget->SetTargetDiskThickness(i,fTargetDiskThicknessCmd[i]->GetNewDoubleValue(newValue));
        }
   }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4String DetectorConstructionMessenger::GetCurrentValue(G4UIcommand * command){
    G4String cv;
    
    if( command==fTempCmd ){
        cv = fTempCmd->ConvertToString(fTarget->GetTemperature(),"kelvin");
    }
    if( command==fTargetMaterialCmd ){
        cv = fTarget->GetTargetMaterialName();
    }
    if( command==fTargetDiskNumberCmd ){
        cv = fTempCmd->ConvertToString(fTarget->GetTargetDiskNumber());
    }
    if( command==fTargetDensityCmd ){
        cv = fTargetDensityCmd->ConvertToString(fTarget->GetTargetDensity(),"g/cm3");
    }
    if( command==fTargetBoxInitCmd ){
        cv = fTargetBoxInitCmd->ConvertToString(fTarget->GetTargetBoxInit(),"cm");
    }
    if( command==fTargetBoxEndCmd ){
        cv = fTargetBoxEndCmd->ConvertToString(fTarget->GetTargetBoxEnd(),"cm");
    }
    if( command==fTargetDiskRadiusCmd ){
        cv = fTargetDiskRadiusCmd->ConvertToString(fTarget->GetTargetDiskRadius(),"cm");
    }
    for(int i = 0;i<MAX_DISK_NUMBER;i++){
        if( command==fTargetDiskPositionCmd[i] ){
            cv = fTargetDiskPositionCmd[i]->ConvertToString(fTarget->GetTargetDiskPosition(i),"cm");
        }
        if( command==fTargetDiskThicknessCmd[i] ){
            cv = fTargetDiskThicknessCmd[i]->ConvertToString(fTarget->GetTargetDiskThickness(i),"cm");
        }
    }

    return cv;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
