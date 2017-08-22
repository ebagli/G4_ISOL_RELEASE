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

#include "G4ios.hh"
#include "DiffusionProcess.hh"

#include "globals.hh"
#include "geomdefs.hh"

#include "G4GeometryTolerance.hh"
#include "G4ParticleChange.hh"
#include "G4VParticleChange.hh"
#include "G4SystemOfUnits.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "G4RandomDirection.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DiffusionProcess::DiffusionProcess(const G4String& processName): G4VDiscreteProcess(processName){
    kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
    fEffusionID = G4PhysicsModelCatalog::GetIndex("effusion");
    if(fEffusionID == -1){
        fEffusionID = G4PhysicsModelCatalog::Register("effusion");
    }
    /*
    for(auto i = 0;i< 100;i++){
        theDiffusionCoefficientMap.insert({GetIndex(i,"TVac"),0.});
        thePorousDiffusionCoefficientMap.insert({GetIndex(i,"TVac"),1.E14 * CLHEP::cm2/CLHEP::second});
    }
    */
    
    theDiffusionCoefficientMap.insert({GetIndex(37,"UC4"),7.943E-13 * CLHEP::cm2/CLHEP::second});
    thePorousDiffusionCoefficientMap.insert({GetIndex(37,"UC4"),1.E10 * CLHEP::cm2/CLHEP::second});

    theDiffusionCoefficientMap.insert({GetIndex(36,"UC4"),7.943E-13 * CLHEP::cm2/CLHEP::second});
    thePorousDiffusionCoefficientMap.insert({GetIndex(36,"UC4"),1.E10 * CLHEP::cm2/CLHEP::second});

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DiffusionProcess::~DiffusionProcess(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EffusionTrackData* DiffusionProcess::GetTrackData(const G4Track& aTrack){
    EffusionTrackData* trackdata =
    (EffusionTrackData*)(aTrack.GetAuxiliaryTrackInformation(fEffusionID));
    if(trackdata == nullptr){
        trackdata = new EffusionTrackData();
        aTrack.SetAuxiliaryTrackInformation(fEffusionID,trackdata);
    }
    return trackdata;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VParticleChange*
DiffusionProcess::PostStepDoIt(const G4Track& aTrack, const G4Step&)
{
    aParticleChange.Initialize(aTrack);

    if(aTrack.GetTrackID()==1 && aTrack.GetCurrentStepNumber()!=1) {
        G4double diff_coeff0  = GetDiffusionCoefficient(aTrack);//cm2/s

        if(diff_coeff0 == 0.){
            return &aParticleChange;
        }
        
        G4double R = 1.9872036E-3;// kcal/mol/K;
        G4double activation_energy = 56.4;// kcal/mol/K;
        G4double a_mean = 1.E-2 * CLHEP::nanometer;
        G4double a_sigma = 1.E-3 * CLHEP::nanometer;
        G4double a = G4RandGauss::shoot(a_mean,a_sigma);

        G4double T = aTrack.GetVolume()->GetLogicalVolume()->GetMaterial()->GetTemperature();
        
        G4double diff_coeff  = diff_coeff0 * exp(-activation_energy/R/T);//cm2/s
        
        // exact solution Fick's equation for a sphere of radius "a"
        // Fujioka, NIM 186, 409 (1981)
        G4double tau = a * a / diff_coeff;
        
        aParticleChange.ProposeGlobalTime(aTrack.GetGlobalTime() + tau ) ;
        GetTrackData(aTrack)->SetTimeSticked(tau);
    }
    
    if(aTrack.GetStepLength()<=kCarTolerance/2){
        return &aParticleChange;
    }
    
    G4ThreeVector newDir = G4RandomDirection();
    aParticleChange.ProposeMomentumDirection(newDir);
    
    return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DiffusionProcess::GetMeanFreePath(const G4Track& aTrack,
                                          G4double ,
                                          G4ForceCondition* condition)
{
    if(bDONT_USE_DIFFUSION==true){
        *condition = Forced;
        return DBL_MAX;
    }
    
    G4double theMFP = DBL_MAX;
    
    G4double diff_coeff0  = GetPorousDiffusionCoefficient(aTrack);
    if(diff_coeff0 == DBL_MAX){
        theMFP = DBL_MAX;
    }
    else{
        // mean free path from diffusion coefficient
        G4double R = 1.9872036E-3;// kcal/mol/K;
        G4double activation_energy = 56.4;// kcal/mol/K;
        G4double T = aTrack.GetVolume()->GetLogicalVolume()->GetMaterial()->GetTemperature();
        G4double diff_coeff  = diff_coeff0 * exp(-activation_energy/R/T);//cm2/s
        G4double velocity = aTrack.GetVelocity();
        theMFP = 2. * diff_coeff / velocity;
    }
    
    G4double minMFP = 0.001 * CLHEP::mm;
    if(theMFP < minMFP){
        theMFP =  minMFP;
    }
    
    return theMFP;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DiffusionProcess::GetDiffusionCoefficient(const G4Track& aTrack)
{
    G4Material* mat = aTrack.GetVolume()->GetLogicalVolume()->GetMaterial();
    const std::string matName = mat->GetName();
    int partZ = int(std::round(aTrack.GetDefinition()->GetPDGCharge()));

    std::unordered_map<int, double>::iterator it =
        theDiffusionCoefficientMap.find(GetIndex(partZ,matName));
    
    G4double diffusionCoefficient = 0.;
    
    if (it != theDiffusionCoefficientMap.end()){
        diffusionCoefficient = it->second;
    }

    return diffusionCoefficient;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DiffusionProcess::GetPorousDiffusionCoefficient(const G4Track& aTrack)
{
    G4Material* mat = aTrack.GetVolume()->GetLogicalVolume()->GetMaterial();
    const std::string matName = mat->GetName();
    int partZ = int(std::round(aTrack.GetDefinition()->GetPDGCharge()));
    
    std::unordered_map<int, double>::iterator it =
    thePorousDiffusionCoefficientMap.find(GetIndex(partZ,matName));
    
    G4double porousDiffusionCoefficient = DBL_MAX;
    
    if (it != thePorousDiffusionCoefficientMap.end()){
        porousDiffusionCoefficient = it->second;
    }
    
    return porousDiffusionCoefficient;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
