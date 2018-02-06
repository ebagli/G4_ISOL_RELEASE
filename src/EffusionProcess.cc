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
#include "EffusionProcess.hh"

#include "globals.hh"
#include "geomdefs.hh"

#include "G4GeometryTolerance.hh"
#include "G4ParticleChange.hh"
#include "G4VParticleChange.hh"
#include "G4SystemOfUnits.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "G4RandomTools.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EffusionProcess::EffusionProcess(const G4String& processName)
: G4VDiscreteProcess(processName),
theSampledKineticEnergy(0.),
theLocalNormal(G4ThreeVector()),
theLocalPoint(G4ThreeVector()),
theGlobalNormal(G4ThreeVector()),
theGlobalPoint(G4ThreeVector()),
validLocalNorm(false){
    kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
    fEffusionID = G4PhysicsModelCatalog::GetIndex("effusion");
    if(fEffusionID == -1){
        fEffusionID = G4PhysicsModelCatalog::Register("effusion");
    }
    
    SetAdsorptionTime(37,73,100);
    SetAdsorptionTime(37,6,4420);

    SetAdsorptionTime(36,73,100);
    SetAdsorptionTime(36,6,4420);

    fAdsorptionTimeMessenger =
    new G4GenericMessenger(this,
                           "/effusion/",
                           "Change Adsorption Time" );
    fAdsorptionTimeMessenger->DeclareMethod("loadAdsTime", &EffusionProcess::LoadAdsorptionTime,
                                            "load adsorption time partZ;matZ;time_ns" );
    fAdsorptionTimeMessenger->SetGuidance("particle_Z;material_Z;time_ns");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EffusionTrackData* EffusionProcess::GetTrackData(const G4Track& aTrack){
    EffusionTrackData* trackdata =
    (EffusionTrackData*)(aTrack.GetAuxiliaryTrackInformation(fEffusionID));
    if(trackdata == nullptr){
        trackdata = new EffusionTrackData();
        aTrack.SetAuxiliaryTrackInformation(fEffusionID,trackdata);
    }
    return trackdata;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EffusionProcess::~EffusionProcess(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VParticleChange*
EffusionProcess::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
    aParticleChange.Initialize(aTrack);

    if(aTrack.GetGlobalTime() > 360000. * CLHEP::second) {
        G4Exception("EffusionProcess::PostStepDoIt",
                    "eff0001",
                    JustWarning,
                    "Particle killed after 100 hours.");
        
        aParticleChange.ProposeEnergy(0.);
        aParticleChange.ProposeTrackStatus(fStopAndKill);
        return &aParticleChange;
    }
    
    // Check Boundaries
    G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
    G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();
    
    if(pPostStepPoint->GetStepStatus() != fGeomBoundary){
        return &aParticleChange;
    }

    // Check StepLength
    if(aTrack.GetStepLength()<=kCarTolerance/2){
        return &aParticleChange;
    }

    // Check Materials of next and previous volumes, if same do nothing.
    G4Material* aMaterialPre  = pPreStepPoint ->GetPhysicalVolume()->GetLogicalVolume()->GetMaterial();
    G4Material* aMaterialPost = pPostStepPoint->GetPhysicalVolume()->GetLogicalVolume()->GetMaterial();
    
    if(aMaterialPre == aMaterialPost){
        return &aParticleChange;
    }

    ///////////////////////////////////////////////////////////////////////////////////
    if(pPostStepPoint->GetPhysicalVolume()->GetMotherLogical() == 0){
        return &aParticleChange;
    }
    
    // If the particle is diffused into the material, it continues its motion
    if(G4UniformRand() < GetDiffusionProbability(aTrack)){
        return &aParticleChange;
    }
    
    // If the particle is adsorbed, the average adsorption time is summed to the
    // particle global time, i.e., the particle re-start to travel after a time
    // equal to the average adsorption time
    if(G4UniformRand() < GetAdsorptionProbability(aTrack)){
        
        G4double adsorptionTime = GetAdsorptionTime(aTrack);
        aParticleChange.ProposeGlobalTime(aTrack.GetGlobalTime() + adsorptionTime ) ;
        GetTrackData(aTrack)->SetTimeSticked(adsorptionTime);

        // If the particle is adsorbed and not released, it is killed
        if(G4UniformRand() < GetFullAdsorptionProbability(aTrack)){
            aParticleChange.ProposeEnergy(0.);
            aParticleChange.ProposeTrackStatus(fStopAndKill);
            return &aParticleChange;
        }
        
        // Get a random kinetic energy following the Maxwell-Boltzmann distribution
        // theSampledKineticEnergy is computed by the SampleMaxwellBoltzmannKineticEnergy function
        //SampleMaxwellBoltzmannKineticEnergy(aTrack);
        //aParticleChange.ProposeEnergy(theSampledKineticEnergy);

        // Compute the outgoing particle direction following the cosine-law (Lambertian) distribution
        // theGlobalNormal is computed by the SampleLambertianDirection() function
        SampleLambertianDirection(aTrack);
        aParticleChange.ProposeMomentumDirection(G4LambertianRand(theGlobalNormal));
 
        return &aParticleChange;
    }

    // Compute the outgoing particle direction following the cosine-law (Lambertian) distribution
    // theGlobalNormal is computed by the SampleLambertianDirection() function
    SampleLambertianDirection(aTrack);
    aParticleChange.ProposeMomentumDirection(G4LambertianRand(theGlobalNormal));
    
    return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EffusionProcess::SampleMaxwellBoltzmannKineticEnergy(const G4Track& aTrack){
    // Ideal Gas case : Maxwell Boltzmann Distribution for Kinetic Energy
    
    G4StepPoint* pPostStepPoint = aTrack.GetStep()->GetPostStepPoint();
    G4Material* aMaterialPost = pPostStepPoint->GetPhysicalVolume()->GetLogicalVolume()->GetMaterial();
    
    sqrtkTm = std::sqrt(CLHEP::k_Boltzmann * aMaterialPost->GetTemperature()*aTrack.GetDefinition()->GetPDGMass());
    
    thePx = G4RandGauss::shoot(0.,sqrtkTm);
    thePy = G4RandGauss::shoot(0.,sqrtkTm);
    thePz = G4RandGauss::shoot(0.,sqrtkTm);

    theSampledKineticEnergy = (thePx*thePx+thePy*thePy+thePz*thePz) * 0.5 / aTrack.GetDefinition()->GetPDGMass();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EffusionProcess::SampleLambertianDirection(const G4Track& aTrack){

    // Get the Global Point
    G4StepPoint* pPostStepPoint = aTrack.GetStep()->GetPostStepPoint();
    theGlobalPoint = pPostStepPoint->GetPosition();
 
    // Get Transport Navigator
    G4Navigator* theNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
    
    // Get the Local Point
    theLocalPoint = theNavigator->GetGlobalToLocalTransform().TransformPoint(theGlobalPoint);
    
    // Get the Local Normal : Normal points back and check if it is possible to get it
    validLocalNorm = false;
    theLocalNormal = theNavigator->GetLocalExitNormal(&validLocalNorm);
    
    // Point the Local Normal back into volume
    if (validLocalNorm == true) {
        theLocalNormal = -theLocalNormal;
    }
    // Transform the Local to the Global Normal
    theGlobalNormal = theNavigator->GetLocalToGlobalTransform().TransformAxis(theLocalNormal);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double EffusionProcess::GetMeanFreePath(const G4Track&,
                                          G4double,
                                          G4ForceCondition* condition)
{
    *condition = Forced;
    return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double EffusionProcess::GetDiffusionProbability(const G4Track&)
{
    return 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double EffusionProcess::GetAdsorptionProbability(const G4Track&)
{
    return 1.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double EffusionProcess::GetAdsorptionTime(const G4Track& aTrack)
{
    G4Material* mat = aTrack.GetVolume()->GetLogicalVolume()->GetMaterial();

    const G4ElementVector* theElementVector = mat->GetElementVector();
    const G4Element* element = (*theElementVector)[0] ;
  
    int matZ  = int(std::round(element->GetZ()));
    int partZ = int(std::round(aTrack.GetDefinition()->GetPDGCharge()));
    
    std::unordered_map<int, double>::iterator it =
        theAdsorptionTimeMap.find(GetIndex(partZ,matZ));
    
    G4double adsorptionTime = 0.;
    
    if (it != theAdsorptionTimeMap.end()){
        adsorptionTime = it->second;
    }
    
    return adsorptionTime;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double EffusionProcess::GetFullAdsorptionProbability(const G4Track&)
{
    return 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......





