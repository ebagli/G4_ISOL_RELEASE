// ------------------------------------------------------------
// Brahim Apr 10th 2003 : Diffusion Process source file for ganilN01
//
// This file is based on lN01DiffusionProcess.cc of li8N01
//-------------------------------------------------------------
////////////////////////////////////////////////////////////////////////
// Particle Boundary Process Class Implementation
////////////////////////////////////////////////////////////////////////
// File:        DiffusionProcess.cc
// Description: Discrete Process -- Adsorption/Desorption of Particles
////////////////////////////////////////////////////////////////////////

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

DiffusionProcess::DiffusionProcess(const G4String& processName)
: G4VDiscreteProcess(processName){
    kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DiffusionProcess::~DiffusionProcess(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VParticleChange*
DiffusionProcess::PostStepDoIt(const G4Track& aTrack, const G4Step&)
{
    aParticleChange.Initialize(aTrack);

    // Check StepLength
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
                                          G4ForceCondition*)
{
    G4double theMFPforPorous = 0.1 *CLHEP::mm;
    G4double theMFPforCrystal = 0.1 *CLHEP::mm;
    G4double thePorosity = 0.5;
    G4bool aDiffusiveMaterial = false;

    G4StepPoint* pPreStepPoint  = aTrack.GetStep()->GetPreStepPoint();
    G4Material* aMaterialPre  = pPreStepPoint ->GetPhysicalVolume()->GetLogicalVolume()->GetMaterial();
    
    
    if(aMaterialPre->GetName() == "Tantalum"){
        theMFPforPorous = 0.1 *CLHEP::mm;
        theMFPforCrystal = 0.1 *CLHEP::mm;
        thePorosity = 0.0;
        aDiffusiveMaterial = true;
    }
    
    if(aMaterialPre->GetName() == "Graphite"){
        theMFPforPorous = 0.1 *CLHEP::mm;
        theMFPforCrystal = 0.1 *CLHEP::mm;
        thePorosity = 0.0;
        aDiffusiveMaterial = true;
    }

    if(aMaterialPre->GetName() == "Target"){
        theMFPforPorous = 0.1 *CLHEP::mm;
        theMFPforCrystal = 0.1 *CLHEP::mm;
        thePorosity = 0.0;
        aDiffusiveMaterial = true;
    }

    if(aDiffusiveMaterial == true){
        if(G4UniformRand() < thePorosity){
            return theMFPforPorous;
        }
        else{
            return theMFPforCrystal;
        }
    }

    return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
