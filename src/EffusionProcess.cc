// ------------------------------------------------------------
// Brahim Apr 10th 2003 : Effusion Process source file for ganilN01
//
// This file is based on lN01EffusionProcess.cc of li8N01
//-------------------------------------------------------------
////////////////////////////////////////////////////////////////////////
// Particle Boundary Process Class Implementation
////////////////////////////////////////////////////////////////////////
// File:        EffusionProcess.cc
// Description: Discrete Process -- Adsorption/Desorption of Particles
////////////////////////////////////////////////////////////////////////

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EffusionProcess::EffusionProcess(const G4String& processName)
: G4VDiscreteProcess(processName){
    kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EffusionProcess::~EffusionProcess(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VParticleChange*
EffusionProcess::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
    aParticleChange.Initialize(aTrack);

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
    
    G4bool bComputeMomentum = false;
    G4double theRediffusionProbability = 0.001;
    
    G4double theAdsorptionTime;
    if(aMaterialPost->GetName() == "Tantalum"){
        theAdsorptionTime = +100*CLHEP::ns;
        bComputeMomentum = true;
        theRediffusionProbability = 0.0;
    }
    
    if(aMaterialPost->GetName() == "Graphite"){
        theAdsorptionTime = +4420*CLHEP::ns;
        bComputeMomentum = true;
        theRediffusionProbability = 0.001;
    }

    if(aMaterialPost->GetName() == "Target"){
        theAdsorptionTime = +5000*CLHEP::ns;
        bComputeMomentum = true;
        theRediffusionProbability = 0.001;
    }
    
    ///////////////////////////////////////////////////////////////////////////////////
    
    // Check for redifussion probability
    if(G4UniformRand() > theRediffusionProbability){
        bComputeMomentum = true;
    }
    else{
        bComputeMomentum = false;
    }
    
    // Register Changes
    if(bComputeMomentum == true){
        aParticleChange.Initialize(aTrack);
        // Get the Global Point
        G4ThreeVector theGlobalPoint = pPostStepPoint->GetPosition();
        
        // Get Transport Navigator
        G4Navigator* theNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
        
        // Get the Local Point
        G4ThreeVector theLocalPoint = theNavigator->GetGlobalToLocalTransform().TransformPoint(theGlobalPoint);
        
        // Get the Local Normal : Normal points back into volume
        G4bool valid;
        G4ThreeVector theLocalNormal = theNavigator->GetLocalExitNormal(&valid);
        
        // Point the Local Normal back into volume
        if (valid) {
            theLocalNormal = -theLocalNormal;
        }
        // Transform the Local to the Global Normal
        G4ThreeVector theGlobalNormal = theNavigator->GetLocalToGlobalTransform().TransformAxis(theLocalNormal);
        
        aParticleChange.ProposeMomentumDirection(RDCosine(theGlobalNormal));
        aParticleChange.ProposeGlobalTime(aTrack.GetGlobalTime() + theAdsorptionTime) ;
        
        return &aParticleChange;
    }
    
    
    return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double EffusionProcess::GetMeanFreePath(const G4Track& ,
                                          G4double ,
                                          G4ForceCondition* condition)
{
    *condition = Forced;
    return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Random Direction according to the Cosine distribution : ncos9
G4ThreeVector EffusionProcess::RDCosine(G4ThreeVector N) const{
    N = N.unit();
    
    //
    //   Define the New Frame (P,Y1,N) where N is the new Z-axis.
    //  P and Y1 are perpendicular to N : parallel to the surface. And they
    //  will be determined using N angles in the Global frame : thn and phn
    //
    G4double phn = N.phi();
    G4double sphn = sin(phn);
    G4double cphn = cos(phn);
    
    // Define Y1 using phn: sphn and cphn
    G4ThreeVector Y1(-sphn,cphn,0.);
    
    // Define P using the cross product of Y1 and N
    G4ThreeVector P = Y1.cross(N);
    
    // Generate the polar angle : theta = N^RD, varies only between 0 and pi/2
    G4double rth = G4UniformRand();
    G4double sth = sqrt(rth);
    G4double cth = sqrt(1.- rth);
    
    // Generate the azimuthal angle : phi = P^RD varies between 0 and 2pi
    G4double rph = G4UniformRand();
    G4double ph = 4*acos(0.)*rph;
    
    // Define the random direction RD directly
    G4ThreeVector RD = sth*cos(ph)*P + sth*sin(ph)*Y1 + cth*N;
    
    
    return RD.unit();
}





