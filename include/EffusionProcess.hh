// ------------------------------------------------------------
// Brahim Apr 10th 2003 : Effusion process header file for ganilN01
//
// This file is based on lN01EffusionProcess.hh of li8N01
// ------------------------------------------------------------
//
////////////////////////////////////////////////////////////////////////
// Effusion Process Class Definition
////////////////////////////////////////////////////////////////////////
// File:        gnlN06EffusionProcess.hh
// Description: Discrete Process -- adsorption/desorption
////////////////////////////////////////////////////////////////////////

#ifndef EffusionProcess_h
#define EffusionProcess_h 1

/////////////
// Includes
/////////////

#include "globals.hh"
#include "templates.hh"
#include "geomdefs.hh"
#include "Randomize.hh"
#include "G4ParticleChange.hh"
#include "G4VParticleChange.hh"

#include "G4VDiscreteProcess.hh"

class EffusionProcess : public G4VDiscreteProcess{
    
public:
    EffusionProcess(const G4String& processName = "Effusion");
    ~EffusionProcess();
    
    G4double GetMeanFreePath(const G4Track& ,
                             G4double ,
                             G4ForceCondition* condition);
    
    G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                    const G4Step&  aStep);
    
private:
    G4double kCarTolerance;
    
    
protected:
    // Random Direction according to the Cosine law
    G4ThreeVector RDCosine(G4ThreeVector) const;
};

#endif /* EffusionProcess_h */
