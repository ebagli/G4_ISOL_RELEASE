// ------------------------------------------------------------
// Brahim Apr 10th 2003 : Diffusion process header file for ganilN01
//
// This file is based on lN01DiffusionProcess.hh of li8N01
// ------------------------------------------------------------
//
////////////////////////////////////////////////////////////////////////
// Diffusion Process Class Definition
////////////////////////////////////////////////////////////////////////
// File:        gnlN06DiffusionProcess.hh
// Description: Discrete Process -- adsorption/desorption
////////////////////////////////////////////////////////////////////////

#ifndef DiffusionProcess_h
#define DiffusionProcess_h 1

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

class DiffusionProcess : public G4VDiscreteProcess{
    
public:
    DiffusionProcess(const G4String& processName = "Diffusion");
    ~DiffusionProcess();
    
    G4double GetMeanFreePath(const G4Track& ,
                             G4double ,
                             G4ForceCondition* condition);
    
    G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                    const G4Step&  aStep);
    
private:
    G4double kCarTolerance;
};

#endif /* DiffusionProcess_h */
