//----------------------------------------------------
// Brahim Apr 10th 2003 : Physics List header file for ganilN01
//
// This file is based on lN01PhysicsList.hh of li8N01
//----------------------------------------------------

#ifndef EffusionPhysicsList_h
#define EffusionPhysicsList_h 1

#include "G4VPhysicsConstructor.hh"

class EffusionPhysicsList: public G4VPhysicsConstructor
{
  public:
    EffusionPhysicsList(G4int verbose =1);
    ~EffusionPhysicsList();

  protected:
    void ConstructParticle();
    void ConstructProcess();
};

#endif

