//----------------------------------------------------
// Brahim Apr 10th 2003 : Physics List header file for ganilN01
//
// This file is based on lN01PhysicsList.hh of li8N01
//----------------------------------------------------

#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class PhysicsList: public G4VUserPhysicsList
{
  public:
  
    PhysicsList();
    ~PhysicsList();

  protected:
  
    // Construct particle and physics process
    void ConstructParticle();
    void ConstructProcess();
    void SetCuts();
    
    // Construct Effusion process
    void ConstructEffusion();

};

#endif

