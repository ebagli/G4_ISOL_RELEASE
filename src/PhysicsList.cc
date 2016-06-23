// ------------------------------------------------------------
// Brahim Apr 10th 2003 : Physics List source file for ganilN01
//
// This code is based on lN01PhysicsList.cc of li8N01
//----------------------------------------------------

// special header
#include "PhysicsList.hh"

// general header
#include "G4ParticleTypes.hh"

// Include other needed files
#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"
#include "G4IonConstructor.hh"
#include "G4GenericIon.hh"


// Include the Effusion Process
#include "EffusionProcess.hh"
#include "DiffusionProcess.hh"

// constructor
PhysicsList::PhysicsList()
{;}

// destructor
PhysicsList::~PhysicsList()
{;}

// member fuction 1 : particle definitions
void PhysicsList::ConstructParticle()
{
    // In this method, static member functions should be called
    // for all particles which you want to use.
    // This ensures that objects of these particle types will be
    // created in the program.
    
    // Define the alpha particle
    // Define the electron and the neutron (needed in G4IonTable).
    // Define the proton (needed to calculate ions cut values).
    
    //  nuclei
    G4Alpha::AlphaDefinition();
    
    // the electron
    G4Electron::ElectronDefinition();
    
    // the proton
    G4Proton::ProtonDefinition();
    
    // the neutron
    G4Neutron::NeutronDefinition();
    
    
    G4GenericIon::GenericIonDefinition();
}

// member funtion 2 : process definition
void PhysicsList::ConstructProcess()
{
    // Define transportation process
    AddTransportation();
    
    // Add the Effusion Process
    ConstructEffusion();
    
}

// Create and register the Effusion Process

void PhysicsList::ConstructEffusion()
{
    // Get the process manager for alpha
    G4ParticleDefinition* particle = G4GenericIon::GenericIonDefinition();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    
    // Construct Effusion process for alpha
    EffusionProcess* theEffusionProcess = new EffusionProcess();
    DiffusionProcess* theDiffusionProcess = new DiffusionProcess();
    
    // Register the Effusion process to alpha's process manager
    pmanager->AddDiscreteProcess(theEffusionProcess);
    pmanager->AddDiscreteProcess(theDiffusionProcess);
    
}

/**** Set Cuts ****/

// member function 3 : cut settings
void PhysicsList::SetCuts()
{
    // uppress error messages even in case e/gamma/proton do not exist
    G4int temp = GetVerboseLevel();
    SetVerboseLevel(0);
    
    //  " G4VUserPhysicsList::SetCutsWithDefault" method sets
    //   the default cut value for all particle types
    SetCutsWithDefault();
    
    // Retrieve verbose level
    SetVerboseLevel(temp);  
}


