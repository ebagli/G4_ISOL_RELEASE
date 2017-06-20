// ------------------------------------------------------------
// Brahim Apr 10th 2003 : Physics List source file for ganilN01
//
// This code is based on lN01PhysicsList.cc of li8N01
//----------------------------------------------------

#include "EffusionPhysicsList.hh"

#include "G4GenericIon.hh"
#include "G4ProcessManager.hh"

#include "EffusionProcess.hh"
#include "DiffusionProcess.hh"

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EffusionPhysicsList::EffusionPhysicsList(G4int)
:G4VPhysicsConstructor("ef10effusion"){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EffusionPhysicsList::~EffusionPhysicsList(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EffusionPhysicsList::ConstructParticle(){
    G4BosonConstructor  pBosonConstructor;
    pBosonConstructor.ConstructParticle();

    G4LeptonConstructor pLeptonConstructor;
    pLeptonConstructor.ConstructParticle();

    G4MesonConstructor pMesonConstructor;
    pMesonConstructor.ConstructParticle();

    G4BaryonConstructor pBaryonConstructor;
    pBaryonConstructor.ConstructParticle();
    
    G4IonConstructor pIonConstructor;
    pIonConstructor.ConstructParticle();
    
    G4ShortLivedConstructor pShortLivedConstructor;
    pShortLivedConstructor.ConstructParticle();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EffusionPhysicsList::ConstructProcess(){
    
    EffusionProcess* effusion = new EffusionProcess();
    DiffusionProcess* diffusion = new DiffusionProcess();
    
    G4ParticleTable::G4PTblDicIterator* aParticleIterator =
    G4ParticleTable::GetParticleTable()->GetIterator();
    aParticleIterator->reset();
    
    while( (*aParticleIterator)() ){
        G4ParticleDefinition* particle = aParticleIterator->value();

        G4ProcessManager* pManager = particle->GetProcessManager();

        pManager->AddDiscreteProcess(effusion);
        pManager->AddDiscreteProcess(diffusion);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
