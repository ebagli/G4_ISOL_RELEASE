#include "G4MTRunManager.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"
#include "G4UIterminal.hh"

#include "DetectorConstruction.hh"
#include "UserActionInitialization.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif


#include "PhysicsList.hh"
#include "QGSP_BERT.hh"
#include "G4GenericBiasingPhysics.hh"

#include <iostream>
#include <fstream>
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main (int argc,char** argv) {
    
    G4MTRunManager* runManager = new G4MTRunManager;
    runManager->SetNumberOfThreads(G4Threading::G4GetNumberOfCores() - 2);
    
    // Set mandatory initialization classes
    runManager->SetUserInitialization(new DetectorConstruction());
    if(argc>2){
        if(atoi(argv[2])==1){
            runManager->SetUserInitialization(new QGSP_BERT());
        }
    }
    else{
        PhysicsList* physlist = new PhysicsList();
        G4GenericBiasingPhysics* biasingPhysics = new G4GenericBiasingPhysics();
        biasingPhysics->PhysicsBiasAllCharged();
        physlist->RegisterPhysics(biasingPhysics);
        runManager->SetUserInitialization(physlist);
    }
    
    runManager->SetUserInitialization(new UserActionInitialization());
    
    G4UImanager* UI = G4UImanager::GetUIpointer();
    
    if (argc!=1)   // batch mode
    {
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UI->ApplyCommand(command+fileName);
    }
    
    else           //define visualization and UI terminal for interactive mode
    {
#ifdef G4VIS_USE
        G4VisManager* visManager = new G4VisExecutive;
        visManager->Initialize();
#endif
        
#ifdef G4UI_USE
        G4UIExecutive * ui = new G4UIExecutive(argc,argv);
        ui->SessionStart();
        delete ui;
#endif
        
#ifdef G4VIS_USE
        delete visManager;
#endif
    }
    
    delete runManager;
    
    return 0;
    
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

