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

#include "EventAction.hh"

#include "G4RunManager.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

#include "SensitiveDetectorHit.hh"
#include "TargetSensitiveDetectorHit.hh"

#include "Analysis.hh"

EventAction::EventAction()
{
    fSD_ID = -1;
    fUCx_ID = -1;
    fVerboseLevel = 0;
    bSaveAllUCx = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

EventAction::~EventAction(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void EventAction::BeginOfEventAction(const G4Event*){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void EventAction::EndOfEventAction(const G4Event* evt){
    
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    
    if(fSD_ID == -1) {
        G4String sdName;
        if(SDman->FindSensitiveDetector(sdName="telescope",0)){
            fSD_ID = SDman->GetCollectionID(sdName="telescope/collection");
        }
    }
    if(fUCx_ID == -1) {
        G4String sdName;
        if(SDman->FindSensitiveDetector(sdName="ucx",0)){
            fUCx_ID = SDman->GetCollectionID(sdName="ucx/collection");
        }
    }

    G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
    TargetSensitiveDetectorHitsCollection* fSD = 0;
    TargetSensitiveDetectorHitsCollection* fUCx = 0;

    if(HCE)
    {
        if(fSD_ID != -1){
            G4VHitsCollection* aHCSD = HCE->GetHC(fSD_ID);
            fSD = (TargetSensitiveDetectorHitsCollection*)(aHCSD);
        }
        if(fUCx_ID != -1){
            G4VHitsCollection* aHCUCx = HCE->GetHC(fUCx_ID);
            fUCx = (TargetSensitiveDetectorHitsCollection*)(aHCUCx);
        }
    }

    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

    if(fSD)
    {
        int n_hit_sd = fSD->entries();
        for(int i1=0;i1<n_hit_sd;i1++)
        {
            TargetSensitiveDetectorHit* aHit = (*fSD)[i1];
            analysisManager->FillNtupleDColumn(0,0, aHit->GetTime()/CLHEP::s);
            analysisManager->FillNtupleDColumn(0,1, aHit->GetA());
            analysisManager->FillNtupleDColumn(0,2, aHit->GetZ());
            analysisManager->AddNtupleRow(0);
        }
    }

    if(fUCx)
    {
        int n_hit_sd = fUCx->entries();
        for(int i1=0;i1<n_hit_sd;i1++)
        {
            TargetSensitiveDetectorHit* aHit = (*fUCx)[i1];
            if(bSaveAllUCx){
                analysisManager->FillNtupleDColumn(1,0, aHit->GetTime()/CLHEP::s);
                analysisManager->FillNtupleDColumn(1,1, aHit->GetA());
                analysisManager->FillNtupleDColumn(1,2, aHit->GetZ());
                analysisManager->FillNtupleDColumn(1,3, aHit->GetAP());
                analysisManager->FillNtupleDColumn(1,4, aHit->GetZP());
                analysisManager->FillNtupleDColumn(1,5, aHit->GetWorldPos().x());
                analysisManager->FillNtupleDColumn(1,6, aHit->GetWorldPos().y());
                analysisManager->FillNtupleDColumn(1,7, aHit->GetWorldPos().z());
                analysisManager->AddNtupleRow(1);
            }
            
            analysisManager->FillH2(0,aHit->GetZ(),aHit->GetA());
            
        }
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
