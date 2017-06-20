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

#include "StackingAction.hh"
#include "G4Track.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::StackingAction(){
    fKillSecondary  = 1;
    
    
    // -- Define messengers:
    fKillSecondaryMessenger =
    new G4GenericMessenger(this, "/stacking/","Biasing control" );
    
    G4GenericMessenger::Command& killSecondaryCmd =
    fKillSecondaryMessenger->DeclareProperty("killSecondary", fKillSecondary,
                                             "Kill secondary particles yes/no." );

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::~StackingAction(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track* aTrack){
    G4ClassificationOfNewTrack status = fUrgent;
    
    if(fKillSecondary == 1){
        if(aTrack->GetParticleDefinition()->GetParticleType()!="nucleus"){
            status = fKill;
            //G4cout << "DEAD  " << aTrack->GetTrackID() << " " << aTrack->GetParticleDefinition()->GetParticleType() << " " << aTrack->GetParticleDefinition()->GetParticleName() << G4endl;

        }
        else{
            //G4cout << "ALIVE " << aTrack->GetTrackID() << " " << aTrack->GetParticleDefinition()->GetParticleType() << " " << aTrack->GetParticleDefinition()->GetParticleName() << G4endl;
        }
    }
    else if(fKillSecondary == 2){
        if(aTrack->GetTrackID()>1){
            status = fKill;
        }
    }

    return status;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StackingAction::SetKillStatus(G4bool value){
    fKillSecondary = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
