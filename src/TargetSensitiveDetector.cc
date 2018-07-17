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
/// \file analysis//src/TargetSensitiveDetector.cc
/// \brief Implementation of the TargetSensitiveDetector class
//
// $Id$
// --------------------------------------------------------------
//
#include "TargetSensitiveDetector.hh"
#include "TargetSensitiveDetectorHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4Navigator.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"

#include "G4TrajectoryContainer.hh"
#include "G4RunManager.hh"

TargetSensitiveDetector::TargetSensitiveDetector(G4String name,G4int type):G4VSensitiveDetector(name),
detType(type)

{
    G4String HCname;
    collectionName.insert(HCname="collection");
    fAParent.clear();
    fZParent.clear();
    fHCID = -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TargetSensitiveDetector::~TargetSensitiveDetector(){
    ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TargetSensitiveDetector::Initialize(G4HCofThisEvent*HCE)
{
    fHitsCollection = new TargetSensitiveDetectorHitsCollection(SensitiveDetectorName,collectionName[0]);
    if(fHCID<0){
        fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
    }
    HCE->AddHitsCollection(fHCID,fHitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

bool TargetSensitiveDetector::ProcessHits(G4Step *aStep,G4TouchableHistory* /*ROhist*/)
{
    G4Track* vTrack = aStep->GetTrack();

    G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
    G4StepPoint* postStepPoint = aStep->GetPostStepPoint();

    std::map<int,int>::const_iterator it = fAParent.find(vTrack->GetTrackID());
    if(it==fAParent.end()){
        fAParent[vTrack->GetTrackID()] = vTrack->GetDefinition()->GetAtomicMass();
        fZParent[vTrack->GetTrackID()] = vTrack->GetDefinition()->GetAtomicNumber();
    }

    
    if(vTrack->GetCurrentStepNumber() != 1 && detType==0){
        return true;
    }
    if(!(postStepPoint->GetStepStatus() == fGeomBoundary) && detType==1){
        return true;
    }


    G4TouchableHistory* theTouchable = (G4TouchableHistory*)(preStepPoint->GetTouchable());

    G4ThreeVector worldPos = preStepPoint->GetPosition();
    G4ThreeVector localPos = theTouchable->GetHistory()->GetTopTransform().TransformPoint(worldPos);
    
    TargetSensitiveDetectorHit* aHit = new TargetSensitiveDetectorHit();
    aHit->SetTrackID(vTrack->GetTrackID());
    
    if(vTrack->GetTrackID() == 1 ){
        aHit->SetAP(-1);
        aHit->SetZP(-1);
    }
    else{
        aHit->SetAP(fAParent[vTrack->GetParentID()]);
        aHit->SetZP(fZParent[vTrack->GetParentID()]);
    }
    aHit->SetA(vTrack->GetDefinition()->GetAtomicMass());
    aHit->SetZ(vTrack->GetDefinition()->GetAtomicNumber());
    aHit->SetWorldPos(worldPos);
    aHit->SetLocalPos(localPos);
    aHit->SetTime(preStepPoint->GetGlobalTime());
    aHit->SetEnergy(preStepPoint->GetKineticEnergy());
    aHit->SetEnergyPrevious(fEnParent);

    G4VPhysicalVolume* thePhysical = theTouchable->GetVolume(0);
    G4int copyNo = thePhysical->GetCopyNo();
    aHit->SetDiskNumber(copyNo);

    fHitsCollection->insert(aHit);

    fEnParent = preStepPoint->GetKineticEnergy();
    return true;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TargetSensitiveDetector::EndOfEvent(G4HCofThisEvent* /*HCE*/)
{
   fAParent.clear();
    fZParent.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
