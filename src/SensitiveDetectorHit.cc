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

#include "SensitiveDetectorHit.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"

#ifdef G4MULTITHREADED
G4ThreadLocal G4Allocator<SensitiveDetectorHit>*
    SensitiveDetectorHitAllocator = 0;
#else
G4Allocator<SensitiveDetectorHit>
    SensitiveDetectorHitAllocator;
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SensitiveDetectorHit::SensitiveDetectorHit(){
    fTime = -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SensitiveDetectorHit::~SensitiveDetectorHit(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SensitiveDetectorHit::SensitiveDetectorHit(const SensitiveDetectorHit &right): G4VHit(){
    fTime = right.fTime;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const SensitiveDetectorHit& SensitiveDetectorHit::operator=(const SensitiveDetectorHit &right){
    fTime = right.fTime;
    return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int SensitiveDetectorHit::operator==(const SensitiveDetectorHit &/*right*/) const{
    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SensitiveDetectorHit::Draw(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const std::map<G4String,G4AttDef>*
    SensitiveDetectorHit::GetAttDefs() const{
    G4bool isNew;
    std::map<G4String,G4AttDef>* store = G4AttDefStore::GetInstance("SensitiveDetectorHit",isNew);
    if (isNew) {
        G4String Time("Time");
        (*store)[Time] = G4AttDef(Time,"Time","Physics","G4BestUnit","G4double");

    }
    return store;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::vector<G4AttValue>*
    SensitiveDetectorHit::CreateAttValues() const{
    std::vector<G4AttValue>* values = new std::vector<G4AttValue>;
    values->push_back(G4AttValue("Time",G4BestUnit(fTime,"Time"),""));
    return values;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SensitiveDetectorHit::Print(){
    G4cout << fTime << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
