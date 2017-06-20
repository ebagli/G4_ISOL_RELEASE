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
/// \file analysis//src/TargetSensitiveDetectorHit.cc
/// \brief Implementation of the TargetSensitiveDetectorHit class
//
// $Id$
// --------------------------------------------------------------
//
#include "TargetSensitiveDetectorHit.hh"
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
G4ThreadLocal G4Allocator<TargetSensitiveDetectorHit>* TargetSensitiveDetectorHitAllocator = 0;
#else
G4Allocator<TargetSensitiveDetectorHit> TargetSensitiveDetectorHitAllocator;
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TargetSensitiveDetectorHit::TargetSensitiveDetectorHit()
{
    fTime = -1.;
    fWorldPos = G4ThreeVector(0.,0.,0.);
    fLocalPos = G4ThreeVector(0.,0.,0.);
    fEnergy = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TargetSensitiveDetectorHit::~TargetSensitiveDetectorHit()
{
    fTime = -1.;
    fWorldPos = G4ThreeVector(0.,0.,0.);
    fLocalPos = G4ThreeVector(0.,0.,0.);
    fEnergy = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TargetSensitiveDetectorHit::TargetSensitiveDetectorHit(const TargetSensitiveDetectorHit &right): G4VHit(){
    fTrackID = right.fTrackID;
    fAP = right.fAP;
    fZP = right.fZP;
    fA = right.fA;
    fZ = right.fZ;
    fWorldPos = right.fWorldPos;
    fLocalPos = right.fLocalPos;
    fTime = right.fTime;
    fEnergy = right.fEnergy;
    fEnergyPrevious = right.fEnergyPrevious;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const TargetSensitiveDetectorHit& TargetSensitiveDetectorHit::operator=(const TargetSensitiveDetectorHit &right)
{
    fTrackID = right.fTrackID;
    fAP = right.fAP;
    fZP = right.fZP;
    fA = right.fA;
    fZ = right.fZ;
    fWorldPos = right.fWorldPos;
    fLocalPos = right.fLocalPos;
    fTime = right.fTime;
    fEnergy = right.fEnergy;
    fEnergyPrevious = right.fEnergyPrevious;
    return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int TargetSensitiveDetectorHit::operator==(const TargetSensitiveDetectorHit &/*right*/) const
{
    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TargetSensitiveDetectorHit::Draw()
{
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    if(pVVisManager)
    {
        G4Circle circle(fWorldPos);
        circle.SetScreenSize(2);
        circle.SetFillStyle(G4Circle::filled);
        G4Colour colour(1.,1.,0.);
        G4VisAttributes attribs(colour);
        circle.SetVisAttributes(attribs);
        pVVisManager->Draw(circle);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const std::map<G4String,G4AttDef>* TargetSensitiveDetectorHit::GetAttDefs() const
{
    G4bool isNew;
    std::map<G4String,G4AttDef>* store = G4AttDefStore::GetInstance("TargetSensitiveDetectorHit",isNew);
    if (isNew) {
        G4String HitType("HitType");
        (*store)[HitType] = G4AttDef(HitType,"Hit Type","Physics","","G4String");
        
        G4String ID("ID");
        (*store)[ID] = G4AttDef(ID,"ID","Physics","","G4int");

        G4String AP("AP");
        (*store)[AP] = G4AttDef(AP,"AP","Physics","","G4int");
        
        G4String ZP("ZP");
        (*store)[ZP] = G4AttDef(ZP,"ZP","Physics","","G4int");

        G4String A("A");
        (*store)[A] = G4AttDef(A,"A","Physics","","G4int");

        G4String Z("Z");
        (*store)[Z] = G4AttDef(Z,"Z","Physics","","G4int");

        G4String Time("Time");
        (*store)[Time] = G4AttDef(Time,"Time","Physics","G4BestUnit","G4double");
        
        G4String Pos("Pos");
        (*store)[Pos] = G4AttDef(Pos, "Position","Physics","G4BestUnit","G4ThreeVector");

        G4String En("En");
        (*store)[En] = G4AttDef(En, "Energy","Physics","G4BestUnit","G4double");

        G4String EnP("EnP");
        (*store)[EnP] = G4AttDef(EnP, "EnergyPrevious","Physics","G4BestUnit","G4double");
    }
    return store;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::vector<G4AttValue>* TargetSensitiveDetectorHit::CreateAttValues() const
{
    std::vector<G4AttValue>* values = new std::vector<G4AttValue>;
    
    values->push_back(G4AttValue("HitType","TargetSensitiveDetectorHit",""));
    
    values->push_back(G4AttValue("ID",G4UIcommand::ConvertToString(fTrackID),""));

    values->push_back(G4AttValue("AP",G4UIcommand::ConvertToString(fAP),""));
    
    values->push_back(G4AttValue("ZP",G4UIcommand::ConvertToString(fZP),""));

    values->push_back(G4AttValue("A",G4UIcommand::ConvertToString(fA),""));

    values->push_back(G4AttValue("Z",G4UIcommand::ConvertToString(fZ),""));

    values->push_back(G4AttValue("Time",G4BestUnit(fTime,"Time"),""));
    
    values->push_back(G4AttValue("Pos",G4BestUnit(fWorldPos,"Length"),""));
  
    values->push_back(G4AttValue("En",G4BestUnit(fEnergy,"Energy"),""));

    values->push_back(G4AttValue("EnP",G4BestUnit(fEnergyPrevious,"Energy"),""));

    return values;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TargetSensitiveDetectorHit::Print()
{
    G4cout << "  Layer[" << fTrackID << "] : time " << fTime/ns
    << " (nsec) --- local (x,y) " << fLocalPos.x()
    << ", " << fLocalPos.y() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
