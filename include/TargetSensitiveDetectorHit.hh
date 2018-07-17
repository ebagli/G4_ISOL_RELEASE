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
/// \file analysis//include/TargetSensitiveDetectorHit.hh
/// \brief Definition of the TargetSensitiveDetectorHit class
//
// $Id$
// --------------------------------------------------------------
//
#ifndef TargetSensitiveDetectorHit_h
#define TargetSensitiveDetectorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class G4AttDef;
class G4AttValue;

class TargetSensitiveDetectorHit : public G4VHit
{
public:
    TargetSensitiveDetectorHit();
    TargetSensitiveDetectorHit(G4int);
    
    virtual ~TargetSensitiveDetectorHit();
    
    TargetSensitiveDetectorHit(const TargetSensitiveDetectorHit &right);
    const TargetSensitiveDetectorHit& operator=(const TargetSensitiveDetectorHit &right);
    
    int operator==(const TargetSensitiveDetectorHit &right) const;
    
    inline void *operator new(size_t);
    inline void operator delete(void *aHit);
    
    inline float x();
    inline float y();
    
    virtual void Draw();
    virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
    virtual std::vector<G4AttValue>* CreateAttValues() const;
    virtual void Print();
    
private:
    G4int fTrackID;
    G4int fAP;
    G4int fZP;
    G4int fA;
    G4int fZ;
    G4double fTime;
    G4ThreeVector fLocalPos;
    G4ThreeVector fWorldPos;
    G4double fEnergy;
    G4double fEnergyPrevious;
    G4int fDisk;

public:
    inline void SetTrackID(G4int z) { fTrackID = z; }
    inline G4int GetTrackID() const { return fTrackID; }
    inline void SetAP(G4int z) { fAP = z; }
    inline G4int GetAP() const { return fAP; }
    inline void SetZP(G4int z) { fZP = z; }
    inline G4int GetZP() const { return fZP; }
    inline void SetA(G4int z) { fA = z; }
    inline G4int GetA() const { return fA; }
    inline void SetZ(G4int z) { fZ = z; }
    inline G4int GetZ() const { return fZ; }
    inline void SetTime(G4double t) { fTime = t; }
    inline G4double GetTime() const { return fTime; }
    inline void SetLocalPos(G4ThreeVector xyz) { fLocalPos = xyz; }
    inline G4ThreeVector GetLocalPos() const { return fLocalPos; }
    inline void SetWorldPos(G4ThreeVector xyz) { fWorldPos = xyz; }
    inline G4ThreeVector GetWorldPos() const { return fWorldPos; }
    inline void SetEnergy(G4double energy) { fEnergy = energy; }
    inline G4double GetEnergy() const { return fEnergy; }
    inline void SetEnergyPrevious(G4double energy) { fEnergyPrevious = energy; }
    inline G4double GetEnergyPrevious() const { return fEnergyPrevious; }
    inline void SetDiskNumber(G4int z) { fDisk = z; }
    inline G4int GetDiskNumber() const { return fDisk; }
};

typedef G4THitsCollection<TargetSensitiveDetectorHit> TargetSensitiveDetectorHitsCollection;

#ifdef G4MULTITHREADED
extern G4ThreadLocal G4Allocator<TargetSensitiveDetectorHit>* TargetSensitiveDetectorHitAllocator;
#else
extern G4Allocator<TargetSensitiveDetectorHit> TargetSensitiveDetectorHitAllocator;
#endif

inline void* TargetSensitiveDetectorHit::operator new(size_t)
{
#ifdef G4MULTITHREADED
    if(!TargetSensitiveDetectorHitAllocator) TargetSensitiveDetectorHitAllocator = new G4Allocator<TargetSensitiveDetectorHit>;
    return (void *) TargetSensitiveDetectorHitAllocator->MallocSingle();
#else
    void* aHit;
    aHit = (void*)TargetSensitiveDetectorHitAllocator.MallocSingle();
    return aHit;
#endif
}

inline void TargetSensitiveDetectorHit::operator delete(void* aHit)
{
#ifdef G4MULTITHREADED
    TargetSensitiveDetectorHitAllocator->FreeSingle((TargetSensitiveDetectorHit*) aHit);
#else
   TargetSensitiveDetectorHitAllocator.FreeSingle((TargetSensitiveDetectorHit*) aHit);
#endif
}

#endif


