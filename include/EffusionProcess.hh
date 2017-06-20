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

#ifndef EffusionProcess_h
#define EffusionProcess_h 1

#include "globals.hh"
#include "templates.hh"
#include "geomdefs.hh"
#include "Randomize.hh"
#include "G4ParticleChange.hh"
#include "G4VParticleChange.hh"

#include "G4VDiscreteProcess.hh"
#include "EffusionMaterialData.hh"
#include "G4ExtendedMaterial.hh"

class EffusionProcess : public G4VDiscreteProcess{
    
public:
    EffusionProcess(const G4String& processName = "effusion");
    ~EffusionProcess();
    
    G4double GetMeanFreePath(const G4Track& ,
                             G4double ,
                             G4ForceCondition* condition);
    
    G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                    const G4Step&  aStep);
    
public:
    void SampleMaxwellBoltzmannKineticEnergy(const G4Track&);
    void SampleLambertianDirection(const G4Track& aTrack);
    
private:
    G4double kCarTolerance;
    
    // Variables for Maxwell-Boltzmann
    G4double theSampledKineticEnergy;
    G4double thePx;
    G4double thePy;
    G4double thePz;
    G4double sqrtkTm;

    // Variables for Lambertian
    G4ThreeVector theLocalNormal;
    G4ThreeVector theLocalPoint;
    G4ThreeVector theGlobalNormal;
    G4ThreeVector theGlobalPoint;
    G4bool validLocalNorm;
    
private:
    EffusionMaterialData* GetMatData(const G4Track& aTrack){
        G4LogicalVolume* aLV = aTrack.GetNextVolume()->GetLogicalVolume();
        if(aLV->IsExtended() == true){
            G4ExtendedMaterial* aEM = (G4ExtendedMaterial*) aLV->GetMaterial();
            return (EffusionMaterialData*) aEM->RetrieveExtension("effusion");
        }
        else{
            return nullptr;
        }
    }
    
private:
    G4int fEffusionID;
    EffusionTrackData* GetTrackData(const G4Track&);


};

#endif /* EffusionProcess_h */
