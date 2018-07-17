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
// $Id: Run.cc 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file Run.cc
/// \brief Implementation of the Run class

#include "Run.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4THitsMap.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "TargetSensitiveDetectorHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run()
 : G4Run(),
fUCx_ID(-1),
fIsotopes(0.)
{ }


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::RecordEvent(const G4Event* event)
{
  //Hit Detection System
    G4SDManager* SDman = G4SDManager::GetSDMpointer();

    if(fUCx_ID == -1) {
        G4String sdName;
        if(SDman->FindSensitiveDetector(sdName="ucx",0)){
            fUCx_ID = SDman->GetCollectionID(sdName="ucx/collection");
        }
    }

    G4HCofThisEvent * HCE = event->GetHCofThisEvent();
    TargetSensitiveDetectorHitsCollection* fUCx = 0;

    if(HCE)
    {
        if(fUCx_ID != -1){
            G4VHitsCollection* aHCUCx = HCE->GetHC(fUCx_ID);
            fUCx = (TargetSensitiveDetectorHitsCollection*)(aHCUCx);
        }
    }
    
    
    if(fUCx)
    {
        int n_hit_sd = fUCx->entries();
        for(int i1=0;i1<n_hit_sd;i1++)
        {
            TargetSensitiveDetectorHit* aHit = (*fUCx)[i1];
            fIsotopes[GetCode(aHit->GetA(),aHit->GetZ(),aHit->GetDiskNumber())] += 1;
//            auto search = fIsotopes.find(GetCode(aHit->GetA(),aHit->GetZ(),aHit->GetDiskNumber()));
//
//            if(search != fIsotopes.end()) {
//                search->second += 1;
//            }
//            else {
//                fIsotopes.insert(GetCode(aHit->GetA(),aHit->GetZ(),aHit->GetDiskNumber()),1)
//            }

        }
    }

   
  G4Run::RecordEvent(event);      
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* aRun)
{
  const Run* localRun = static_cast<const Run*>(aRun);
    
    for (auto it : localRun->fIsotopes){
        fIsotopes[it.first] += it.second;
    }

  G4Run::Merge(aRun); 
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
