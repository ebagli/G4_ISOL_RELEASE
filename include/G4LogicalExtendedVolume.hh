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
//
// $Id: G4CrystalMaterial.hh 94016 2015-11-05 10:14:49Z gcosmo $
//

//---------------------------------------------------------------------------
//
// ClassName:   G4LogicalExtendedVolume
//
// Description: XXX
//
// Class description:
//
// XXX
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// 21-04-16, created by E.Bagli

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4LogicalExtendedVolume_HH
#define G4LogicalExtendedVolume_HH 1

#include "G4LogicalVolume.hh"
#include "G4ExtendedMaterial.hh"
#include <algorithm>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4LogicalExtendedVolume : public G4LogicalVolume
{
public:
    G4LogicalExtendedVolume(G4VSolid* pSolid,
                           G4ExtendedMaterial* pMaterial,
                           const G4String& name,
                           G4FieldManager* pFieldMgr=0,
                           G4VSensitiveDetector* pSDetector=0,
                           G4UserLimits* pULimits=0,
                           G4bool optimise=true);

    ~G4LogicalExtendedVolume();
    
public:
    virtual G4bool IsExtended() const {return true;};
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
