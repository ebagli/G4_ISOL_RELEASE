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

#ifndef EffusionMaterialData_h
#define EffusionMaterialData_h 1

#include "G4ios.hh"
#include "globals.hh"

#include "G4VMaterialExtension.hh"
#include "EffusionTrackData.hh"

class EffusionMaterialData : public G4VMaterialExtension{
public:
    
    EffusionMaterialData(const G4String&);
    virtual ~EffusionMaterialData();
    
public:
    void Print() const {G4cout << "Effusion Material Data" << G4endl;};
    
public:
    void SetAdsorptionTime(G4double aDouble){
        theAdsorptionTime=aDouble;};
    G4double GetAdsorptionTime(){
        return theAdsorptionTime;};
    
    void SetDiffusionProbability(G4double aDouble){
        theDiffusionProbability=aDouble;};
    G4double GetDiffusionProbability(){
        return theDiffusionProbability;};

    void SetAdsorptionProbability(G4double aDouble){
        theAdsorptionProbability=aDouble;};
    G4double GetAdsorptionProbability(){
        return theAdsorptionProbability;};

    void SetFullAdsorptionProbability(G4double aDouble){
        theFullAdsorptionProbability=aDouble;};
    G4double GetFullAdsorptionProbability(){
        return theFullAdsorptionProbability;};

    void SetDiffusionCoefficient(G4double aDouble){
        theDiffusionCoefficient=aDouble;};
    G4double GetDiffusionCoefficient(){
        return theDiffusionCoefficient;};

    void SetPorousDiffusionCoefficient(G4double aDouble){
        thePorousDiffusionCoefficient=aDouble;};
    G4double GetPorousDiffusionCoefficient(){
        return thePorousDiffusionCoefficient;};
    
private:
    // theAdsorptionTime is the average sticking time on the surface
    G4double theAdsorptionTime;
    
    // theDiffusionProbability is the probability that the atom enters the material
    // range [0,1]
    G4double theDiffusionProbability;

    // theAdsorptionProbability is the probability that, if the atom has not entered
    // in the material, the atom sticks to the surface
    // range [0,1]
    G4double theAdsorptionProbability;
    
    // theFullAdsorptionProbability is the probability that, once the atom has sticked
    // to the surface, it remains ther
    // range [0,1] relative to the adsorption probability
    G4double theFullAdsorptionProbability;

    // theDiffusionCoefficient is the diffusion coefficient of the material in the molecule cm2/s
    G4double theDiffusionCoefficient;

    // thePorousDiffusionCoefficient is the diffusion coefficient of the material in the porous cm2/s
    G4double thePorousDiffusionCoefficient;
};

#endif
