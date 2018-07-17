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

#include "RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "Run.hh"

#include "Analysis.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(): G4UserRunAction(){
    G4RunManager::GetRunManager()->SetPrintProgress(100);
    
    auto analysisManager = G4AnalysisManager::Instance();
    G4cout << "Using " << analysisManager->GetType() << G4endl;
    //analysisManager->SetNtupleMerging(true);

    // Create directories
    analysisManager->SetVerboseLevel(1);
    
    // Creating ntuple
    analysisManager->CreateNtuple("detector","Detector hits");
    analysisManager->CreateNtupleDColumn("t");
    analysisManager->CreateNtupleDColumn("A");
    analysisManager->CreateNtupleDColumn("Z");
    analysisManager->FinishNtuple();

    if(bSAVEALLPRIMARIES){
        analysisManager->CreateNtuple("ucx","UCx hits");
        analysisManager->CreateNtupleDColumn("t");
        analysisManager->CreateNtupleDColumn("A");
        analysisManager->CreateNtupleDColumn("Z");
        analysisManager->CreateNtupleDColumn("AP");
        analysisManager->CreateNtupleDColumn("AZ");
        analysisManager->CreateNtupleDColumn("x");
        analysisManager->CreateNtupleDColumn("y");
        analysisManager->CreateNtupleDColumn("z");
        analysisManager->FinishNtuple();
    }
    
    analysisManager->CreateH2("IT","Isotopes Table",120,-0.5,119.5,300,-0.5,299.5);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction(){
    delete G4AnalysisManager::Instance();
}

G4Run* RunAction::GenerateRun()
{ return new Run; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* /*run*/){
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->OpenFile("output");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* run){
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->Write();
    analysisManager->CloseFile();
    
    const Run* run_spes = static_cast<const Run*>(run);

    if (IsMaster())
    {
        std::ofstream fFileOut;
        fFileOut.open("isotope_table.dat",std::ofstream::out | std::ofstream::app);
        for (auto it : run_spes->fIsotopes){
            fFileOut << it.first << " , " << it.second << std::endl;
        }
        fFileOut.close();
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
