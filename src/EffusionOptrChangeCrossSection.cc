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
/// \file GB01/src/EffusionOptrChangeCrossSection.cc
/// \brief Implementation of the EffusionOptrChangeCrossSection class
//
#include "EffusionOptrChangeCrossSection.hh"
#include "G4BiasingProcessInterface.hh"
#include "G4BOptnChangeCrossSection.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4VProcess.hh"

#include "Randomize.hh"

#include "G4InteractionLawPhysical.hh"
#include "EffusionTrackData.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EffusionOptrChangeCrossSection::EffusionOptrChangeCrossSection(G4String particleName,
                                                         G4String         name)
: G4VBiasingOperator(name),
fSetup(true),
fEffusionID(-1)
{
    fParticleToBias = G4ParticleTable::GetParticleTable()->FindParticle(particleName);
    
    if ( fParticleToBias == 0 )
    {
        G4ExceptionDescription ed;
        ed << "Particle `" << particleName << "' not found !" << G4endl;
        G4Exception("EffusionOptrChangeCrossSection(...)",
                    "exGB01.01",
                    JustWarning,
                    ed);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EffusionOptrChangeCrossSection::~EffusionOptrChangeCrossSection()
{
    for ( std::map< const G4BiasingProcessInterface*, G4BOptnChangeCrossSection* >::iterator
         it = fChangeCrossSectionOperations.begin() ;
         it != fChangeCrossSectionOperations.end() ;
         it++ ) delete (*it).second;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EffusionOptrChangeCrossSection::StartRun()
{
    // --------------
    // -- Setup stage:
    // ---------------
    // -- Start by collecting processes under biasing, create needed biasing
    // -- operations and associate these operations to the processes:
    if ( fSetup )
    {
        const G4ProcessManager* processManager = fParticleToBias->GetProcessManager();
        const G4BiasingProcessSharedData* sharedData =
        G4BiasingProcessInterface::GetSharedData( processManager );
        if ( sharedData ) // -- sharedData tested, as is can happen a user attaches an operator to a
            // -- volume but without defined BiasingProcessInterface processes.
        {
            for ( size_t i = 0 ; i < (sharedData->GetPhysicsBiasingProcessInterfaces()).size(); i++ )
            {
                const G4BiasingProcessInterface* wrapperProcess =
                (sharedData->GetPhysicsBiasingProcessInterfaces())[i];
                G4String operationName = "EffusionXSchange-" +
                wrapperProcess->GetWrappedProcess()->GetProcessName();
                fChangeCrossSectionOperations[wrapperProcess] =
                new G4BOptnChangeCrossSection(operationName);
            }
        }
        fSetup = false;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VBiasingOperation*
EffusionOptrChangeCrossSection::ProposeOccurenceBiasingOperation(const G4Track*            track,
                                                              const G4BiasingProcessInterface*
                                                              callingProcess)
{
    
    G4double analogInteractionLength = callingProcess->GetWrappedProcess()->GetCurrentInteractionLength();
    if ( analogInteractionLength > DBL_MAX/10. ) return 0;

    if(fEffusionID==-1){
        fEffusionID = G4PhysicsModelCatalog::GetIndex("effusion");
    }
    
    EffusionTrackData* trackdata =
    (EffusionTrackData*)(track->GetAuxiliaryTrackInformation(fEffusionID));
    if(trackdata==nullptr) return 0;
    
    const G4DynamicParticle* aParticle = track->GetDynamicParticle();
    G4double stickedTime = trackdata->GetTimeSticked() ;
    G4double aMass = aParticle->GetMass();
    G4double pathlength = 0.0;

    if(stickedTime > 0.){
        G4double betaGamma = aParticle->GetTotalMomentum()/aMass;
        pathlength = CLHEP::c_light*stickedTime*betaGamma;
        if (pathlength < 0.0) {
            pathlength = 0.0;
        }
    }
    trackdata->SetTimeSticked(0.) ;

    G4double analogXS = 1./(analogInteractionLength);
    if (pathlength > analogInteractionLength){
        analogXS = DBL_MAX;
    }
    else{
        analogXS = 1./(analogInteractionLength - pathlength);
    }

    G4BOptnChangeCrossSection*   operation = fChangeCrossSectionOperations[callingProcess];
    G4VBiasingOperation* previousOperation = callingProcess->GetPreviousOccurenceBiasingOperation();
    
    if ( previousOperation == 0 )
    {
        operation->SetBiasedCrossSection( analogXS );
        operation->Sample();
    }
    else
    {
        if (  previousOperation != operation )
        {
            // -- should not happen !
            G4ExceptionDescription ed;
            ed << " Logic problem in operation handling !" << G4endl;
            G4Exception("GB01BOptrChangeCrossSection::ProposeOccurenceBiasingOperation(...)",
                        "exGB01.02",
                        JustWarning,
                        ed);
            return 0;
        }
        if ( operation->GetInteractionOccured() )
        {
            operation->SetBiasedCrossSection( analogXS );
            operation->Sample();
        }
        else
        {
            // -- update the 'interaction length' and underneath 'number of interaction lengths'
            // -- for past step  (this takes into accout the previous step cross-section value)
            operation->UpdateForStep( callingProcess->GetPreviousStepSize() );
            // -- update the cross-section value:
            operation->SetBiasedCrossSection( analogXS );
            // -- forces recomputation of the 'interaction length' taking into account above
            // -- new cross-section value [tricky : to be improved]
            operation->UpdateForStep( 0.0 );
        }
    }
return operation;
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EffusionOptrChangeCrossSection::
OperationApplied(const G4BiasingProcessInterface*           callingProcess,
                 G4BiasingAppliedCase,
                 G4VBiasingOperation*             occurenceOperationApplied,
                 G4double,
                 G4VBiasingOperation*,
                 const G4VParticleChange*                                  )
{
    G4BOptnChangeCrossSection* operation = fChangeCrossSectionOperations[callingProcess];
    if ( operation ==  occurenceOperationApplied ) operation->SetInteractionOccured();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
