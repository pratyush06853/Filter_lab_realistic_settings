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
// $Id: IronFilterSteppingAction.cc 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file IronFilterSteppingAction.cc
/// \brief Implementation of the IronFilterSteppingAction class

#include "IronFilterSteppingAction.hh"
#include "IronFilterEventAction.hh"
#include "IronFilterDetectorConstruction.hh"
#include "IronFilterAnalysis.hh"

#include "G4Neutron.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IronFilterSteppingAction::IronFilterSteppingAction(
                      const IronFilterDetectorConstruction* detectorConstruction,
                      IronFilterEventAction* eventAction)
  : G4UserSteppingAction(),
    fDetConstruction(detectorConstruction),
    fEventAction(eventAction)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IronFilterSteppingAction::~IronFilterSteppingAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void IronFilterSteppingAction::UserSteppingAction(const G4Step* step)
{
// Collect energy and number of scatters step by step

  G4StepPoint* preStep = step->GetPreStepPoint();
  G4StepPoint* postStep = step->GetPostStepPoint();

  G4Track* track = step->GetTrack();

  // get volume of the current step
  const G4VPhysicalVolume* volume = postStep->GetTouchableHandle()->GetVolume();
  const G4ParticleDefinition* particle = track->GetDefinition();
  const G4String processName = postStep->GetProcessDefinedStep()->GetProcessName();
  G4double energy = preStep->GetKineticEnergy();

//if ( (volume == fDetConstruction->GetboratedwatersolidPV() && processName == "Transportation" ))
//{
        // Produce code for which direction particle left GetboratedwatersolidPV()
        G4int test_volumeID;
        G4int flag=0;

        if(volume == fDetConstruction->GetTestCENTERPOINTPV()){
          test_volumeID = 0;
          flag=1;
         }
        else if(volume == fDetConstruction->GetTestFRONTSIDEPV()){
           test_volumeID = 1;
           flag=1;
        }
        else if(volume == fDetConstruction->GetTestREARSIDEPV()){
           test_volumeID = 2;
           flag=1;
        }
        else if(volume == fDetConstruction->GetTestRIGHTSIDEPV()){
         test_volumeID = 3;
         flag=1;
        }
        else if(volume == fDetConstruction->GetTestLEFTSIDEPV()){
          test_volumeID = 4;
          flag=1;
        }
        else if(volume == fDetConstruction->GetTestBOTTOMPV()){
         test_volumeID = 5;
         flag=1;
        }
        else if(volume == fDetConstruction->GetTestTOPPV()){
          test_volumeID = 6;
          flag=1;
        }
        else {
          test_volumeID =8;
          flag=0;
        }

        //way to kill particle less than a predefined energy to speed up thre simulation.
        //G4int flag=1;

        //if(energy*1000 < 1 ){
        //    track->SetTrackStatus(fStopAndKill);
        //    flag = 0;
          //cout<<"Track is killed as the energy of the particle is than 1 KeV"<<endl;
        //  }
        if(flag){
               // get analysis manager
               G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
               G4int eventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

               G4int trackID = track->GetTrackID();
               G4int stepID = track->GetCurrentStepNumber();
               G4int particle_ID = track->GetParticleDefinition()->GetPDGEncoding();

               // TODO: turn this into a tree that gets filled 2112 neutrons
               if( (particle_ID == 2112 || particle_ID == 22) ){
               //if( particle_ID == 2112 && (test_volumeID == 2 || test_volumeID == 3)  ){
                   //analysisManager->FillNtupleIColumn(0, eventID);
                   //analysisManager->FillNtupleIColumn(0, trackID);
                   //analysisManager->FillNtupleIColumn(1, stepID);
                   //analysisManager->FillNtupleIColumn(0, particle_ID);
                   analysisManager->FillNtupleDColumn(0, energy);
                   analysisManager->FillNtupleDColumn(1, track->GetPosition().x());
                   analysisManager->FillNtupleDColumn(2, track->GetPosition().y());
                   analysisManager->FillNtupleDColumn(3, track->GetPosition().z());
                   analysisManager->FillNtupleDColumn(4, track->GetGlobalTime());
                   analysisManager->FillNtupleDColumn(5, track->GetMomentumDirection().x());
                   analysisManager->FillNtupleDColumn(6, track->GetMomentumDirection().y());
                   analysisManager->FillNtupleDColumn(7, track->GetMomentumDirection().z());
                   analysisManager->FillNtupleIColumn(8, test_volumeID);
                   analysisManager->FillNtupleIColumn(9, particle_ID);
                   analysisManager->FillNtupleIColumn(10, eventID);
                   analysisManager->FillNtupleIColumn(11, trackID);
                   analysisManager->FillNtupleIColumn(12, stepID);
                   analysisManager->AddNtupleRow();
               }

         }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
