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
// $Id: IronFilterRunAction.cc $
//
/// \file IronFilterRunAction.cc
/// \brief Implementation of the IronFilterRunAction class

#include "IronFilterRunAction.hh"
#include "IronFilterAnalysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IronFilterRunAction::IronFilterRunAction()
 : G4UserRunAction()
{ 
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);     

  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in IronFilterAnalysis.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Default Settings
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetFileName("IronFilter");

  analysisManager->SetVerboseLevel(1);

  // Creating ntuple
  //
  analysisManager->CreateNtuple("IronFilter", "Particles");
  //analysisManager->CreateNtupleIColumn("EventID");
  //analysisManager->CreateNtupleIColumn("TrackID");
  //analysisManager->CreateNtupleIColumn("StepID");
  //analysisManager->CreateNtupleIColumn("ParticleType");
  analysisManager->CreateNtupleDColumn("ParticleE");
  analysisManager->CreateNtupleDColumn("Xpos");
  analysisManager->CreateNtupleDColumn("Ypos");
  analysisManager->CreateNtupleDColumn("Zpos");
  analysisManager->CreateNtupleDColumn("Time");
  analysisManager->CreateNtupleDColumn("Xmom");
  analysisManager->CreateNtupleDColumn("Ymom");
  analysisManager->CreateNtupleDColumn("Zmom");
  analysisManager->CreateNtupleIColumn("TestVolume");
  analysisManager->CreateNtupleIColumn("ParticleType");
  analysisManager->CreateNtupleIColumn("EventID");
  analysisManager->CreateNtupleIColumn("TrackID");
  analysisManager->CreateNtupleIColumn("StepID");
  analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IronFilterRunAction::~IronFilterRunAction()
{
  delete G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void IronFilterRunAction::BeginOfRunAction(const G4Run* /*run*/)
{   
  // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  analysisManager->OpenFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void IronFilterRunAction::EndOfRunAction(const G4Run* /*run*/)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // save ntuple
  analysisManager->Write();
  analysisManager->CloseFile();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
