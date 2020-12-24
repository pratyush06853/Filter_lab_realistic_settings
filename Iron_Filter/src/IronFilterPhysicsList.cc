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
// $Id: IronFilterPhysicsList.cc 88427 2015-02-19 08:19:38Z gcosmo $
//
/// \file IronFilterPhysicsList.cc
/// \brief Implementation of the IronFilterPhysicsList class

#include "IronFilterPhysicsList.hh"
#include "G4HadronicInteractionRegistry.hh"
#include "G4HadronicProcess.hh"
#include "G4HadronicProcessType.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4Neutron.hh"
#include "G4NeutronHPThermalScattering.hh"
#include "G4NeutronHPThermalScatteringData.hh"
#include <assert.h>

#include "G4HadronElasticProcess.hh"
#include "G4HadronInelasticProcess.hh"

#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPInelastic.hh"

#include "G4NeutronHPThermalScattering.hh"

void IronFilterPhysicsList::ConstructProcess() {


  G4VModularPhysicsList::ConstructProcess();
  //Shielding::ConstructProcess();
  // Get the elastic scattering process for neutrons


  G4ParticleDefinition* nPD = G4Neutron::Definition();
  G4ProcessVector* pvec = nPD->GetProcessManager()->GetProcessList();


  G4HadronicProcess* nElastic = 0;
  for (G4int i=0; i<pvec->size(); i++) {
    if ((*pvec)[i]->GetProcessSubType() != fHadronElastic) continue;  //fHadronInelastic
    nElastic = dynamic_cast<G4HadronicProcess*>((*pvec)[i]);
    break;
  }
  assert(nElastic != 0);
  // Get the "regular" HP elastic scattering model, exclude thermal scattering region
  G4HadronicInteraction* nElasticHP = G4HadronicInteractionRegistry::Instance()->FindModel("NeutronHPElastic");
  assert(nElasticHP != 0);
  nElasticHP->SetMinEnergy(4.*CLHEP::eV);
  // Attach HP thermal scattering model and data files to process
  nElastic->RegisterMe(new G4NeutronHPThermalScattering);
  nElastic->AddDataSet(new G4NeutronHPThermalScatteringData);

/*
  G4HadronicProcess* inElastic = 0;
  for (G4int i=0; i<pvec->size(); i++) {
    if ((*pvec)[i]->GetProcessSubType() != fHadronInelastic) continue;
    inElastic = dynamic_cast<G4HadronicProcess*>((*pvec)[i]);
    break;
  }
  assert(inElastic != 0);
  // Get the "regular" HP inelastic scattering model, exclude thermal scattering region
  G4HadronicInteraction* inElasticHP = G4HadronicInteractionRegistry::Instance()->FindModel("NeutronHPInelastic");
  assert(inElasticHP != 0);
  inElasticHP->SetMinEnergy(4.*CLHEP::eV);
  // Attach HP thermal scattering model and data files to process
  inElastic->RegisterMe(new G4NeutronHPThermalScattering);
  inElastic->AddDataSet(new G4NeutronHPThermalScatteringData);
*/



}
