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
// $Id: IronFilterDTGeneratorAction.cc $
//
/// \file IronFilterDTGeneratorAction.cc
/// \brief Implementation of the IronFilterDTGeneratorAction class

#include "IronFilterDTGeneratorAction.hh"

#include "G4RunManager.hh"
//#include "G4Navigator.hh"
//#include "G4PhysicalVolumeStore.hh"
//#include "G4VPhysicalVolume.hh"
//#include "G4SolidStore.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4RandomDirection.hh"
//#include "G4Neutron.hh"
//#include "G4TransportationManager.hh"
//#include "G4Navigator.hh"
#//include "G4GenericIon.hh"
//#include "G4IonTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IronFilterDTGeneratorAction::IronFilterDTGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleSource(0)
{
  fParticleSource = new G4ParticleGun();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IronFilterDTGeneratorAction::~IronFilterDTGeneratorAction()
{
  delete fParticleSource;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void IronFilterDTGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //Parameters of Energy_vs_Angle fitted into a cubic polynomial
  G4double e1 =   4.852e-7;
  G4double e2 =  -0.0001289;
  G4double e3 =  0.0008434;
  G4double e4 = 14.67;

  //Parameters of the Differential_crossection_vs_Angle fitted into a cubic polynomial
  G4double d1 = 3.042e-08 ;
  G4double d2 = -8.194e-06;
  G4double d3 = 2.906e-05;
  G4double d4 = 1.0;

  G4double phi = G4UniformRand()*2*3.14159265358979323846*radian;
  G4double angle = DT_dist(d1,d2, d3, d4);
  G4double theta = angle*(3.14159265358979323846/180)*radian;
  G4ThreeVector neutronDirection;
  neutronDirection.setRhoPhiTheta(-1.0,phi,theta);

  // set particle parameters
  fParticleSource->SetParticleMomentumDirection(neutronDirection);
  fParticleSource->SetParticleEnergy(e1*angle*angle*angle+e2*angle*angle+e3*angle+e4);
  G4ParticleDefinition* particleDefinition
    = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
  fParticleSource->SetParticleDefinition(particleDefinition);

  // Set source position
  fParticleSource->SetParticlePosition(G4ThreeVector(0., 0., 45.0*cm));
  fParticleSource->GeneratePrimaryVertex(anEvent);

}
G4double IronFilterDTGeneratorAction::DT_dist(G4double w1, G4double w2, G4double w3, G4double w4)
{
  G4bool flag = FALSE;
  G4double x;
  while(flag == FALSE){
    G4double weight = G4UniformRand();
    x = G4UniformRand()*180.0;
    G4double actual_weight=w1*x*x+w2*x*x+w3*x+w4;
    if(weight<actual_weight){
       flag = TRUE;
    }
  }
  return x;
}
