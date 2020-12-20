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
// $Id: IronFilterPrimaryGeneratorAction.cc $
//
/// \file IronFilterPrimaryGeneratorAction.cc
/// \brief Implementation of the IronFilterPrimaryGeneratorAction class

#include "IronFilterPrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4Navigator.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4SolidStore.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4RandomDirection.hh"
#include "G4Neutron.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4GenericIon.hh"
#include "G4IonTable.hh"
#include "G4VisExtent.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IronFilterPrimaryGeneratorAction::IronFilterPrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleSource(0),
   gamma_energy(1691.0*keV),
   generator_mode(""),
   first_event(true),
   disk_global_position(0., 0., 0.),
   source_PV(0),
   boundingSphereRadius(0),
   process_threshold(-1666*keV),
   target_mass(8.394794 * GeV),
   particleDefinition(0),
   event_position(0., 0., 0.),
   accepted(false),
   gamma_direction(0., 0., 0.),
   neutron_direction(0., 0., 0.),
   neutron_angle(0),
   neutron_energy(0)

{
  fParticleSource = new G4ParticleGun();

  gNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();

  // pick whether to generate photoneutrons or Sb-124 decay gammas here
  generator_mode = "gammas";
  //generator_mode = "neutrons";

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IronFilterPrimaryGeneratorAction::~IronFilterPrimaryGeneratorAction()
{
  delete fParticleSource;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void IronFilterPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  // Do the setup stuff here because detector construction initialized after
  // primary generator action
  if(first_event){
    G4VPhysicalVolume* source_disk_PV = G4PhysicalVolumeStore::GetInstance()->GetVolume("source_disk");
    G4VPhysicalVolume* current_PV = source_disk_PV;

    // here we get the global coordinates of the disk by looping up through mother volumes
    while( current_PV->GetName() != "World" ){
      disk_global_position += current_PV->GetTranslation();
      // This is really ugly... only one name per PV and LV and they must match!!
      current_PV = G4PhysicalVolumeStore::GetInstance()->GetVolume(current_PV->GetMotherLogical()->GetName());
    }

    gNavigator->LocateGlobalPointAndSetup(disk_global_position);
    G4TouchableHistoryHandle gTouchable = gNavigator->CreateTouchableHistoryHandle();

    // This gets the parent volume of the source disk (BeO)
    source_PV = gTouchable->GetSolid(1);

    // The extent is a bounding rectangular solid of the source volume
    G4VisExtent source_extent = source_PV->GetExtent();
    G4double xmin = source_extent.GetXmin();
    G4double xmax = source_extent.GetXmax();
    G4double ymin = source_extent.GetYmin();
    G4double ymax = source_extent.GetYmax();
    G4double zmin = source_extent.GetZmin();
    G4double zmax = source_extent.GetZmax();

    // Use maximum ray in source as radius of bounding sphere to ensure total envelopment
    boundingSphereRadius = sqrt((xmax-xmin)*(xmax-xmin) + (ymax-ymin)*(ymax-ymin) + (zmax-zmin)*(zmax-zmin));

    first_event = false;
  }

  if(generator_mode == "gammas"){
    particleDefinition = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(51, 124, 0);

    fParticleSource->SetParticlePosition(disk_global_position);
    fParticleSource->SetParticleDefinition(particleDefinition);

    fParticleSource->SetParticleMomentumDirection( G4RandomDirection() );
  } else if (generator_mode == "neutrons" ){
    particleDefinition = G4Neutron::Definition();

    accepted = false;
    gamma_direction = G4ThreeVector(0., 0., 0.);

    while( !accepted ){
      accepted = GenerateNeutronPoint(boundingSphereRadius, disk_global_position,
                                      gNavigator, source_PV, gamma_direction, event_position);
    }

    neutron_direction = G4RandomDirection();
    neutron_angle = gamma_direction.angle(neutron_direction);

    neutron_energy = PhotoneutronEnergy(process_threshold, neutron_angle,
                                        gamma_energy, target_mass);

    fParticleSource->SetParticleDefinition(particleDefinition);

    fParticleSource->SetParticlePosition(event_position);
    fParticleSource->SetParticleEnergy(neutron_energy);
    fParticleSource->SetParticleMomentumDirection(neutron_direction);
  } else {
    std::cout << "ERROR! Bad generator mode definition: " << generator_mode << std::endl;
  }

  fParticleSource->GeneratePrimaryVertex(anEvent);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool IronFilterPrimaryGeneratorAction::GenerateNeutronPoint(G4double sphereRadius,
                G4ThreeVector sourcePosition, G4Navigator* navigator, G4VSolid* parentVolume,
                G4ThreeVector& gammaDirection, G4ThreeVector& interactionPoint)
{
  G4double interactionLength = G4UniformRand() * sphereRadius;

  gammaDirection = G4RandomDirection();

  // This is the baccarat way:
  /*G4ThreeVector BaccSource::GetRandomDirection()
{
    cosTheta = -1. + 2. * G4UniformRand();
    sinTheta = sqrt(1. - cosTheta * cosTheta);
    phi = 2. * 3.141592653589793238 * G4UniformRand();
    xDir = sinTheta * cos(phi);
    yDir = sinTheta * sin(phi);
    zDir = cosTheta;
    G4ThreeVector direction(xDir, yDir, zDir);
    return direction;
}
*/
  //gammaDirection = GetRandomDirection();


  interactionPoint = sourcePosition + gammaDirection * interactionLength;

  navigator->LocateGlobalPointAndSetup(interactionPoint);

  return navigator->CreateTouchableHistoryHandle()->GetSolid() == parentVolume;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double IronFilterPrimaryGeneratorAction::PhotoneutronEnergy(G4double threshold,
                G4double scatteringAngle, G4double gammaEnergy, G4double nucleusMass)
{
  G4double neutronMass = G4Neutron::Definition()->GetPDGMass();
  G4double firstTerm = nucleusMass * (gammaEnergy + threshold) / (neutronMass + nucleusMass);
  G4double secondTerm = gammaEnergy * sqrt((2 * neutronMass * nucleusMass) *
       (neutronMass + nucleusMass) * (gammaEnergy + threshold)) * cos(scatteringAngle) /
      ((neutronMass + nucleusMass) * (neutronMass + nucleusMass));

  return firstTerm + secondTerm;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
