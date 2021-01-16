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
#include "TFile.h"

#include "TTree.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "THStack.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TMinuit.h"
#include "TColor.h"
#include "TLine.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TApplication.h"


#include "IronFilterDTGeneratorAction.hh"
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

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
//#include "G4GenericIon.hh"
//#include "G4IonTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IronFilterDTGeneratorAction::IronFilterDTGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleSource(0)
{
  fParticleSource = new G4ParticleGun();

  //H_Eout = new TH1F("H_Eout","Energy",7500,0,15000);

  //H_XYout = new TH2D("H_XYout","XY",4000,-200,200,4000,-200,200); //1500

  //H_Yout = new TH1F("H_Yout","Energy",4000,-200,200); //1500

  //H_PXPYPZout = new TH3D("H_PXPYPZout","PXPYPZ",100,-1,1,100,-1,1,100,-1,1); //1500

  //H_PYout = new TH1F("H_PYout","Energy",100,-1,1); //1500

  //H_PZout = new TH1F("H_PZout","Energy",100,-1,1); //1500
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IronFilterDTGeneratorAction::~IronFilterDTGeneratorAction()
{
  delete fParticleSource;
  //delete H_Eout;
  //delete H_Xout;
  //delete H_Yout;
  //delete H_PXout;
  //delete H_PYout;
  //delete H_PZout;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void IronFilterDTGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
//  TH1F *H_Eout = new TH1F("H_Eout","Energy",1500000,0,15); //1500

//  TH1F *H_Xout = new TH1F("H_Xout","Energy",4000,-200,200); //1500

//  TH1F *H_Yout = new TH1F("H_Yout","Energy",4000,-200,200); //1500

//  TH1F *H_PXout = new TH1F("H_PXout","Energy",100,-1,1); //1500

//  TH1F *H_PYout = new TH1F("H_PYout","Energy",100,-1,1); //1500

//  TH1F *H_PZout = new TH1F("H_PZout","Energy",100,-1,1); //1500

  G4int eventID = anEvent->GetEventID();
  //cout<<eventID<<endl;
  //TFile* f;
  //TTree *t1;
  //Int_t nentries;
  //Double_t testE,testx,testy,testz,testxmom,testymom,testzmom;
  //Int_t testParticleID;

  if(eventID ==0)
  {

  f = TFile::Open("/Volumes/Seagate Backup Plus Drive/Simulation IronFilter/Shielding_DT/GEANT4_build_files/Inputfiles/Geom11_Doseage_1pulse_input.root");
  t1 = (TTree*)f->Get("IronFilter");
  nentries = (Int_t)t1->GetEntries();

  t1->SetBranchAddress("ParticleE",&testE);
  t1->SetBranchAddress("Xpos",&testx);
  t1->SetBranchAddress("Ypos",&testy);
  t1->SetBranchAddress("Zpos",&testz);
  t1->SetBranchAddress("Xmom",&testxmom);
  t1->SetBranchAddress("Ymom",&testymom);
  t1->SetBranchAddress("Zmom",&testzmom);
  t1->SetBranchAddress("ParticleType",&testParticleID);

}

  G4double Energy,X,Y,Z,Px,Py,Pz,angle,ParticleID;
  G4String particlename="gamma";
  G4ThreeVector neutronDirection;


  //Energy= 0.024*MeV;
  //cout<<Energy<<endl;

  t1->GetEntry(eventID);


  Energy=testE*MeV;
  X=testx*mm;
  Y=testy*mm;
  Z=testz*mm;;
  Px=testxmom;
  Py=testymom;
  Pz=testzmom;
  ParticleID=testParticleID;
  //ParticleID=2112;


  if(ParticleID==2112){
    particlename="neutron";
  }

  if (eventID==1959536-1) {
    f->Close();
  }


  neutronDirection.set(Px,Py,Pz);//negative sign becuase the ROOT file of the filter was incorrectly aligned

  // set particle parameters
  fParticleSource->SetParticleMomentumDirection(neutronDirection);
  fParticleSource->SetParticleEnergy(Energy);
  //cout<<Energy<<endl;
  //G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
  G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle(particlename);
  //cout<<ParticleID<<"  "<<particlename<<endl;
  fParticleSource->SetParticleDefinition(particleDefinition);

  // Set source position
  //fParticleSource->SetParticlePosition(G4ThreeVector(X, Z, Y));//please be careful about the coordinate in this case its rotated
  fParticleSource->SetParticlePosition(G4ThreeVector(X, Y, Z));
  fParticleSource->GeneratePrimaryVertex(anEvent);



}
