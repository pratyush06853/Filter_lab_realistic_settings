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
// $Id: IronFilterDTGeneratorAction.hh $
//
/// \file IronFilterDTGeneratorAction.hh
/// \brief Definition of the IronFilterDTGeneratorAction class

#include "TFile.h"

#include "TTree.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH3D.h"
#include "THStack.h"
#include "TH2D.h"
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

#ifndef IronFilterDTGeneratorAction_h
#define IronFilterDTGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"

class G4ParticleGun;
class G4Event;

class IronFilterDTGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  IronFilterDTGeneratorAction();
  virtual ~IronFilterDTGeneratorAction();

  virtual void GeneratePrimaries(G4Event* event);

  // set methods
  void SetRandomFlag(G4bool value);

private:
  G4ParticleGun*  fParticleSource;
  TFile* f;
  TTree *t1;
  Int_t nentries;
  Double_t testE,testx,testy,testz,testxmom,testymom,testzmom;
  Int_t testParticleID;
  //TH1F *H_Eout;
  //G4double energy_beam[1000], X_beam[1000],Y_beam[1000],Px_beam[1000],Py_beam[1000],Pz_beam[1000];//857 events so 1000 works
  //G4double energy_beam[3500], X_beam[3500],Y_beam[3500],Px_beam[3500],Py_beam[3500],Pz_beam[3500];//3500 events so 3428works
  //TH2D *H_XYout;
  //TH1F *H_Yout;
  //TH3D *H_PXPYPZout;
  //TH1F *H_PYout;
  //TH1F *H_PZout;

  G4double DT_dist(G4double w1, G4double w2, G4double w3, G4double w4);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
