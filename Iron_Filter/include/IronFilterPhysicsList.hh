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
// $Id: IronFilterPhysicsList.hh $
//
/// \file IronFilterPhysicsList.hh
/// \brief Definition of the IronFilterPhysicsList class

#ifndef IronFilterPhysicsList_h
#define IronFilterPhysicsList_h 1


#include "Shielding.hh"
//#include "QGSP_BIC_HP.hh"
#include "globals.hh"

class IronFilterPhysicsList : public Shielding
//class IronFilterPhysicsList : public QGSP_BIC_HP
{
  public:
    IronFilterPhysicsList(G4int verbose=1) : Shielding(verbose) {};
    //IronFilterPhysicsList(G4int verbose=1) : QGSP_BIC_HP(verbose) {};
    virtual ~IronFilterPhysicsList() {;}
    virtual void ConstructProcess();
};

#endif
