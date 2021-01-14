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
/// \file IronFilterDetectorMessenger.hh
/// \brief Definition of the IronFilterDetectorMessenger class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef IronFilterDetectorMessenger_h
#define IronFilterDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class IronFilterDetectorConstruction;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class IronFilterDetectorMessenger: public G4UImessenger
{
  public:
    IronFilterDetectorMessenger(IronFilterDetectorConstruction* );
   ~IronFilterDetectorMessenger();

    virtual void SetNewValue(G4UIcommand*, G4String);

  private:
    IronFilterDetectorConstruction*      IronFilterDetector;

    G4UIdirectory*             IronFilterDir;
    G4UIdirectory*             DetDir;

    G4UIcmdWithADoubleAndUnit* PolyHeightcmd;
    G4UIcmdWithADoubleAndUnit* FilterSpacingcmd;
    G4UIcmdWithADoubleAndUnit* MultiplierLeadRadiuscmd;
    G4UIcmdWithADoubleAndUnit* ModeratorAluminumRadiuscmd;
    G4UIcmdWithADoubleAndUnit* MultiplierLeadHeightRearcmd;
    G4UIcmdWithADoubleAndUnit* FilterCellSpacingcmd;
    G4UIcmdWithADoubleAndUnit* ModeratorTitaniumHeightcmd;
    G4UIcmdWithADoubleAndUnit* ModeratorAluminumHeightcmd;
    G4UIcmdWithADoubleAndUnit* MultiplierLeadHeightFrontcmd;
    G4UIcmdWithADoubleAndUnit* ModeratorTitaniumRadiuscmd;
    G4UIcmdWithADoubleAndUnit* TestXcmd;
    G4UIcmdWithADoubleAndUnit* TestYcmd;
    G4UIcmdWithADoubleAndUnit* TestZcmd;


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
