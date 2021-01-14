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
/// \file IronFilterDetectorMessenger.cc
/// \brief Implementation of the IronFilterDetectorMessenger class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "IronFilterDetectorMessenger.hh"


#include "IronFilterDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IronFilterDetectorMessenger::IronFilterDetectorMessenger(IronFilterDetectorConstruction * Det)
:G4UImessenger(),IronFilterDetector(Det),
 IronFilterDir(0),
 DetDir(0),
 PolyHeightcmd(0),
 FilterSpacingcmd(0),
 MultiplierLeadRadiuscmd(0),
 ModeratorAluminumRadiuscmd(0),
 MultiplierLeadHeightRearcmd(0),
 FilterCellSpacingcmd(0),
 ModeratorTitaniumHeightcmd(0),
 ModeratorAluminumHeightcmd(0),
 MultiplierLeadHeightFrontcmd(0),
 ModeratorTitaniumRadiuscmd(0),
 TestXcmd(0),
 TestYcmd(0),
 TestZcmd(0)


{
  IronFilterDir = new G4UIdirectory("/IronFilter/");
  IronFilterDir->SetGuidance("UI commands specific to this example");

  DetDir = new G4UIdirectory("/IronFilter/det/");
  DetDir->SetGuidance("detector construction commands");

  PolyHeightcmd= new G4UIcmdWithADoubleAndUnit("/IronFilter/det/setPolyHeight",this);
  PolyHeightcmd->SetGuidance("Set extra amount of Poly after Sc rod");
  PolyHeightcmd->SetParameterName("Radius1",false);
  PolyHeightcmd->SetRange("Radius1>0.");
  PolyHeightcmd->SetUnitCategory("Length");
  PolyHeightcmd->SetDefaultUnit("cm");
  PolyHeightcmd->AvailableForStates(G4State_PreInit);


  FilterSpacingcmd= new G4UIcmdWithADoubleAndUnit("/IronFilter/det/setFilterSpacing",this);
  FilterSpacingcmd->SetGuidance("Set Spacing between Filter and Helium Cell");
  FilterSpacingcmd->SetParameterName("Radius1",false);
  FilterSpacingcmd->SetRange("Radius1>0.");
  FilterSpacingcmd->SetUnitCategory("Length");
  FilterSpacingcmd->SetDefaultUnit("cm");
  FilterSpacingcmd->AvailableForStates(G4State_PreInit);

  MultiplierLeadRadiuscmd= new G4UIcmdWithADoubleAndUnit("/IronFilter/det/setMultiplierLeadRadius",this);
  MultiplierLeadRadiuscmd->SetGuidance("Set height Lead cylinder that encloses the DT ");
  MultiplierLeadRadiuscmd->SetParameterName("Height1",false);
  MultiplierLeadRadiuscmd->SetRange("Height1>0.");
  MultiplierLeadRadiuscmd->SetUnitCategory("Length");
  MultiplierLeadRadiuscmd->SetDefaultUnit("cm");
  MultiplierLeadRadiuscmd->AvailableForStates(G4State_PreInit);

  ModeratorAluminumRadiuscmd= new G4UIcmdWithADoubleAndUnit("/IronFilter/det/setModeratorAluminumRadius",this);
  ModeratorAluminumRadiuscmd->SetGuidance("Set radius of the Aluminum Moderator:");
  ModeratorAluminumRadiuscmd->SetParameterName("Radius2",false);
  ModeratorAluminumRadiuscmd->SetRange("Radius2>0.");
  ModeratorAluminumRadiuscmd->SetUnitCategory("Length");
  ModeratorAluminumRadiuscmd->SetDefaultUnit("cm");
  ModeratorAluminumRadiuscmd->AvailableForStates(G4State_PreInit);

  MultiplierLeadHeightRearcmd= new G4UIcmdWithADoubleAndUnit("/IronFilter/det/setMultiplierLeadHeightRear",this);
  MultiplierLeadHeightRearcmd->SetGuidance("Set Height of rear Lead cylinder:");
  MultiplierLeadHeightRearcmd->SetParameterName("Radius3",false);
  MultiplierLeadHeightRearcmd->SetRange("Radius3>0.");
  MultiplierLeadHeightRearcmd->SetUnitCategory("Length");
  MultiplierLeadHeightRearcmd->SetDefaultUnit("cm");
  MultiplierLeadHeightRearcmd->AvailableForStates(G4State_PreInit);

  FilterCellSpacingcmd= new G4UIcmdWithADoubleAndUnit("/IronFilter/det/setFilterCellSpacing",this);
  FilterCellSpacingcmd->SetGuidance("Set Spacing between Filter and Cell");
  FilterCellSpacingcmd->SetParameterName("Height2",false);
  FilterCellSpacingcmd->SetRange("Height2>0.");
  FilterCellSpacingcmd->SetUnitCategory("Length");
  FilterCellSpacingcmd->SetDefaultUnit("cm");
  FilterCellSpacingcmd->AvailableForStates(G4State_PreInit);

  ModeratorTitaniumHeightcmd= new G4UIcmdWithADoubleAndUnit("/IronFilter/det/setModeratorTitaniumHeight",this);
  ModeratorTitaniumHeightcmd->SetGuidance("Set Height of Titanium Moderator");
  ModeratorTitaniumHeightcmd->SetParameterName("Height3",false);
  ModeratorTitaniumHeightcmd->SetRange("Height3>0.");
  ModeratorTitaniumHeightcmd->SetUnitCategory("Length");
  ModeratorTitaniumHeightcmd->SetDefaultUnit("cm");
  ModeratorTitaniumHeightcmd->AvailableForStates(G4State_PreInit);

  ModeratorAluminumHeightcmd= new G4UIcmdWithADoubleAndUnit("/IronFilter/det/setModeratorAluminumHeight",this);
  ModeratorAluminumHeightcmd->SetGuidance("Set Height of the Fluental Moderator");
  ModeratorAluminumHeightcmd->SetParameterName("Height4",false);
  ModeratorAluminumHeightcmd->SetRange("Height4>0.");
  ModeratorAluminumHeightcmd->SetUnitCategory("Length");
  ModeratorAluminumHeightcmd->SetDefaultUnit("cm");
  ModeratorAluminumHeightcmd->AvailableForStates(G4State_PreInit);

  MultiplierLeadHeightFrontcmd= new G4UIcmdWithADoubleAndUnit("/IronFilter/det/setMultiplierLeadHeightFront",this);
  MultiplierLeadHeightFrontcmd->SetGuidance("Set Height of front Lead cylinder:");
  MultiplierLeadHeightFrontcmd->SetParameterName("Height5",false);
  MultiplierLeadHeightFrontcmd->SetRange("Height5>0.");
  MultiplierLeadHeightFrontcmd->SetUnitCategory("Length");
  MultiplierLeadHeightFrontcmd->SetDefaultUnit("cm");
  MultiplierLeadHeightFrontcmd->AvailableForStates(G4State_PreInit);

  ModeratorTitaniumRadiuscmd= new G4UIcmdWithADoubleAndUnit("/IronFilter/det/setModeratorTitaniumRadius",this);
  ModeratorTitaniumRadiuscmd->SetGuidance("Set Radius of the Titanium Moderator");
  ModeratorTitaniumRadiuscmd->SetParameterName("Radius4",false);
  ModeratorTitaniumRadiuscmd->SetRange("Radius4>0.");
  ModeratorTitaniumRadiuscmd->SetUnitCategory("Length");
  ModeratorTitaniumRadiuscmd->SetDefaultUnit("cm");
  ModeratorTitaniumRadiuscmd->AvailableForStates(G4State_PreInit);


  TestXcmd= new G4UIcmdWithADoubleAndUnit("/IronFilter/det/setTestX",this);
  TestXcmd->SetGuidance("Set X position of the Human");
  TestXcmd->SetParameterName("XPos",false);
  TestXcmd->SetRange("XPos>0.");
  TestXcmd->SetUnitCategory("Length");
  TestXcmd->SetDefaultUnit("cm");
  TestXcmd->AvailableForStates(G4State_PreInit);

  TestYcmd= new G4UIcmdWithADoubleAndUnit("/IronFilter/det/setTestY",this);
  TestYcmd->SetGuidance("Set Y position of the Human");
  TestYcmd->SetParameterName("YPos",false);
  TestYcmd->SetRange("YPos>0.");
  TestYcmd->SetUnitCategory("Length");
  TestYcmd->SetDefaultUnit("cm");
  TestYcmd->AvailableForStates(G4State_PreInit);

  TestZcmd= new G4UIcmdWithADoubleAndUnit("/IronFilter/det/setTestZ",this);
  TestZcmd->SetGuidance("Set Z position of the Human");
  TestZcmd->SetParameterName("ZPos",false);
  TestZcmd->SetRange("ZPos>0.");
  TestZcmd->SetUnitCategory("Length");
  TestZcmd->SetDefaultUnit("cm");
  TestZcmd->AvailableForStates(G4State_PreInit);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IronFilterDetectorMessenger::~IronFilterDetectorMessenger()
{
  delete MultiplierLeadRadiuscmd;
  delete PolyHeightcmd;
  delete FilterSpacingcmd;
  delete ModeratorAluminumRadiuscmd;
  delete MultiplierLeadHeightRearcmd;
  delete FilterCellSpacingcmd;
  delete ModeratorTitaniumHeightcmd;
  delete ModeratorAluminumHeightcmd;
  delete MultiplierLeadHeightFrontcmd;
  delete ModeratorTitaniumRadiuscmd;
  delete TestXcmd;
  delete TestYcmd;
  delete TestZcmd;
  delete DetDir;
  delete IronFilterDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void IronFilterDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if( command == PolyHeightcmd)
    { IronFilterDetector->SetPolyHeight(PolyHeightcmd->GetNewDoubleValue(newValue));}

  if( command == FilterSpacingcmd)
    { IronFilterDetector->SetFilterSpacing(FilterSpacingcmd->GetNewDoubleValue(newValue));}

  if( command == MultiplierLeadRadiuscmd)
    { IronFilterDetector->SetMultiplierLeadRadius(MultiplierLeadRadiuscmd->GetNewDoubleValue(newValue));}

  if( command == ModeratorAluminumRadiuscmd)
    { IronFilterDetector->SetModeratorAluminumRadius(ModeratorAluminumRadiuscmd->GetNewDoubleValue(newValue));}

  if( command == MultiplierLeadHeightRearcmd)
    { IronFilterDetector->SetMultiplierLeadHeightRear(MultiplierLeadHeightRearcmd->GetNewDoubleValue(newValue));}

  if( command == FilterCellSpacingcmd)
    { IronFilterDetector->SetFilterCellSpacing(FilterCellSpacingcmd->GetNewDoubleValue(newValue));}

  if( command == ModeratorTitaniumHeightcmd)
    { IronFilterDetector->SetModeratorTitaniumHeight(ModeratorTitaniumHeightcmd->GetNewDoubleValue(newValue));}

  if( command == ModeratorAluminumHeightcmd)
    { IronFilterDetector->SetModeratorAluminumHeight(ModeratorAluminumHeightcmd->GetNewDoubleValue(newValue));}

  if( command == MultiplierLeadHeightFrontcmd)
    { IronFilterDetector->SetMultiplierLeadHeightFront(MultiplierLeadHeightFrontcmd->GetNewDoubleValue(newValue));}

  if( command == ModeratorTitaniumRadiuscmd)
    { IronFilterDetector->SetModeratorTitaniumRadius(ModeratorTitaniumRadiuscmd->GetNewDoubleValue(newValue));}

  if( command == TestXcmd)
    { IronFilterDetector->SetTestX(TestXcmd->GetNewDoubleValue(newValue));}

  if( command == TestYcmd)
    { IronFilterDetector->SetTestY(TestYcmd->GetNewDoubleValue(newValue));}

  if( command == TestZcmd)
    { IronFilterDetector->SetTestZ(TestZcmd->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
