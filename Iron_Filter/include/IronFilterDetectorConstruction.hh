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
// $Id: IronFilterDetectorConstruction.hh $
//
/// \file IronFilterDetectorConstruction.hh
/// \brief Definition of the IronFilterDetectorConstruction class

#ifndef IronFilterDetectorConstruction_h
#define IronFilterDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4RotationMatrix.hh"

class G4VPhysicalVolume;
//for messenger
class IronFilterDetectorMessenger;
//
class IronFilterDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    IronFilterDetectorConstruction();
    virtual ~IronFilterDetectorConstruction();

  public:
    virtual G4VPhysicalVolume* Construct();

    //for messenger
    void SetPolyHeight   (G4double);
    void SetFilterSpacing   (G4double);
    void SetMultiplierLeadRadius   (G4double);
    void SetModeratorAluminumRadius   (G4double);
    void SetMultiplierLeadHeightRear   (G4double);
    void SetFilterCellSpacing   (G4double);
    void SetModeratorTitaniumHeight   (G4double);
    void SetModeratorAluminumHeight   (G4double);
    void SetMultiplierLeadHeightFront   (G4double);
    void SetModeratorTitaniumRadius   (G4double);
    void SetTestX   (G4double);
    void SetTestY   (G4double);
    void SetTestZ   (G4double);
    ///


    // get methods
    //
    const G4VPhysicalVolume* GetboratedwaterPV() const;
    const G4VPhysicalVolume* GetcollimationholePV() const;
    const G4VPhysicalVolume* GetPhantomPV() const;
    //const G4VPhysicalVolume* GetPhantom2PV() const;
    //const G4VPhysicalVolume* GetPhantom3PV() const;
    //const G4VPhysicalVolume* GetPhantom4PV() const;
    ////const G4VPhysicalVolume* GetPhantom5PV() const;
    //const G4VPhysicalVolume* GetPhantom6PV() const;
    //const G4VPhysicalVolume* GetPhantom7PV() const;
    //const G4VPhysicalVolume* GetPhantom8PV() const;
    const G4VPhysicalVolume* GetLabFloorWallsolidPV() const;
    const G4VPhysicalVolume* GetLabFloorExtendedsolidPV() const;
    const G4VPhysicalVolume* GetfrontglassdoorPV() const;
    const G4VPhysicalVolume* GetfrontdoorPV() const;
    const G4VPhysicalVolume* GetglasswindowPV() const;
    const G4VPhysicalVolume* GetreardoorPV() const;
    const G4VPhysicalVolume* GetTestSurfacesolidPV() const;
    const G4VPhysicalVolume* GetmultiplierleadPV() const;
    const G4VPhysicalVolume* GetTestRIGHTSIDEPV() const;
    const G4VPhysicalVolume* GetTestTOPPV() const;
    const G4VPhysicalVolume* GetTestBOTTOMPV() const;
    const G4VPhysicalVolume* GetTestLEFTSIDEPV() const;
    const G4VPhysicalVolume* GetmoderatoraluminumPV() const;
    const G4VPhysicalVolume* GetmoderatortitaniumPV() const;
    const G4VPhysicalVolume* GetTestREARSIDEPV() const;
    const G4VPhysicalVolume* GetinnerBPolyPV() const;
    const G4VPhysicalVolume* GetTestFRONTSIDEPV() const;
    //const G4VPhysicalVolume* GetTestCENTERPOINTPV() const;
    const G4VPhysicalVolume* GetfilterscandiumPV() const; //Test_LEFTSIDE
    //const G4VPhysicalVolume* GetLiFsolidPV() const;

    //for Messenger
    G4double GetPolyHeight() const    {return fPolyHeight;};
    G4double GetFilterSpacing() const    {return fFilterSpacing;};
    G4double GetMultiplierLeadRadius() const    {return fMultiplierLeadRadius;};
    G4double GetModeratorAluminumRadius() const    {return fModeratorAluminumRadius;};
    G4double GetMultiplierLeadHeightRear() const    {return fMultiplierLeadHeightRear;};
    G4double GetFilterCellSpacing() const    {return fFilterCellSpacing;};
    G4double GetModeratorTitaniumHeight() const    {return fModeratorTitaniumHeight;};
    G4double GetModeratorAluminumHeight() const    {return fModeratorAluminumHeight;};
    G4double GetMultiplierLeadHeightFront() const    {return fMultiplierLeadHeightFront;};
    G4double GetModeratorTitaniumRadius() const    {return fModeratorTitaniumRadius;};
    G4double GetTestX() const    {return ftestx;};
    G4double GetTestY() const    {return ftesty;};
    G4double GetTestZ() const    {return ftestz;};
    ///
  private:
    // methods

    //for messenger
    IronFilterDetectorMessenger* fDetectorMessenger;
    //


    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();

    G4VPhysicalVolume* boratedwater_PV;
    G4VPhysicalVolume* collimation_hole_PV;
    G4VPhysicalVolume* Phantom_PV;
    //G4VPhysicalVolume* Phantom2_PV;
    //G4VPhysicalVolume* Phantom3_PV;
    //G4VPhysicalVolume* Phantom4_PV;
    ////G4VPhysicalVolume* Phantom5_PV;
    //G4VPhysicalVolume* Phantom6_PV;
    //G4VPhysicalVolume* Phantom7_PV;
    //G4VPhysicalVolume* Phantom8_PV;
    G4VPhysicalVolume* LabFloorWall_solid_PV;
    G4VPhysicalVolume* LabFloorExtended_solid_PV;
    G4VPhysicalVolume* frontglassdoor_PV;
    G4VPhysicalVolume* frontdoor_PV;
    G4VPhysicalVolume* glasswindow_PV;
    G4VPhysicalVolume* reardoor_PV;
    G4VPhysicalVolume* TestSurface_solid_PV;
    G4VPhysicalVolume* multiplier_lead_PV;
    G4VPhysicalVolume* Test_RIGHTSIDE_PV;
    G4VPhysicalVolume* Test_TOP_PV;
    G4VPhysicalVolume* Test_BOTTOM_PV;
    G4VPhysicalVolume* Test_LEFTSIDE_PV;
    G4VPhysicalVolume* moderator_aluminum_PV;
    G4VPhysicalVolume* moderator_titanium_PV;
    //G4VPhysicalVolume* LiF_solid_PV;
    G4VPhysicalVolume* Test_REARSIDE_PV;
    G4VPhysicalVolume* inner_BPoly_PV;
    G4VPhysicalVolume* Test_FRONTSIDE_PV;
    //G4VPhysicalVolume* Test_CENTERPOINT_PV;
    G4VPhysicalVolume* filter_scandium_PV;


    //for messenger, these are declared here but values are set in src file
    //constant parameters
    G4double delta;// 1.0cm parameter of catchment area
    G4double zeroRadius;
    G4double startAngle;
    G4double spanningAngle;
    G4double DD_Height;
    G4double RoomLength;


    //variables that you can vary using messenger
    G4double  fPolyHeight;
    G4double  fFilterSpacing;
    G4double  fMultiplierLeadRadius;
    G4double  fModeratorAluminumRadius;
    G4double  fMultiplierLeadHeightRear;
    G4double  fFilterCellSpacing;
    G4double  fModeratorTitaniumHeight;
    G4double  fModeratorAluminumHeight;
    G4double  fMultiplierLeadHeightFront;
    G4double  fModeratorTitaniumRadius;
    G4double  ftestx;
    G4double  ftesty;
    G4double  ftestz;

    //

    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps
};

// inline functions

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetLabFloorWallsolidPV() const {
  return LabFloorWall_solid_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetLabFloorExtendedsolidPV() const {
  return LabFloorExtended_solid_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetfrontglassdoorPV() const {
  return frontglassdoor_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetfrontdoorPV() const {
  return frontdoor_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetglasswindowPV() const {
  return glasswindow_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetreardoorPV() const {
  return reardoor_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetTestSurfacesolidPV() const {
  return TestSurface_solid_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetboratedwaterPV() const {
  return boratedwater_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetinnerBPolyPV() const {
  return inner_BPoly_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetmultiplierleadPV() const {
  return multiplier_lead_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetmoderatoraluminumPV() const {
  return moderator_aluminum_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetmoderatortitaniumPV() const {
  return moderator_titanium_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetfilterscandiumPV() const {
  return filter_scandium_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetcollimationholePV() const {
  return collimation_hole_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetPhantomPV() const {
  return Phantom_PV;
}

//inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetPhantom2PV() const {
//  return Phantom2_PV;
//}

//inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetPhantom3PV() const {
//  return Phantom3_PV;
//}

//inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetPhantom4PV() const {
//  return Phantom4_PV;
//}

////inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetPhantom5PV() const {
////  return Phantom5_PV;
////}

//inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetPhantom6PV() const {
//  return Phantom6_PV;
//}

//inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetPhantom7PV() const {
//  return Phantom7_PV;
//}

//inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetPhantom8PV() const {
//  return Phantom8_PV;
//}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetTestTOPPV() const {
  return Test_TOP_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetTestBOTTOMPV() const {
  return Test_BOTTOM_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetTestLEFTSIDEPV() const {
  return Test_LEFTSIDE_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetTestRIGHTSIDEPV() const {
  return Test_RIGHTSIDE_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetTestREARSIDEPV() const {
  return Test_REARSIDE_PV;
}

inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetTestFRONTSIDEPV() const {
  return Test_FRONTSIDE_PV;
}

//inline const G4VPhysicalVolume* IronFilterDetectorConstruction::GetTestCENTERPOINTPV() const {
//  return Test_CENTERPOINT_PV;
//}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
