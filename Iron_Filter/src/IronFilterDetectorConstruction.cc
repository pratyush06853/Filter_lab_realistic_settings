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
// $Id: IronFilterDetectorConstruction.cc $
//
/// \file IronFilterDetectorConstruction.cc
/// \brief Implementation of the IronFilterDetectorConstruction class

#include "IronFilterDetectorConstruction.hh"
#include "IronFilterDetectorMessenger.hh"//pratyush
#include "G4Material.hh"
#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4UnitsTable.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4NeutronHPThermalScatteringNames.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4AutoDelete.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4ios.hh"
#include "G4GeometryManager.hh"//pratyush
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IronFilterDetectorConstruction::IronFilterDetectorConstruction()
 : G4VUserDetectorConstruction(),
   //LabFloorWall_solid_PV(0),
   boratedwater_PV(0),
   multiplier_lead_PV(0),
   moderator_aluminum_PV(0),
   moderator_titanium_PV(0),
   filter_scandium_PV(0),
   collimation_hole_PV(0),
   Test_RIGHTSIDE_PV(0),
   Test_REARSIDE_PV(0),
   inner_BPoly_PV(0),
   Test_FRONTSIDE_PV(0),
   Test_CENTERPOINT_PV(0),
   Test_TOP_PV(0),
   Test_BOTTOM_PV(0),
   Test_LEFTSIDE_PV(0),
   fCheckOverlaps(true)
{

  delta= 1.0*cm;// 1.0cm parameter of catchment area
  zeroRadius = 0.*cm;
  startAngle = 0.*deg;
  spanningAngle = 360.*deg;
  DD_Height = 20.0*cm;

  //for messenger, you vary these variables using the macro
  fMultiplierLeadHeightRear = 20.0*cm;//20.0*cm
  fMultiplierLeadHeightFront=  15.0*cm;//30.0*cm

  fModeratorAluminumHeight= 30.0*cm;//39.0*cm
  fModeratorTitaniumHeight= (5.0)*cm; //(34.0)*cm

  fMultiplierLeadRadius = 15.0*cm;//20.0*cm
  fModeratorAluminumRadius = 15.0*cm;//15.0*cm without Ti
  fModeratorTitaniumRadius = 15.0*cm;

  fPolyHeight = 41.0*cm;//

  fFilterCellSpacing= 50.0*cm;//5



  fDetectorMessenger = new IronFilterDetectorMessenger(this);
  //
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IronFilterDetectorConstruction::~IronFilterDetectorConstruction()
{
  delete fDetectorMessenger;//Pratyush
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* IronFilterDetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void IronFilterDetectorConstruction::DefineMaterials()
{
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;
  G4double density, fractionMass;
  G4String symbol, name;
  G4int nComponents, nAtoms;
  G4double temp;
  G4NistManager* NistMgr = G4NistManager::Instance();

  //G4Element* elBe = new G4Element(name = "Beryllium", symbol = "Be", z = 4.0, a = 9.012*g/mole);
  G4Element* elF = new G4Element(name= "Fluorine", symbol = "F", z = 9.0, a= 18.998403*g/mole); //pratyush
  G4Element* elLi6  = new G4Element(name = "Lithium", symbol = "Li", z = 3.0, a = 6.015*g/mole);
  G4Element* elO  = new G4Element(name = "Oxygen", symbol = "O", z = 8.0, a = 15.999*g/mole);
  //G4Element* elCr  = new G4Element(name = "Chromium", symbol = "Cr", z = 24.0, a = 51.996*g/mole);
  G4Element* elFe  = new G4Element(name = "Iron", symbol = "Fe", z = 26.0, a = 55.845*g/mole);
  //G4Element* elNi  = new G4Element(name = "Nickel", symbol = "Ni", z = 28.0, a = 58.693*g/mole);
  //G4Element* elMo  = new G4Element(name = "Molybdenum", symbol = "Mo", z = 42.0, a = 95.94*g/mole);
  G4Element* elAl  = new G4Element(name = "Aluminum", symbol = "Al", z = 13.0, a = 26.982*g/mole);


/* //////////////////////////////////////////////////////////////////////////////////////////////////////////
  G4Element* elH  = new G4Element(name = "Hydrogen", symbol = "H", z = 1.0, a = 1.008*g/mole);
  G4Element* elC  = new G4Element(name = "Carbon", symbol = "C", z = 6.0, a = 12.011*g/mole);
  G4Element* elNa  = new G4Element(name = "Sodium", symbol = "Na", z = 11.0, a = 22.990*g/mole);
  G4Element* elSi  = new G4Element(name = "Silicon", symbol = "Si", z = 14.0, a = 28.085*g/mole);
  G4Element* elP  = new G4Element(name = "Phosphorus", symbol = "P", z = 15.0, a = 30.974*g/mole);
  G4Element* elK  = new G4Element(name = "Potassium", symbol = "K", z = 19.0, a = 39.098*g/mole);
  G4Element* elCa  = new G4Element(name = "Calcium", symbol = "Ca", z = 20.0, a = 40.078*g/mole);
  G4Element* elMn  = new G4Element(name = "Maganese", symbol = "Mn", z = 25.0, a = 54.938*g/mole);
  G4Element* elTi  = new G4Element(name = "Titanium", symbol = "Ti", z = 22.0, a = 47.867*g/mole);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/
  //Vacuum
  new G4Material("galactic", z = 1.0, a = 1.01*g/mole, density = universe_mean_density, kStateGas, 2.73*kelvin, 3.e-18*pascal);

  // Lead
  new G4Material(name = "Pb", z = 82.0, a = 207.2*g/mole, density = 11.34*g/cm3);

  // Iron
  new G4Material("NatIron", z = 26.0, a = 55.845*g/mole, density = 7.874*g/cm3);


  // Iron
  new G4Material("NatScandium", z = 21.0, a = 44.95*g/mole, density = 2.985*g/cm3);

  // Titanium
  new G4Material("NatTi", z = 22.0, a = 47.867*g/mole, density = 4.507*g/cm3);

  //Natural Boron
  //G4Material* NatB = new G4Material("NatB", z = 5.0, a = 10.811*g/mole, density = 2.37*g/cm3);
  G4Element* NatB=NistMgr->FindOrBuildElement("B");

  //Natural Carbon
  //G4Material* NatC = new G4Material("NatC", z = 6.0, a = 12.0107*g/mole, density = 1.95*g/cm3);
  G4Element* NatC=NistMgr->FindOrBuildElement("C");

  //Natural Hydrogen
  //G4Material* NatH = new G4Material("NatH", z = 1.0, a = 1.00784*g/mole, density = 0.0008999*g/cm3);
  G4Element* NatH=NistMgr->FindOrBuildElement("H");

  //Natural Oxygen
  //G4Material* NatO = new G4Material("NatO", z = 8.0, a = 15.999*g/mole, density = 0.001141*g/cm3);
  G4Element* NatO=NistMgr->FindOrBuildElement("O");

  //Natural Oxygen
  //G4Material* NatNa = new G4Material("NatNa", z = 11.0, a = 22.989*g/mole, density = 0.968*g/cm3);
  G4Element* NatNa=NistMgr->FindOrBuildElement("Na");

   //Natural Sulphur
   //G4Material* NatS = new G4Material("NatS", z = 16.0, a = 32.065*g/mole, density = 1.92*g/cm3);
   G4Element* NatS=NistMgr->FindOrBuildElement("S");

   //Ni60
   G4Material* Ni60 = new G4Material("Ni60", z = 28.0, a = 59.9307864*g/mole, density = 8.9*g/cm3);
   //Fe54
   G4Material* Fe54 = new G4Material("Fe54", z = 26.0, a = 	53.9396090*g/mole, density = 7.874*g/cm3);
   //Co
   G4Material* NatCo = new G4Material("NatCo", z = 27.0, a = 	58.933*g/mole, density = 8.86*g/cm3);

   G4Element *TS_H_P = new G4Element("TS_H_of_Polyethylene", "H", 1, 1.007*g/mole);
   G4Element *TS_H_W = new G4Element("TS_H_of_Water", "H", 1, 1.007*g/mole);

   G4Material* air=NistMgr->FindOrBuildMaterial("G4_AIR");
   //G4Material* water=NistMgr->FindOrBuildMaterial("G4_WATER");

   //Water
   G4Material* water = new G4Material("water", density= 1.00 * g / cm3,nComponents= 2, kStateLiquid, 296*kelvin);
   water->AddElement(TS_H_W,2);//pratyush
   water->AddElement(elO,1);//pratyush

  //Lithium6_Fluoride
  G4Material* Li6F = new G4Material("Li6F", density= 2.54 * g / cm3,nComponents= 2);
  Li6F->AddElement(elLi6, 1); //pratyush
  Li6F->AddElement(elF,1);   //pratyush

  //Aluminum Fluoride
  G4Material* AlF3 = new G4Material("AlF3", density= 3.10 * g / cm3,nComponents= 2); //pratyush
  AlF3->AddElement(elAl, 1);  //pratyush
  AlF3->AddElement(elF, 3);   //pratyush

  //Aluminum Fluoride
  G4Material* Borax = new G4Material("Borax", density= 0.76* g / cm3,nComponents= 4,kStateSolid, 296*kelvin); //pratyush
  Borax->AddElement(NatNa,12.06*perCent);//2
  Borax->AddElement(NatB,11.34*perCent);//4
  Borax->AddElement(NatH,5.29*perCent);//20
  Borax->AddElement(NatO,71.32*perCent);//17

  //Boric Acid https://www.convertunits.com/molarmass/Boric+Acid
  G4Material* boric_acid = new G4Material("boric_acid", density= 1.44* g / cm3,nComponents= 3,kStateSolid, 296*kelvin); //pratyush
  boric_acid->AddElement(NatH,4.890*perCent);//3
  boric_acid->AddElement(NatB,17.484*perCent);//1
  boric_acid->AddElement(NatO,77.626*perCent);//3

  //Borax_water_Mixture(5.8% solubity of borax https://omsi.edu/sites/all/FTP/files/kids/Borax-msds.pdf)
  //mixture of 5.5% Borax and 94.5% of Water
  G4Material* borax_water = new G4Material( "borax_water",density= 0.9868*g/cm3, nComponents= 2,kStateLiquid, 296*kelvin); //pratyush
  borax_water->AddMaterial( Borax, 5.5*perCent );  //pratyush
  borax_water->AddMaterial( water, 94.5*perCent ); //pratyush


  //Borax_BoricAcid_buffer(https://www.researchgate.net/publication/244069630_Preparation_of_highly_concentrated_aqueous_solution_of_sodium_borate)
  //mixture of 20g BoricAcid, 25g of Borax and 100g of water
  G4Material* borax_boricacid_buffer = new G4Material( "borax_boricacid_buffer",density= 1.019*g/cm3, nComponents= 3,kStateLiquid, 296*kelvin);
  borax_boricacid_buffer->AddMaterial( boric_acid, 13.7*perCent );//pratyush
  borax_boricacid_buffer->AddMaterial( Borax, 17.2*perCent );//pratyush
  borax_boricacid_buffer->AddMaterial( water, 69.1*perCent );//pratyush

  //Fluental
  //mixture of 40% Al and 60% of AlF_3
  G4Material* fluental = new G4Material( "fluental",density= 2.94*g/cm3, nComponents= 2); //pratyush
  fluental->AddMaterial( AlF3, 60.*perCent );  //pratyush
  fluental->AddElement( elAl, fractionMass = 40.*perCent ); //pratyush



  //polyethyleneBorated
  G4Material* boratedPoly = new G4Material( "boratedPoly", density=1.19*g/cm3, nComponents=3,kStateSolid, 296*kelvin);
  //boratedPoly->AddElement( NatB, 3.*perCent );
  //boratedPoly->AddElement( NatC, 81.4*0.97*perCent );
  //boratedPoly->AddElement(NatH, 18.6*0.97*perCent );
  boratedPoly->AddElement( NatB, 3.*perCent );
  boratedPoly->AddElement( NatC, 14.424*perCent );
  boratedPoly->AddElement(TS_H_P, 82.576*perCent );

  // Polyethylene with boron at 3% - Has borated polyethylene any oxygen elements?
  //G4Material* boratedPoly = new G4Material("boratedPoly ", density=0.96*g/cm3, nComponents=3);
  //boratedPoly->AddElement(NatH,  14.424*perCent);
  //boratedPoly->AddElement(NatC,  82.576*perCent);
  //boratedPoly->AddElement(NatB,  3.00*perCent);




  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* IronFilterDetectorConstruction::DefineVolumes()
{

  // Get materials
  G4Material* Vacuum = G4Material::GetMaterial("galactic");
  G4Material* Iron = G4Material::GetMaterial("NatIron");
  G4Material* Scandium = G4Material::GetMaterial("NatScandium"); //
  G4Material* Titanium = G4Material::GetMaterial("NatTi");
  G4Material* Sulphur = G4Material::GetMaterial("NatS");
  G4Material* Fluental = G4Material::GetMaterial("fluental");
  G4Material* BoraxWater = G4Material::GetMaterial("borax_water");
  G4Material* BoraxBoricAcidBuffer = G4Material::GetMaterial("borax_boricacid_buffer");
  G4Material* Lithium6_Fluoride = G4Material::GetMaterial("Li6F");
  G4Material* Lead = G4Material::GetMaterial("Pb");
  G4Material* BoratedPoly = G4Material::GetMaterial("boratedPoly");
  G4Material* Cobalt = G4Material::GetMaterial("NatCo");
  G4Material* Nickel60 = G4Material::GetMaterial("Ni60");
  G4Material* Iron54 = G4Material::GetMaterial("Fe54");
  G4Material* Air=G4Material::GetMaterial("air");

  if ( ! Vacuum ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined.";
    G4Exception("IronFilterDetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }

//
//  calculate the Sizes of some derived units.
//
G4double Scandium_diameter_limited=3.5*cm;//3.5*cm;
G4double Scandium_height_limited=35*cm;//3.5*cm;
//G4double Pb_radius = fSource_radius + 5.0*cm ;


G4double NeutronFilter_length = fMultiplierLeadHeightRear+fMultiplierLeadHeightFront+fModeratorAluminumHeight+fModeratorTitaniumHeight+fPolyHeight;


G4double Water_rear_side = 20.0*cm ;
G4double Water_x = 200.0*cm;
G4double Water_y = NeutronFilter_length+Water_rear_side;
G4double Water_z = 200.0*cm;
G4double Poly_a = 40.0*cm;
G4double hole_length = NeutronFilter_length-(fMultiplierLeadHeightRear+fMultiplierLeadHeightFront)-fModeratorAluminumHeight-fModeratorTitaniumHeight-Scandium_height_limited;

G4double lab68_wall_thickness = 25.0*cm ;
//distance from the outer egdes
G4double lab68_wall_x = 6.6*m ;
G4double lab68_wall_y = 9.3*m ;
G4double lab68_wall_z = 5.6*m ;

G4ThreeVector position_of_origin = {2.7*m, -2.45*m, 1.3*m}; //with repect to the inner upper left corner of the room


//
// Rotations
//

  G4RotationMatrix* NO_ROT = new G4RotationMatrix;

  G4RotationMatrix* turnAlongX = new G4RotationMatrix;
  turnAlongX->rotateX(90*deg);

  G4RotationMatrix* turnAlongXY = new G4RotationMatrix;
  turnAlongXY->rotateZ(90*deg);
  turnAlongXY->rotateZ(-100*deg);

  G4RotationMatrix* turnAlong = new G4RotationMatrix;
  turnAlong->rotateZ(10*deg);
  //turnAlong->rotateY(10*deg);

  G4RotationMatrix* turnAlongZ = new G4RotationMatrix;
  turnAlongZ->rotateZ(90*deg);
  turnAlongZ->rotateZ(-90*deg);



//
// GEOMETRY
//

  // Start with a  vacuum layer. This will be the mother volume

  G4VSolid* vacuum_solid =new G4Box("vacuum_solid", 150.0*m, 150.0*m, 150.0*m);
  G4LogicalVolume* vacuum_solid_LV = new G4LogicalVolume(vacuum_solid, Vacuum, "vacuum_solid");
  G4VPhysicalVolume* vacuum_solid_PV = new G4PVPlacement(NO_ROT,G4ThreeVector{0.,0.,0.}, vacuum_solid_LV, "Vacuum_solid", 0, false, 0, fCheckOverlaps);
  vacuum_solid_LV->SetVisAttributes(G4VisAttributes::Invisible);

  //out_neutron donot include ceiling
  G4VSolid* Main_2_S = new G4Box("Main_2_solid", lab68_wall_x/2.0, lab68_wall_y/2.0 , lab68_wall_z/2.0);
  G4VSolid* hole_2_S = new G4Box("hole_2_solid", (lab68_wall_x-2*lab68_wall_thickness)/2.0, (lab68_wall_y-2*lab68_wall_thickness)/2.0, (lab68_wall_z-2*lab68_wall_thickness)/2.0);
  G4SubtractionSolid* LabFloorWall_solid_S= new G4SubtractionSolid("LabFloorWall_solid", Main_2_S, hole_2_S, NO_ROT, G4ThreeVector(0.,0., 0.));
  G4LogicalVolume* LabFloorWall_solid_LV = new G4LogicalVolume(LabFloorWall_solid_S, Vacuum, "LabFloorWall_solid");
  //LabFloorWall_solid_PV = new G4PVPlacement(turnAlong, G4ThreeVector{lab68_wall_x/2.0,-lab68_wall_y/2.0,lab68_wall_z/2.0}-position_of_origin, LabFloorWall_solid_LV, "OutSpacer", vacuum_solid_LV, false, 0, fCheckOverlaps);
  //LabFloorWall_solid_LV->SetVisAttributes(G4VisAttributes(G4Colour::Grey()));
  //LabFloorWall_solid_LV->SetVisAttributes(G4VisAttributes::Invisible);


  //layer of air at the bottom of the filter
  G4VSolid* Main_S = new G4Box("Main_solid", Water_x/2.0, Water_y/2.0, Water_z/2.0);
  G4VSolid* hole_S = new G4Box("hole_solid", Poly_a/2.0 , Poly_a/2.0, NeutronFilter_length/2.0);
  G4SubtractionSolid* boratedwater_S= new G4SubtractionSolid("boratedwater", Main_S, hole_S, turnAlongX, G4ThreeVector(0., Water_y/2.0-Water_rear_side-NeutronFilter_length/2.0, 0.));
  //G4LogicalVolume* boratedwater_LV = new G4LogicalVolume(boratedwater_S, BoraxWater, "boratedwater");
  //G4LogicalVolume* boratedwater_LV = new G4LogicalVolume(boratedwater_S, BoraxBoricAcidBuffer, "boratedwater");
  G4LogicalVolume* boratedwater_LV = new G4LogicalVolume(boratedwater_S, BoratedPoly, "boratedwater");
  boratedwater_PV = new G4PVPlacement(turnAlongZ, G4ThreeVector(0., fFilterCellSpacing+Water_y/2.0, 0), boratedwater_LV, "BoratedWater", vacuum_solid_LV, false, 0, fCheckOverlaps);
  boratedwater_LV->SetVisAttributes(G4VisAttributes(G4Colour::Cyan()));


  // inner Shield made of of Borated Poly
  G4VSolid* Main_1_S = new G4Box("Main_1_solid", Poly_a/2.0 , Poly_a/2.0, NeutronFilter_length/2.0);
  G4VSolid* hole_1_S = new G4Tubs("hole_1_solid", 0 , Scandium_diameter_limited/2.0, hole_length/2.0,startAngle, spanningAngle);
  G4SubtractionSolid* inner_BPoly_S= new G4SubtractionSolid("inner_BPoly_solid", Main_1_S, hole_1_S, NO_ROT, G4ThreeVector(0., 0, -NeutronFilter_length/2.0+hole_length/2.0));
  G4LogicalVolume *inner_BPoly_LV = new G4LogicalVolume(inner_BPoly_S, BoratedPoly,"inner_BPoly" );
  inner_BPoly_PV = new G4PVPlacement( turnAlongX, G4ThreeVector(0., fFilterCellSpacing+NeutronFilter_length/2.0, 0.), inner_BPoly_LV, "inner_BPoly", vacuum_solid_LV, false, 0, fCheckOverlaps);
  inner_BPoly_LV->SetVisAttributes(G4VisAttributes(G4Colour::Yellow()));


  //lead in the form of cylinder as outer shield
  G4VSolid* multiplier_lead_S = new G4Tubs("multiplier_lead", zeroRadius, fMultiplierLeadRadius , (fMultiplierLeadHeightRear+fMultiplierLeadHeightFront)/2.0, startAngle, spanningAngle);
  G4LogicalVolume* multiplier_lead_LV = new G4LogicalVolume(multiplier_lead_S, Lead, "multiplier_lead");
  multiplier_lead_PV = new G4PVPlacement(NO_ROT, G4ThreeVector(0., 0.,NeutronFilter_length/2.0-(fMultiplierLeadHeightRear+fMultiplierLeadHeightFront)/2.0), multiplier_lead_LV, "Pb_shield", inner_BPoly_LV, false, 0, fCheckOverlaps);
  multiplier_lead_LV->SetVisAttributes(G4VisAttributes(G4Colour::Blue()));

  //Fluental (Aluminum) Base
  G4VSolid* moderator_aluminum_S = new G4Tubs("moderator_aluminum", zeroRadius, fModeratorAluminumRadius, (fModeratorAluminumHeight /2.0), startAngle, spanningAngle);
  G4LogicalVolume* moderator_aluminum_LV = new G4LogicalVolume(moderator_aluminum_S, Fluental, "moderator_aluminum");
  moderator_aluminum_PV =new G4PVPlacement(NO_ROT, G4ThreeVector(0,0,NeutronFilter_length/2.0-(fMultiplierLeadHeightRear+fMultiplierLeadHeightFront)-fModeratorAluminumHeight /2.0),moderator_aluminum_LV, "Fluental_Moderator", inner_BPoly_LV, false, 0, fCheckOverlaps);
  moderator_aluminum_LV->SetVisAttributes(G4VisAttributes(G4Colour::Grey()));


  //Iron final  Moderator in the form of cylinder
  G4VSolid* moderator_titanium_S = new G4Tubs( "moderator_titanium", zeroRadius, fModeratorTitaniumRadius,(fModeratorTitaniumHeight/2.0), startAngle, spanningAngle);
  G4LogicalVolume* moderator_titanium_LV = new G4LogicalVolume( moderator_titanium_S, Titanium, "moderator_titanium" );//to check with fluental
  moderator_titanium_PV = new G4PVPlacement( NO_ROT, G4ThreeVector(0,0,NeutronFilter_length/2.0-(fMultiplierLeadHeightRear+fMultiplierLeadHeightFront)-fModeratorAluminumHeight-fModeratorTitaniumHeight/2.0),moderator_titanium_LV, "Moderator_titanium",inner_BPoly_LV, false, 0, fCheckOverlaps);
  moderator_titanium_LV->SetVisAttributes(G4VisAttributes(G4Colour::Green()));

  //Iron as the first moderator
  G4VSolid* filter_scandium_S = new G4Tubs("filter_scandium", zeroRadius,  Scandium_diameter_limited/2.0,(Scandium_height_limited/2.0), startAngle, spanningAngle);
  G4LogicalVolume *filter_scandium_LV = new G4LogicalVolume(filter_scandium_S, Scandium,"filter_scandium" );
  filter_scandium_PV = new G4PVPlacement( NO_ROT, G4ThreeVector(0,0,NeutronFilter_length/2.0-(fMultiplierLeadHeightRear+fMultiplierLeadHeightFront)-fModeratorAluminumHeight-fModeratorTitaniumHeight-Scandium_height_limited/2.0), filter_scandium_LV, "Filter_scandium",inner_BPoly_LV,false, 0, fCheckOverlaps);
  filter_scandium_LV->SetVisAttributes(G4VisAttributes(G4Colour::Red()));


  // Poly need to change
  G4VSolid* collimation_hole_S = new G4Tubs("collimation_hole", zeroRadius, Scandium_diameter_limited/2.0, hole_length/2.0, startAngle, spanningAngle);
  //G4LogicalVolume* collimation_hole_LV = new G4LogicalVolume(collimation_hole_S, Air, "collimation_hole");
  G4LogicalVolume* collimation_hole_LV = new G4LogicalVolume(collimation_hole_S, Vacuum, "collimation_hole");
  collimation_hole_PV = new G4PVPlacement(turnAlongX, G4ThreeVector(0., fFilterCellSpacing+hole_length/2.0, 0.0), collimation_hole_LV, "collimation_hole", vacuum_solid_LV, false, 0, fCheckOverlaps);
  collimation_hole_LV->SetVisAttributes(G4VisAttributes(G4Colour::Blue()));

  //TopSide; Daughter volume of Boraxwater Vessel
  G4VSolid* Test_TOP_S = new G4Box("Test_TOP_solid", Water_x/2.0, Water_y/2.0, delta/2.0);
  G4LogicalVolume *Test_TOP_LV = new G4LogicalVolume(Test_TOP_S, Vacuum,"Test_TOP" );
  Test_TOP_PV = new G4PVPlacement( NO_ROT, G4ThreeVector(0., 0, Water_z/2.0-delta/2.0), Test_TOP_LV, "Test_TOP", boratedwater_LV, false, 0, fCheckOverlaps);
  Test_TOP_LV->SetVisAttributes(G4VisAttributes(G4Colour::Red()));

  //Bottomside; Daughter volume of Boraxwater Vessel
  G4VSolid* Test_BOTTOM_S = new G4Box("Test_BOTTOM_solid", Water_x/2.0, Water_y/2.0, delta/2.0);
  G4LogicalVolume *Test_BOTTOM_LV = new G4LogicalVolume(Test_BOTTOM_S, Vacuum,"Test_BOTTOM" );
  Test_BOTTOM_PV = new G4PVPlacement( NO_ROT, G4ThreeVector(0., 0, -Water_z/2.0+delta/2.0), Test_BOTTOM_LV, "Test_BOTTOM", boratedwater_LV, false, 0, fCheckOverlaps);
  Test_BOTTOM_LV->SetVisAttributes(G4VisAttributes(G4Colour::Red()));


  //Sides
  //Leftside; Daughter volume of Boraxwater Vessel
  G4VSolid* Test_LEFTSIDE_S = new G4Box("Test_LEFTSIDE_solid",delta/2.0, Water_y/2.0, Water_z/2.0-delta);
  G4LogicalVolume *Test_LEFTSIDE_LV = new G4LogicalVolume(Test_LEFTSIDE_S, Vacuum,"Test_LEFTSIDE" );
  Test_LEFTSIDE_PV = new G4PVPlacement( NO_ROT, G4ThreeVector( Water_x/2.0-delta/2.0, 0,0), Test_LEFTSIDE_LV, "Test_LEFTSIDE", boratedwater_LV, false, 0, fCheckOverlaps);
  Test_LEFTSIDE_LV->SetVisAttributes(G4VisAttributes(G4Colour::Red()));


  //Rightside; Daughter volume of Boraxwater Vessel
  G4VSolid* Test_RIGHTSIDE_S = new G4Box("Test_RIGHTSIDE_solid", delta/2.0, Water_y/2.0, Water_z/2.0-delta);
  G4LogicalVolume *Test_RIGHTSIDE_LV = new G4LogicalVolume(Test_RIGHTSIDE_S, Vacuum,"Test_RIGHTSIDE" );
  Test_RIGHTSIDE_PV = new G4PVPlacement( NO_ROT, G4ThreeVector( -Water_x/2.0+delta/2.0, 0,0), Test_RIGHTSIDE_LV, "Test_RIGHTSIDE", boratedwater_LV, false, 0, fCheckOverlaps);
  Test_RIGHTSIDE_LV->SetVisAttributes(G4VisAttributes(G4Colour::Red()));

  //RearSide; Daughter Volume of Boraxwater vessel
  G4VSolid* Test_REARSIDE_S = new G4Box("Test_REARSIDE_solid", Water_x/2.0-delta, delta/2.0, Water_z/2.0-delta);
  G4LogicalVolume *Test_REARSIDE_LV = new G4LogicalVolume(Test_REARSIDE_S, Vacuum,"Test_REARSIDE" );
  Test_REARSIDE_PV = new G4PVPlacement( NO_ROT, G4ThreeVector( 0, Water_y/2.0-delta/2.0,0), Test_REARSIDE_LV, "Test_REARSIDE", boratedwater_LV, false, 0, fCheckOverlaps);
  Test_REARSIDE_LV->SetVisAttributes(G4VisAttributes(G4Colour::Red()));


  //FrontSide; just outside the Boraxwater Vessel
  G4VSolid* Test_FRONTSIDE_S = new G4Box("Test_FRONTSIDE_solid", Water_x/2.0, delta/2.0, Water_z/2.0);
  G4LogicalVolume *Test_FRONTSIDE_LV = new G4LogicalVolume(Test_FRONTSIDE_S, Vacuum,"Test_FRONTSIDE" );
  new G4PVPlacement(turnAlongZ, G4ThreeVector(0., fFilterCellSpacing-delta/2.0,0), Test_FRONTSIDE_LV, "Test_FRONTSIDE",vacuum_solid_LV, false, 0, fCheckOverlaps);
  Test_FRONTSIDE_LV->SetVisAttributes(G4VisAttributes(G4Colour::Red()));

  //At the center
  G4VSolid* Test_CENTERPOINT_S = new G4Box("Test_CENTERPOINT_solid", Water_x/2.0, delta/2.0, Water_z/2.0);
  G4LogicalVolume *Test_CENTERPOINT_LV = new G4LogicalVolume(Test_CENTERPOINT_S, Vacuum,"Test_CENTERPOINT" );
  new G4PVPlacement(turnAlongZ, G4ThreeVector(0., 0., 0.), Test_CENTERPOINT_LV, "Test_CENTERPOINT", vacuum_solid_LV, false, 0, fCheckOverlaps);
  Test_CENTERPOINT_LV->SetVisAttributes(G4VisAttributes(G4Colour::Red()));


  // Always return the physical World

  return vacuum_solid_PV;
}


//method SetPolyHeight is defined here
void IronFilterDetectorConstruction::SetPolyHeight(G4double ival)
{
  if (ival < 1)
    { G4cout << "\n --->warning from SetfPolyHeight: "
             << ival << " must be at least 1. Command refused" << G4endl;
      return;
    }
  fPolyHeight = ival;
}


//method SetFilterSpacing is defined here
void IronFilterDetectorConstruction::SetFilterSpacing(G4double ival)
{
  if (ival < 1)
    { G4cout << "\n --->warning from SetfFilterSpacing: "
             << ival << " must be at least 1. Command refused" << G4endl;
      return;
    }
  fFilterSpacing = ival;
}

//method SetMultiplierLeadRadius is defined here
void IronFilterDetectorConstruction::SetMultiplierLeadRadius(G4double ival)
{
  if (ival < 1)
    { G4cout << "\n --->warning from SetfMultiplierLeadRadius: "
             << ival << " must be at least 1. Command refused" << G4endl;
      return;
    }
  fMultiplierLeadRadius = ival;
}
//method SetModeratorAluminumRadius is defined here
void IronFilterDetectorConstruction::SetModeratorAluminumRadius(G4double ival)
{
  if (ival < 1)
    { G4cout << "\n --->warning from SetfModeratorAluminumRadius: "
             << ival << " must be at least 1. Command refused" << G4endl;
      return;
    }
  fModeratorAluminumRadius = ival;
}
//method SetMultiplierLeadHeightRear is defined here
void IronFilterDetectorConstruction::SetMultiplierLeadHeightRear(G4double ival)
{
  if (ival < 1)
    { G4cout << "\n --->warning from SetfMultiplierLeadHeightRear: "
             << ival << " must be at least 1. Command refused" << G4endl;
      return;
    }
  fMultiplierLeadHeightRear = ival;
}
//method SetFilterCellSpacing is defined here
void IronFilterDetectorConstruction::SetFilterCellSpacing(G4double ival)
{
  if (ival < 1)
    { G4cout << "\n --->warning from SetfFilterCellSpacing: "
             << ival << " must be at least 1. Command refused" << G4endl;
      return;
    }
  fFilterCellSpacing = ival;
}
//method SetModeratorTitaniumHeight is defined here
void IronFilterDetectorConstruction::SetModeratorTitaniumHeight(G4double ival)
{
  if (ival < 1)
    { G4cout << "\n --->warning from SetfModeratorTitaniumHeight: "
             << ival << " must be at least 1. Command refused" << G4endl;
      return;
    }
  fModeratorTitaniumHeight = ival;
}
//method SetModeratorAluminumHeight is defined here
void IronFilterDetectorConstruction::SetModeratorAluminumHeight(G4double ival)
{
  if (ival < 1)
   { G4cout << "\n --->warning from SetfModeratorAluminumHeight: "
             << ival << " must be at least 1. Command refused" << G4endl;
      return;
    }
  fModeratorAluminumHeight = ival;
}
//method SetMultiplierLeadHeightFront is defined here
void IronFilterDetectorConstruction::SetMultiplierLeadHeightFront(G4double ival)
{
  if (ival < 1)
    { G4cout << "\n --->warning from SetfMultiplierLeadHeightFront: "
             << ival << " must be at least 1. Command refused" << G4endl;
      return;
    }
  fMultiplierLeadHeightFront = ival;
}
//method SetModeratorTitaniumRadius is defined here
void IronFilterDetectorConstruction::SetModeratorTitaniumRadius(G4double ival)
{
  if (ival < 1)
    { G4cout << "\n --->warning from SetfModeratorTitaniumRadius: "
             << ival << " must be at least 1. Command refused" << G4endl;
      return;
    }
  fModeratorTitaniumRadius = ival;
}
