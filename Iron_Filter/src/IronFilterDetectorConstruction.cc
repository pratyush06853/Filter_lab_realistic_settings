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
#include <TMath.h>

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

IronFilterDetectorConstruction::IronFilterDetectorConstruction()
 : G4VUserDetectorConstruction(),
   LabFloorWall_solid_PV(0),
   LabFloorExtended_solid_PV(0),
   frontglassdoor_PV(0),
   frontdoor_PV(0),
   glasswindow_PV(0),
   reardoor_PV(0),
   lab66FloorWall_solid_PV(0),
   lab70FloorWall_solid_PV(0),
   Lab32FloorWall_solid_PV(0),
   lab30FloorWall_solid_PV(0),
   lab28FloorWall_solid_PV(0),
   TestSurface_solid_PV(0),
   boratedwater_PV(0),
   multiplier_lead_PV(0),
   moderator_aluminum_PV(0),
   moderator_titanium_PV(0),
   filter_scandium_PV(0),
   collimation_hole_PV(0),
   Phantom_PV(0),
   //Phantom2_PV(0),
   //Phantom3_PV(0),
   //Phantom4_PV(0),
   ////Phantom5_PV(0),
   //Phantom6_PV(0),
   //Phantom7_PV(0),
   //Phantom8_PV(0),
   Test_RIGHTSIDE_PV(0),
   Test_REARSIDE_PV(0),
   inner_BPoly_PV(0),
   Test_FRONTSIDE_PV(0),
   //Test_CENTERPOINT_PV(0),
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

  ftestx = 8*m;
  ftesty = 8*m;
  ftestz = 8*m;



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
  //G4Element* elO  = new G4Element(name = "Oxygen", symbol = "O", z = 8.0, a = 15.999*g/mole);
  //G4Element* elCr  = new G4Element(name = "Chromium", symbol = "Cr", z = 24.0, a = 51.996*g/mole);
  //G4Element* elFe  = new G4Element(name = "Iron", symbol = "Fe", z = 26.0, a = 55.845*g/mole);
  //G4Element* elNi  = new G4Element(name = "Nickel", symbol = "Ni", z = 28.0, a = 58.693*g/mole);
  //G4Element* elMo  = new G4Element(name = "Molybdenum", symbol = "Mo", z = 42.0, a = 95.94*g/mole);
  G4Element* elAl  = new G4Element(name = "Aluminum", symbol = "Al", z = 13.0, a = 26.982*g/mole);
  //G4Element* elZn  = new G4Element(name = "zinc", symbol = "Zn", z = 30.0, a = 65.38*g/mole);
  G4Element* elS  = new G4Element(name = "sulphur", symbol = "S", z = 16.0, a = 32.065*g/mole);
  G4Element* elH  = new G4Element(name = "Hydrogen", symbol = "H", z = 1.0, a = 1.008*g/mole);
  G4Element* elC  = new G4Element(name = "Carbon", symbol = "C", z = 6.0, a = 12.011*g/mole);
  G4Element* elN = new G4Element("Nitrogen", symbol = "N",z = 7.,a = 14.01*g/mole);
  G4Element* elO = new G4Element("Oxygen", symbol = "O",z = 8.,a = 16.00*g/mole);

  G4Element* elNa = new G4Element("Sodium",symbol ="Na",z = 11.,a= 22.99*g/mole);

  G4Element* elMg = new G4Element("Magnesium",symbol ="Mg",z = 12.,a= 24.305*g/mole);

  G4Element* elP = new G4Element("Phosphorus",symbol ="P",z = 15.,a= 30.974*g/mole);

  G4Element* elCl = new G4Element("Chlorine",symbol ="Cl",z = 17.,a= 35.453*g/mole);

  G4Element* elK = new G4Element("Potassium",symbol ="K",z = 19.,a= 39.098*g/mole);

  G4Element* elCa = new G4Element("Calcium",symbol ="Ca",z = 20.,a= 40.08*g/mole);

  G4Element* elFe  = new G4Element("Iron",symbol ="Fe",z = 26.,a= 55.85*g/mole);

  G4Element* elZn = new G4Element("Zinc",symbol ="Zn",z = 30.,a= 65.38*g/mole);

  G4Element* elRb = new G4Element("Rb",symbol ="Rb",z = 37.,a= 85.47 *g/mole);

  G4Element* elSr = new G4Element("Sr",symbol ="Sr",z = 38.,a= 87.62 *g/mole);

  G4Element* elZr = new G4Element("Zr",symbol ="Zr",z = 40.,a= 91.22 *g/mole);

  G4Element* elPb = new G4Element("Lead",symbol ="Pb", z = 82.,a= 207.19 *g/mole);

  G4Element* elSi  = new G4Element("Silicon", symbol = "Si", z = 14.0, a = 28.085*g/mole);

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

   G4Material* air = NistMgr->FindOrBuildMaterial("G4_AIR");
   //G4Material* concrete = NistMgr->FindOrBuildMaterial("G4_CONCRETE");

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
  boratedPoly->AddElement( NatC, 82.576*perCent );
  boratedPoly->AddElement(TS_H_P, 14.424*perCent );

  // wood
  G4Material* wood = new G4Material("wood", density=0.9*g/cm3, nComponents=3);
  wood->AddElement(TS_H_P , 4);
  wood->AddElement(elO , 1);
  wood->AddElement(elC, 2);

  G4Material* quartz = new G4Material("quartz", density=2.200*g/cm3, nComponents=2);
  quartz->AddElement(elSi, 1);
  quartz->AddElement(elO , 2);

  // Polyethylene with boron at 3% - Has borated polyethylene any oxygen elements?
  //G4Material* boratedPoly = new G4Material("boratedPoly ", density=0.96*g/cm3, nComponents=3);
  //boratedPoly->AddElement(NatH,  14.424*perCent);
  //boratedPoly->AddElement(NatC,  82.576*perCent);
  //boratedPoly->AddElement(NatB,  3.00*perCent);
  G4Material*soft_tissue = new G4Material("soft_tissue",density= 0.9869*g/cm3,nComponents=9);
  //soft_tissue->AddElement(TS_H_P,0.1047*perCent);
  //soft_tissue->AddElement(elC,0.2302*perCent);
  //soft_tissue->AddElement(elN,0.0234*perCent);
  //soft_tissue->AddElement(elO,0.6321*perCent);
  //soft_tissue->AddElement(elNa,0.0013*perCent);
  //soft_tissue->AddElement(elMg,0.00015*perCent);
  //soft_tissue->AddElement(elP,0.0024*perCent);
  //soft_tissue->AddElement(elS,0.0022*perCent);
  //soft_tissue->AddElement(elCl,0.0014*perCent);
  //soft_tissue->AddElement(elK,0.0021*perCent);
  //soft_tissue->AddElement(elFe,0.000063*perCent);
  //soft_tissue->AddElement(elZn,0.000032*perCent);
  //soft_tissue->AddElement(elRb,0.0000057*perCent);
  //soft_tissue->AddElement(elSr,0.00000034*perCent);
  //soft_tissue->AddElement(elZr,0.000008*perCent);
  //soft_tissue->AddElement(elPb,0.00000016*perCent);
  soft_tissue->AddElement(elH,0.105);
  soft_tissue->AddElement(elC,0.256);
  soft_tissue->AddElement(elN,0.027);
  soft_tissue->AddElement(elO,0.602);
  soft_tissue->AddElement(elNa,0.001);
  soft_tissue->AddElement(elP,0.002);
  soft_tissue->AddElement(elS,0.003);
  soft_tissue->AddElement(elCl,0.002);
  soft_tissue->AddElement(elK,0.002);


  //soil
  G4Material*soil = new G4Material("soil",density= 1.50*g/cm3,nComponents=8);
  soil->AddElement(elH,0.021);
  soil->AddElement(elC,0.016);
  soil->AddElement(elO,0.577);
  soil->AddElement(elAl,0.050);
  soil->AddElement(elSi,0.271);
  soil->AddElement(elK,0.013);
  soil->AddElement(elCa,0.041);
  soil->AddElement(elFe,0.011);

  //concrete
  G4Material*concrete = new G4Material("concrete",density= 2.3*g/cm3,nComponents=10);
  concrete->AddElement(elH,0.01);
  concrete->AddElement(elC,0.001);
  concrete->AddElement(elO,0.529107);
  concrete->AddElement(elNa,0.016);
  concrete->AddElement(elMg,0.002);
  concrete->AddElement(elAl,0.033872);
  concrete->AddElement(elSi,0.337021);
  concrete->AddElement(elK,0.013);
  concrete->AddElement(elCa,0.044);
  concrete->AddElement(elFe,0.014);


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
  G4Material*  Soft_Tissue=G4Material::GetMaterial("soft_tissue");
  G4Material*  Concrete = G4Material::GetMaterial("concrete");
  G4Material*  Wood = G4Material::GetMaterial("wood");
  G4Material*  Quartz = G4Material::GetMaterial("quartz");
  G4Material*  Soil = G4Material::GetMaterial("soil");
  G4Material* Water = G4Material::GetMaterial("water");
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

G4double Pump_chase_y = 2.5*m ;

G4double lab68_frontdoor_glass_height = 2.9*m;
G4double lab68_frontdoor_glass_width = 2.2*m;

G4double lab68_frontdoor_wood_height = 2.3*m;
G4double lab68_frontdoor_wood_width = 1.5*m;

G4double lab68_glasswindow_height = 1.08*m;
G4double lab68_glasswindow_width = 0.57*m;

G4double lab68_reardoor_height = 2.3*m;
G4double lab68_reardoor_width = 0.91*m;

G4double lab68_frontdoor_x_coordinate = lab68_wall_x/2 -lab68_frontdoor_glass_width/2.0-1.2*m;
G4double lab68_reardoor_x_coordinate = -lab68_wall_x/2 +lab68_reardoor_width/2.0+2.6*m;

G4ThreeVector position_of_origin = {2.7*m, -2.45*m, 1.3*m}; //with repect to the inner upper left corner of the room(Doug's corner)

G4ThreeVector xyposition_of_origin = {2.7*m, -2.45*m, 0};

G4double Phantom_Radius=0.4*m;
G4double Phantom_Height=2.0*m;
G4double Phantom_Size=0.25*m;


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

  G4RotationMatrix* turnAlong190 = new G4RotationMatrix;
  turnAlong190->rotateZ(190*deg);

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

  //Lab_68 include ceiling
  G4VSolid* Main_2_S = new G4Box("Main_2_solid", lab68_wall_x/2.0, lab68_wall_y/2.0 , lab68_wall_z/2.0);
  //G4VSolid* hole_2_S = new G4Box("hole_2_solid", (lab68_wall_x-2*lab68_wall_thickness)/2.0, (lab68_wall_y-2*lab68_wall_thickness)/2.0, (lab68_wall_z-2*lab68_wall_thickness)/2.0);
  G4VSolid* hole_2_S = new G4Box("hole_2_solid", (lab68_wall_x-lab68_wall_thickness)/2.0, (lab68_wall_y-2*lab68_wall_thickness)/2.0, (lab68_wall_z-2*lab68_wall_thickness)/2.0);
  G4SubtractionSolid* LabFloorWall_solid_S= new G4SubtractionSolid("LabFloorWall_solid", Main_2_S, hole_2_S, NO_ROT, G4ThreeVector(0.,0., 0.));
  //G4SubtractionSolid* Main_2a_S= new G4SubtractionSolid("Main_2a_solid", Main_2_S, hole_2_S, NO_ROT, G4ThreeVector(0.,0., 0.));
  G4LogicalVolume* LabFloorWall_solid_LV = new G4LogicalVolume(LabFloorWall_solid_S, Concrete, "LabFloorWall_solid");
  LabFloorWall_solid_PV = new G4PVPlacement(turnAlong, G4ThreeVector{lab68_wall_x/2.0,-lab68_wall_y/2.0,lab68_wall_z/2.0}-position_of_origin, LabFloorWall_solid_LV, "LabFloorWall", vacuum_solid_LV, false, 0, fCheckOverlaps);
  LabFloorWall_solid_LV->SetVisAttributes(G4VisAttributes(G4Colour::Grey()));



  ////Lab donot include ceiling
  //G4VSolid* LabFloorExtended_solid_S=  new G4Box("LabFloorExtended_solid", 15.0*m, 15.0*m , 15.0*m);
  //G4LogicalVolume* LabFloorExtended_solid_LV = new G4LogicalVolume(LabFloorExtended_solid_S, Soil, "LabFloorExtended_solid");
  //LabFloorExtended_solid_PV = new G4PVPlacement(turnAlong, G4ThreeVector{lab68_wall_x/2.0,-lab68_wall_y/2.0,lab68_wall_z/2.0}-position_of_origin-G4ThreeVector(0., 0., lab68_wall_z/2.0+15.0*m), LabFloorExtended_solid_LV, "LabFloor_extended", vacuum_solid_LV, false, 0, fCheckOverlaps);
  //LabFloorExtended_solid_LV->SetVisAttributes(G4VisAttributes(G4Colour::Brown()));
  ////LabFloorExtended_solid_LV->SetVisAttributes(G4VisAttributes::Invisible);

  G4VSolid* frontglassdoor_S = new G4Box("frontglassdoor_solid", lab68_frontdoor_glass_width/2.0, lab68_wall_thickness/2.0, lab68_frontdoor_glass_height/2.0);
  //G4SubtractionSolid* Main_2b_S= new G4SubtractionSolid("Main_2b_solid", Main_2a_S, frontdoor_S, NO_ROT, G4ThreeVector(lab68_frontdoor_x_coordinate ,-lab68_wall_y/2.0, -(lab68_wall_z-2*lab68_wall_thickness-lab68_frontdoor_glass_height)/2.0));
  G4LogicalVolume* frontglassdoor_LV = new G4LogicalVolume(frontglassdoor_S, Quartz, "frontglassdoor");
  frontglassdoor_PV = new G4PVPlacement(NO_ROT, G4ThreeVector(lab68_frontdoor_x_coordinate ,-(lab68_wall_y-lab68_wall_thickness)/2.0, -(lab68_wall_z-2*lab68_wall_thickness-lab68_frontdoor_glass_height)/2.0), frontglassdoor_LV, "Front Glass Door", LabFloorWall_solid_LV, false, 0, fCheckOverlaps);
  frontglassdoor_LV->SetVisAttributes(G4VisAttributes(G4Colour::Cyan()));


  G4VSolid* frontdoor_S = new G4Box("frontdoor_solid", lab68_frontdoor_wood_width/2.0, lab68_wall_thickness/2.0, lab68_frontdoor_wood_height/2.0);
  //G4SubtractionSolid* Main_2b_S= new G4SubtractionSolid("Main_2b_solid", Main_2a_S, frontdoor_S, NO_ROT, G4ThreeVector(lab68_frontdoor_x_coordinate ,-lab68_wall_y/2.0, -(lab68_wall_z-2*lab68_wall_thickness-lab68_frontdoor_glass_height)/2.0));
  G4LogicalVolume* frontdoor_LV = new G4LogicalVolume(frontdoor_S, Wood, "frontdoor");
  frontdoor_PV = new G4PVPlacement(NO_ROT, G4ThreeVector(lab68_frontdoor_glass_width/2.0-lab68_frontdoor_wood_width/2.0 ,0, -lab68_frontdoor_glass_height/2.0+lab68_frontdoor_wood_height/2.0), frontdoor_LV, "Front Wood Door", frontglassdoor_LV, false, 0, fCheckOverlaps);
  frontdoor_LV->SetVisAttributes(G4VisAttributes(G4Colour::Brown()));


  G4VSolid* glasswindow_S = new G4Box("glasswindow_solid", lab68_glasswindow_width/2.0, lab68_wall_thickness/2.0, lab68_glasswindow_height/2.0);
  //G4SubtractionSolid* Main_2b_S= new G4SubtractionSolid("Main_2b_solid", Main_2a_S, frontdoor_S, NO_ROT, G4ThreeVector(lab68_frontdoor_x_coordinate ,-lab68_wall_y/2.0, -(lab68_wall_z-2*lab68_wall_thickness-lab68_frontdoor_glass_height)/2.0));
  G4LogicalVolume* glasswindow_LV = new G4LogicalVolume(glasswindow_S, Quartz, "glasswindow");
  glasswindow_PV = new G4PVPlacement(NO_ROT, G4ThreeVector(lab68_frontdoor_wood_width/2.0-0.16*m-lab68_glasswindow_width/2.0, 0, -lab68_frontdoor_wood_height/2.0+1.11*m+lab68_glasswindow_height/2.0), glasswindow_LV, "Front Glass Window", frontdoor_LV, false, 0, fCheckOverlaps);
  glasswindow_LV->SetVisAttributes(G4VisAttributes(G4Colour::Cyan()));


  G4VSolid* reardoor_S = new G4Box("reardoor_solid", lab68_reardoor_width/2.0, lab68_wall_thickness/2.0, lab68_reardoor_height/2.0);
  //G4SubtractionSolid* LabFloorWall_solid_S= new G4SubtractionSolid("LabFloorWall_solid", Main_2b_S, reardoor_S, NO_ROT, G4ThreeVector(lab68_reardoor_x_coordinate ,lab68_wall_y/2.0, -(lab68_wall_z-2*lab68_wall_thickness-lab68_reardoor_height)/2.0));
  G4LogicalVolume* reardoor_LV = new G4LogicalVolume(reardoor_S, Wood, "reardoor");
  reardoor_PV = new G4PVPlacement(NO_ROT, G4ThreeVector(lab68_reardoor_x_coordinate ,(lab68_wall_y-lab68_wall_thickness)/2.0, -(lab68_wall_z-2*lab68_wall_thickness-lab68_reardoor_height)/2.0), reardoor_LV, "Rear Wooden Door", LabFloorWall_solid_LV, false, 0, fCheckOverlaps);
  reardoor_LV->SetVisAttributes(G4VisAttributes(G4Colour::Brown()));

  //Lab_69 include ceiling
  G4ThreeVector position_lab66 = G4ThreeVector{-lab68_wall_x/2.0,-lab68_wall_y/2.0,lab68_wall_z/2.0}-position_of_origin
                                  +G4ThreeVector{0.0,TMath::Sin(TMath::DegToRad() * 10)*lab68_wall_x,0.0};
  lab70FloorWall_solid_PV = new G4PVPlacement(turnAlong, position_lab66, LabFloorWall_solid_LV, "lab66FloorWall", vacuum_solid_LV, false, 0, fCheckOverlaps);

  //Pump_chase_y
  //Lab_67 include ceiling
  G4ThreeVector position_lab70 = G4ThreeVector{3*lab68_wall_x/2.0,-lab68_wall_y/2.0,lab68_wall_z/2.0}-position_of_origin
                                  +G4ThreeVector{0.0,-TMath::Sin(TMath::DegToRad() * 10)*lab68_wall_x,0.0};
  lab70FloorWall_solid_PV = new G4PVPlacement(turnAlong, position_lab70, LabFloorWall_solid_LV, "lab70FloorWall", vacuum_solid_LV, false, 0, fCheckOverlaps);


  //Lab_30 include ceiling
  G4ThreeVector position_lab30 = G4ThreeVector{lab68_wall_x/2.0,lab68_wall_y/2.0+Pump_chase_y,lab68_wall_z/2.0}-position_of_origin
                                  +G4ThreeVector{TMath::Sin(TMath::DegToRad() * 10)*(lab68_wall_y+Pump_chase_y),0.0,0.0};
  lab30FloorWall_solid_PV = new G4PVPlacement(turnAlong190, position_lab30, LabFloorWall_solid_LV, "lab30FloorWall", vacuum_solid_LV, false, 0, fCheckOverlaps);

  //Lab_30 include ceiling
  G4ThreeVector position_lab28 = G4ThreeVector{3*lab68_wall_x/2.0,lab68_wall_y/2.0+Pump_chase_y,lab68_wall_z/2.0}-position_of_origin
                                  +G4ThreeVector{TMath::Sin(TMath::DegToRad() * 10)*(lab68_wall_y+Pump_chase_y),-TMath::Sin(TMath::DegToRad() * 10)*lab68_wall_x,0.0};
  lab28FloorWall_solid_PV = new G4PVPlacement(turnAlong190, position_lab28, LabFloorWall_solid_LV, "lab28FloorWall", vacuum_solid_LV, false, 0, fCheckOverlaps);

  //Lab_32 include ceiling
  G4ThreeVector position_lab32 = G4ThreeVector{-lab68_wall_x/2.0,lab68_wall_y/2.0+Pump_chase_y,lab68_wall_z/2.0}-position_of_origin
                                  +G4ThreeVector{TMath::Sin(TMath::DegToRad() * 10)*(lab68_wall_y+Pump_chase_y),TMath::Sin(TMath::DegToRad() * 10)*lab68_wall_x,0.0};
  Lab32FloorWall_solid_PV = new G4PVPlacement(turnAlong190, position_lab32, LabFloorWall_solid_LV, "Lab32FloorWall", vacuum_solid_LV, false, 0, fCheckOverlaps);



  //G4LogicalVolume* LabFloorWall_solid_LV = new G4LogicalVolume(LabFloorWall_solid_S, Concrete, "LabFloorWall_solid");
  //LabFloorWall_solid_PV = new G4PVPlacement(turnAlong, G4ThreeVector{lab68_wall_x/2.0,-lab68_wall_y/2.0,lab68_wall_z/2.0}-position_of_origin, LabFloorWall_solid_LV, "OutSpacer", vacuum_solid_LV, false, 0, fCheckOverlaps);
  //LabFloorWall_solid_LV->SetVisAttributes(G4VisAttributes(G4Colour::Grey()));
  //LabFloorWall_solid_LV->SetVisAttributes(G4VisAttributes::Invisible);



  //out_neutron box  delta/2.0
  G4VSolid* Main_3_S = new G4Box("Main_3_solid", (1.02*Water_x+delta)/2.0, (1.02*Water_y+delta)/2.0 , (1.02*Water_z+delta)/2.0);
  G4VSolid* hole_3_S = new G4Box("hole_3_solid", 1.02*Water_x/2.0, 1.02*Water_y/2.0, 1.02*Water_z/2.0);
  G4SubtractionSolid* TestSurface_solid_S= new G4SubtractionSolid("TestSurface_solid", Main_3_S, hole_3_S, NO_ROT, G4ThreeVector(0.,0., 0.));
  G4LogicalVolume* TestSurface_solid_LV = new G4LogicalVolume(TestSurface_solid_S, Vacuum, "TestSurface_solid");
  TestSurface_solid_PV = new G4PVPlacement(turnAlongZ, G4ThreeVector(0., fFilterCellSpacing+Water_y/2.0, 0), TestSurface_solid_LV, "TestSurface", vacuum_solid_LV, false, 0, fCheckOverlaps);
  TestSurface_solid_LV->SetVisAttributes(G4VisAttributes(G4Colour::Grey()));
  //TestSurface_solid_LV->SetVisAttributes(G4VisAttributes::Invisible);


  //layer of air at the bottom of the filter
  G4VSolid* Main_S = new G4Box("Main_solid", Water_x/2.0, Water_y/2.0, Water_z/2.0);
  G4VSolid* hole_S = new G4Box("hole_solid", Poly_a/2.0 , Poly_a/2.0, NeutronFilter_length/2.0);
  G4SubtractionSolid* boratedwater_S= new G4SubtractionSolid("boratedwater", Main_S, hole_S, turnAlongX, G4ThreeVector(0., Water_y/2.0-Water_rear_side-NeutronFilter_length/2.0, 0.));
  //G4LogicalVolume* boratedwater_LV = new G4LogicalVolume(boratedwater_S, BoraxWater, "boratedwater");
  G4LogicalVolume* boratedwater_LV = new G4LogicalVolume(boratedwater_S, BoraxBoricAcidBuffer, "boratedwater");
  //G4LogicalVolume* boratedwater_LV = new G4LogicalVolume(boratedwater_S, BoratedPoly, "boratedwater");
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
  Test_FRONTSIDE_PV = new G4PVPlacement(turnAlongZ, G4ThreeVector(0., fFilterCellSpacing-delta/2.0,0), Test_FRONTSIDE_LV, "Test_FRONTSIDE",vacuum_solid_LV, false, 0, fCheckOverlaps);
  Test_FRONTSIDE_LV->SetVisAttributes(G4VisAttributes(G4Colour::Yellow()));

  //At the center
  //G4VSolid* Test_CENTERPOINT_S = new G4Box("Test_CENTERPOINT_solid", Water_x/2.0, delta/2.0, Water_z/2.0);
  //G4LogicalVolume *Test_CENTERPOINT_LV = new G4LogicalVolume(Test_CENTERPOINT_S, Vacuum,"Test_CENTERPOINT" );
  //Test_CENTERPOINT_PV = new G4PVPlacement(turnAlongZ, G4ThreeVector(0., 0., 0.), Test_CENTERPOINT_LV, "Test_CENTERPOINT", vacuum_solid_LV, false, 0, fCheckOverlaps);
  //Test_CENTERPOINT_LV->SetVisAttributes(G4VisAttributes(G4Colour::Yellow()));

  G4ThreeVector Origin_DT=G4ThreeVector(0., fFilterCellSpacing+NeutronFilter_length/2.0, 0.)
                  +G4ThreeVector(0., NeutronFilter_length/2.0-(fMultiplierLeadHeightRear+fMultiplierLeadHeightFront)/2.0,0.);

  G4double assymetric_factor_positivey=Water_rear_side+fMultiplierLeadHeightRear;
  G4double assymetric_factor_negativey=NeutronFilter_length-fMultiplierLeadHeightRear;

  G4ThreeVector Phantom_Placement=  Origin_DT + G4ThreeVector(ftestx, ftesty, 0);
  //G4ThreeVector Phantom_Placement=  Origin_DT + G4ThreeVector(ftestx+Water_x/2.0, 0, ftestz);
  ////G4ThreeVector Phantom_Placement_2= Origin_DT + G4ThreeVector(-ftestx, -ftesty, 0);
  //G4ThreeVector Phantom_Placement_2= Origin_DT + G4ThreeVector(-ftestx-Water_x/2.0, 0, ftestz);
  ////G4ThreeVector Phantom_Placement_3= Origin_DT + G4ThreeVector(ftestx, -ftesty, 0);
  //G4ThreeVector Phantom_Placement_3= Origin_DT + G4ThreeVector(0, -ftesty-assymetric_factor_negativey, ftestz);
  ////G4ThreeVector Phantom_Placement_4= Origin_DT + G4ThreeVector(-ftestx, ftesty, 0);
  //G4ThreeVector Phantom_Placement_4= Origin_DT + G4ThreeVector(0, ftesty+assymetric_factor_positivey, ftestz);
  ////G4ThreeVector Phantom_Placement_5= G4ThreeVector(lab68_wall_x/2.0,-lab68_wall_y/2.0, 0) - xyposition_of_origin;
  //G4ThreeVector Phantom_Placement_5= Origin_DT + G4ThreeVector(ftestx+Water_x/2.0, ftesty+assymetric_factor_positivey, ftestz);
  //G4ThreeVector Phantom_Placement_6= Origin_DT + G4ThreeVector(-ftestx-Water_x/2.0, ftesty+assymetric_factor_positivey, ftestz);
  //G4ThreeVector Phantom_Placement_7= Origin_DT + G4ThreeVector(ftestx+Water_x/2.0, -ftesty-assymetric_factor_negativey, ftestz);
  //G4ThreeVector Phantom_Placement_8= Origin_DT + G4ThreeVector(-ftestx-Water_x/2.0, -ftesty-assymetric_factor_negativey, ftestz);

  // Poly need to change
  //G4VSolid* Phantom_S = new G4Tubs("Phantom", zeroRadius, Phantom_Radius/2.0, Phantom_Height/2.0, startAngle, spanningAngle);
  //G4VSolid* Phantom_S = new G4Box("Phantom",Phantom_Size/2.0, Phantom_Size/2.0, Phantom_Size/2.0-delta);
  G4VSolid* Phantom_S = new G4Sphere("Phantom", zeroRadius,Phantom_Size, startAngle, spanningAngle, startAngle, spanningAngle);
  //G4LogicalVolume* Phantom_LV = new G4LogicalVolume(Phantom_S, Air, "Phantom");
  //G4LogicalVolume* Phantom_LV = new G4LogicalVolume(Phantom_S, Soft_Tissue, "Phantom");
  G4LogicalVolume* Phantom_LV = new G4LogicalVolume(Phantom_S, Vacuum, "Phantom");
  Phantom_PV = new G4PVPlacement(NO_ROT, Phantom_Placement, Phantom_LV, "Phantom", vacuum_solid_LV, false, 0, fCheckOverlaps);

  //Phantom2_PV = new G4PVPlacement(NO_ROT, Phantom_Placement_2, Phantom_LV, "Phantom2", vacuum_solid_LV, false, 0, fCheckOverlaps);

  //Phantom3_PV = new G4PVPlacement(NO_ROT, Phantom_Placement_3, Phantom_LV, "Phantom3", vacuum_solid_LV, false, 0, fCheckOverlaps);

  //Phantom4_PV = new G4PVPlacement(NO_ROT, Phantom_Placement_4, Phantom_LV, "Phantom4", vacuum_solid_LV, false, 0, fCheckOverlaps);

  ////Phantom5_PV = new G4PVPlacement(NO_ROT, Phantom_Placement_5, Phantom_LV, "Phantom5", vacuum_solid_LV, false, 0, fCheckOverlaps);

  //Phantom6_PV = new G4PVPlacement(NO_ROT, Phantom_Placement_6, Phantom_LV, "Phantom6", vacuum_solid_LV, false, 0, fCheckOverlaps);

  //Phantom7_PV = new G4PVPlacement(NO_ROT, Phantom_Placement_7, Phantom_LV, "Phantom7", vacuum_solid_LV, false, 0, fCheckOverlaps);

  //Phantom8_PV = new G4PVPlacement(NO_ROT, Phantom_Placement_8, Phantom_LV, "Phantom8", vacuum_solid_LV, false, 0, fCheckOverlaps);



  Phantom_LV->SetVisAttributes(G4VisAttributes(G4Colour::Green()));


  //Lab donot include ceiling
  G4VSolid* LabFloorExtended_solid_S=  new G4Box("LabFloorExtended_solid", 25.0*m, 25.0*m , 15.0*m);
  //G4SubtractionSolid* Main_2a_S= new G4SubtractionSolid("Main_2a_solid", Main_2_S, hole_2_S, NO_ROT, G4ThreeVector(0.,0., 0.));
  G4LogicalVolume* LabFloorExtended_solid_LV = new G4LogicalVolume(LabFloorExtended_solid_S, Soil, "LabFloorExtended_solid");
  //G4LogicalVolume* LabFloorExtended_solid_LV = new G4LogicalVolume(LabFloorExtended_solid_S, Vacuum, "LabFloorExtended_solid");
  LabFloorExtended_solid_PV = new G4PVPlacement(turnAlong, G4ThreeVector{lab68_wall_x/2.0,-lab68_wall_y/2.0,lab68_wall_z/2.0}-position_of_origin-G4ThreeVector(0., 0., lab68_wall_z/2.0+15.0*m), LabFloorExtended_solid_LV, "LabFloor_extended", vacuum_solid_LV, false, 0, fCheckOverlaps);
  LabFloorExtended_solid_LV->SetVisAttributes(G4VisAttributes(G4Colour::Brown()));
  //LabFloorExtended_solid_LV->SetVisAttributes(G4VisAttributes::Invisible);




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

void IronFilterDetectorConstruction::SetTestX(G4double ival)
{
  ftestx = ival;
}

void IronFilterDetectorConstruction::SetTestY(G4double ival)
{
  ftesty = ival;
}

void IronFilterDetectorConstruction::SetTestZ(G4double ival)
{
  ftestz = ival;
}
