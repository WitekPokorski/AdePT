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
#include "Par04DetectorConstruction.hh"
#include "Par04DetectorMessenger.hh"
#include "Par04SensitiveDetector.hh"
#include "Par04EMShowerModel.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4VisAttributes.hh"
#include "G4RunManager.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

#include "G4SDManager.hh"

#include "G4UnitsTable.hh"

#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"
#include <G4ProductionCutsTable.hh>

#include <VecGeom/base/Config.h>
#include <VecGeom/base/Transformation3D.h>
#include <VecGeom/management/GeoManager.h>
#include <VecGeom/volumes/PlacedVolume.h>
#include <VecGeom/volumes/UnplacedBox.h>

#include "Scoring.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04DetectorConstruction::Par04DetectorConstruction()
  : G4VUserDetectorConstruction()
{
  fDetectorMessenger = new Par04DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04DetectorConstruction::~Par04DetectorConstruction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* Par04DetectorConstruction::Construct()
{
  constexpr double CalorSizeYZ       = 40 * cm;
  constexpr int NbOfLayers           = 50;
  constexpr int NbOfAbsorbers        = 2;
  constexpr double GapThickness      = 2.3 * mm;
  constexpr double AbsorberThickness = 5.7 * mm;

  constexpr double LayerThickness = GapThickness + AbsorberThickness;
  constexpr double CalorThickness = NbOfLayers * LayerThickness;

  constexpr double WorldSizeX  = 1.2 * CalorThickness;
  constexpr double WorldSizeYZ = 1.2 * CalorSizeYZ;

  //--------- Material definition ---------
  const char *WorldMaterial    = "G4_Galactic";
  const char *GapMaterial      = "G4_Pb";
  const char *AbsorberMaterial = "G4_lAr";

  G4Material *default_mat  = G4NistManager::Instance()->FindOrBuildMaterial(WorldMaterial);
  G4Material *gap_mat      = G4NistManager::Instance()->FindOrBuildMaterial(GapMaterial);
  G4Material *absorber_mat = G4NistManager::Instance()->FindOrBuildMaterial(AbsorberMaterial);

  auto SolidWorld = new G4Box("World", WorldSizeX / 2., WorldSizeYZ / 2., WorldSizeYZ / 2.);
  auto LogicWorld = new G4LogicalVolume(SolidWorld, default_mat, "World");
  auto PhysiWorld = new G4PVPlacement(0, G4ThreeVector(), LogicWorld, "World", 0, false, 0);

  auto SolidCalor = new G4Box("Calorimeter", CalorThickness / 2., CalorSizeYZ / 2., CalorSizeYZ / 2.);
  auto LogicCalor = new G4LogicalVolume(SolidCalor, default_mat, "Calorimeter");
  auto PhysiCalor = new G4PVPlacement(0, G4ThreeVector(), LogicCalor, "Calorimeter", LogicWorld, false, 0);

  //
  // Layers
  //

  auto SolidLayer = new G4Box("Layer", LayerThickness / 2, CalorSizeYZ / 2, CalorSizeYZ / 2);
  fLogicLayer     = new G4LogicalVolume(SolidLayer, default_mat, "Layer");
  auto PhysiLayer = new G4PVReplica("Layer", fLogicLayer, LogicCalor, kXAxis, NbOfLayers, LayerThickness);

  //
  // Absorbers
  //

  G4double xfront = -0.5 * LayerThickness;
  auto SolidGap = new G4Box("Gap", GapThickness / 2, CalorSizeYZ / 2, CalorSizeYZ / 2);
  fLogicGap     = new G4LogicalVolume(SolidGap, gap_mat, "Gap_Pb");
  G4double xcenter = xfront + 0.5 * GapThickness;
  xfront += GapThickness;
  auto PhysiGap = new G4PVPlacement(0, G4ThreeVector(xcenter, 0., 0.), fLogicGap, "Gap_Pb", fLogicLayer, false, 0);

  auto SolidAbsorber = new G4Box("Gap", AbsorberThickness / 2, CalorSizeYZ / 2, CalorSizeYZ / 2);
  fLogicAbsorber     = new G4LogicalVolume(SolidAbsorber, absorber_mat, "Absorber_LAr");
  xcenter = xfront + 0.5 * AbsorberThickness;
  xfront += AbsorberThickness;
  auto PhysiAbsorber = new G4PVPlacement(0, G4ThreeVector(xcenter, 0., 0.), fLogicAbsorber, "Absorber_LAr", fLogicLayer, false, 0);

  // Region for fast simulation
  auto detectorRegion = new G4Region("DetectorRegion");
  detectorRegion->AddRootLogicalVolume(LogicCalor);
  detectorRegion->UsedInMassGeometry(true);

  G4ProductionCuts *productionCuts = new G4ProductionCuts();
  productionCuts->SetProductionCut(ProductionCut);
  //
  detectorRegion->SetProductionCuts(productionCuts);
  //
  
  G4ProductionCutsTable *theCoupleTable = G4ProductionCutsTable::GetProductionCutsTable();
  theCoupleTable->UpdateCoupleTable(PhysiWorld);

  Print();
  CreateVecGeomWorld();
  return PhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04DetectorConstruction::ConstructSDandField()
{
  if (fMagFieldVector.mag() > 0.0) {
    // Apply a global uniform magnetic field along the Z axis.
    // Notice that only if the magnetic field is not zero, the Geant4
    // transportion in field gets activated.
    auto uniformMagField     = new G4UniformMagField(fMagFieldVector);
    G4FieldManager *fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    fieldMgr->SetDetectorField(uniformMagField);
    fieldMgr->CreateChordFinder(uniformMagField);
    G4cout << G4endl << " *** SETTING MAGNETIC FIELD : fieldValue = " << fMagFieldVector / kilogauss
           << " [kilogauss] *** " << G4endl << G4endl;

  } else {
    G4cout << G4endl << " *** NO MAGNETIC FIELD SET  *** " << G4endl << G4endl;
  }

  Par04SensitiveDetector* caloSD_gap = new Par04SensitiveDetector(
    "sensitiveDetectorGap", fNbOfLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(caloSD_gap);
  SetSensitiveDetector(fLogicGap, caloSD_gap);

  Par04SensitiveDetector* caloSD_absorber = new Par04SensitiveDetector(
    "sensitiveDetectorAbsorber", fNbOfLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(caloSD_absorber);
  SetSensitiveDetector(fLogicAbsorber, caloSD_absorber);

  auto detectorRegion =
    G4RegionStore::GetInstance()->GetRegion("DetectorRegion");
  Par04EMShowerModel* showermodel = new Par04EMShowerModel("model", detectorRegion);
  showermodel->Initialize();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04DetectorConstruction::Print() const
{
  G4cout << "\n------------------------------------------------------"
         << "\n--- Number of layers:\t" << fNbOfLayers;
  G4cout << "-----------------------------------------------------" << G4endl;
}

void Par04DetectorConstruction::CreateVecGeomWorld()
{
  auto worldSolid = new vecgeom::UnplacedBox(0.5 * WorldSizeX, 0.5 * WorldSizeYZ, 0.5 * WorldSizeYZ);
  auto worldLogic = new vecgeom::LogicalVolume("World", worldSolid);
  vecgeom::VPlacedVolume *worldPlaced = worldLogic->Place();

  //
  // Calorimeter
  //
  auto calorSolid = new vecgeom::UnplacedBox(0.5 * CalorThickness, 0.5 * CalorSizeYZ, 0.5 * CalorSizeYZ);
  auto calorLogic = new vecgeom::LogicalVolume("Calorimeter", calorSolid);
  vecgeom::Transformation3D origin;
  worldLogic->PlaceDaughter(calorLogic, &origin);

  //
  // Layers
  //
  auto layerSolid = new vecgeom::UnplacedBox(0.5 * LayerThickness, 0.5 * CalorSizeYZ, 0.5 * CalorSizeYZ);

  //
  // Absorbers
  //
  auto gapSolid = new vecgeom::UnplacedBox(0.5 * GapThickness, 0.5 * CalorSizeYZ, 0.5 * CalorSizeYZ);
  auto gapLogic = new vecgeom::LogicalVolume("Gap", gapSolid);
  vecgeom::Transformation3D gapPlacement(-0.5 * LayerThickness + 0.5 * GapThickness, 0, 0);

  auto absorberSolid = new vecgeom::UnplacedBox(0.5 * AbsorberThickness, 0.5 * CalorSizeYZ, 0.5 * CalorSizeYZ);
  auto absorberLogic = new vecgeom::LogicalVolume("Absorber", absorberSolid);
  vecgeom::Transformation3D absorberPlacement(0.5 * LayerThickness - 0.5 * AbsorberThickness, 0, 0);

  // Create a new LogicalVolume per layer, we need unique IDs for scoring.
  double xCenter = -0.5 * CalorThickness + 0.5 * LayerThickness;
  for (int i = 0; i < NbOfLayers; i++) {
    auto layerLogic = new vecgeom::LogicalVolume("Layer", layerSolid);
    vecgeom::Transformation3D placement(xCenter, 0, 0);
    calorLogic->PlaceDaughter(layerLogic, &placement);

    layerLogic->PlaceDaughter(gapLogic, &gapPlacement);
    layerLogic->PlaceDaughter(absorberLogic, &absorberPlacement);

    xCenter += LayerThickness;
  }

  vecgeom::GeoManager::Instance().SetWorldAndClose(worldPlaced);
}
