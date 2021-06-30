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
#ifndef PAR03DETECTORCONSTRUCTION_H
#define PAR03DETECTORCONSTRUCTION_H

#include "G4VUserDetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4ThreeVector.hh"

class Par04DetectorMessenger;
class G4LogicalVolume;

/**
 * @brief Detector construction.
 *
 * Creates a cylindrical detector, with cylinder axis along Z-axis. It is placed
 * in the world volume so that its bases are located at z=0 and z=Length.
 * Dimensions of the detector (Radius and Length) and material can be set using
 * the UI commands.
 * Readout geometry of the detector is created, and can be set by UI commands.
 * Cells are created along z-axis, azimuthal angle, and radius (cylindrical
 * segmentation).
 * Sensitive detector Par04SensitiveDetector is attached to the
 * cell volume.
 * Region for the detector is created as an envelope of the fast simulation.
 *
 */

class Par04DetectorConstruction : public G4VUserDetectorConstruction
{
 public:
  Par04DetectorConstruction();
  virtual ~Par04DetectorConstruction();

  virtual G4VPhysicalVolume* Construct() final;
  void CreateVecGeomWorld();
  virtual void ConstructSDandField() final;

  // Set number of readout cells along z-axis
  inline void SetNbOfLayers(G4int aNumber) { fNbOfLayers = aNumber; };
  // Get number of readout cells along z-axis
  inline G4int GetNbOfLayers() const { return fNbOfLayers; };
  // Set uniform magnetic field
  inline void SetMagField(const G4ThreeVector &fv) { fMagFieldVector = fv; }

  // Print detector information
  void Print() const;

 private:
  /// Messenger that allows to modify geometry
  Par04DetectorMessenger* fDetectorMessenger;
  /// Logical volume of replicated cell
  G4LogicalVolume* fLogicLayer = nullptr;
   /// Logical volume of gap
  G4LogicalVolume* fLogicGap = nullptr;
   /// Logical volume of absorber
  G4LogicalVolume* fLogicAbsorber = nullptr;
  /// Number of layers = slices along z axis
  G4int fNbOfLayers = 50;
  // field related members
  G4ThreeVector fMagFieldVector;

};

#endif /* PAR03DETECTORCONSTRUCTION_H */
