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

#include "Par04EMShowerMessenger.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4FastHit.hh"
#include "Randomize.hh"
#include "G4FastSimHitMaker.hh"

#include "RunShower.h"

#include "Par04EMShowerModel.hh"

#include <VecGeom/base/Config.h>
#include <VecGeom/management/GeoManager.h>

enum MaterialCutCouples {
  WorldMC = 0,
  GapMC,
  AbsorberMC,
};

Par04EMShowerModel::Par04EMShowerModel(G4String aModelName, G4Region* aEnvelope)
  : G4VFastSimulationModel(aModelName, aEnvelope)
  , fMessenger(new Par04EMShowerMessenger(this))
  , fHitMaker(new G4FastSimHitMaker)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04EMShowerModel::Par04EMShowerModel(G4String aModelName)
  : G4VFastSimulationModel(aModelName)
  , fMessenger(new Par04EMShowerMessenger(this))
  , fHitMaker(new G4FastSimHitMaker)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04EMShowerModel::~Par04EMShowerModel() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool Par04EMShowerModel::IsApplicable(
  const G4ParticleDefinition& aParticleType)
{
  return &aParticleType == G4Electron::ElectronDefinition() ||
         &aParticleType == G4Positron::PositronDefinition() ||
         &aParticleType == G4Gamma::GammaDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool Par04EMShowerModel::ModelTrigger(const G4FastTrack& aFastTrack)
{
  /*
  // Check energy
  if(aFastTrack.GetPrimaryTrack()->GetKineticEnergy() < 1 * GeV)
  {
    return false;
  }
  // Check length of detector
  // Calculate depth of the detector along shower axis to verify if shower
  // will fit inside. Required max shower depth is defined by fLongMaxDepth, and
  // can be changed with UI command `/Par04/fastSim/longitudinalProfile/maxDepth
  G4double X0 = aFastTrack.GetPrimaryTrack()->GetMaterial()->GetRadlen();
  auto particleDirection     = aFastTrack.GetPrimaryTrackLocalDirection();
  auto particlePosition      = aFastTrack.GetPrimaryTrackLocalPosition();
  G4double detectorDepthInMM = aFastTrack.GetEnvelopeSolid()->DistanceToOut(
    particlePosition, particleDirection);
  G4double detectorDepthInX0 = detectorDepthInMM / X0;
 */
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04EMShowerModel::DoIt(const G4FastTrack& aFastTrack,
                              G4FastStep& aFastStep)
{
  // Remove particle from further processing by G4
  aFastStep.KillPrimaryTrack();
  aFastStep.SetPrimaryTrackPathLength(0.0);
  G4double energy = aFastTrack.GetPrimaryTrack()->GetKineticEnergy();
  // No need to create any deposit, it will be handled by this model (and
  // G4FastSimHitMaker that will call the sensitive detector)
  aFastStep.SetTotalEnergyDeposited(0);
  auto particlePosition  = aFastTrack.GetPrimaryTrackLocalPosition();
  auto particleDirection = aFastTrack.GetPrimaryTrackLocalDirection();

std::cout << "Here I call AdePT to generate the shower" << std::endl;

const vecgeom::VPlacedVolume *world = vecgeom::GeoManager::Instance().GetWorld();

  // Place particles between the world boundary and the calorimeter.
  double startX = particlePosition.z();
//  double chargedTrackLength[NumVolumes];
//  double energyDeposit[NumVolumes];
  ScoringPerVolume scoringPerVolume;
//  scoringPerVolume.chargedTrackLength = chargedTrackLength;
//  scoringPerVolume.energyDeposit      = energyDeposit;
  GlobalScoring globalScoring;

int particles = 1;
energy *= copcore::units::GeV;
int batch = -1;

// I need to pass the particle from Geant4 to AdePT and simulate the shower
Shower(world, particles, energy, batch, startX, MCIndex, &scoringPerVolume, NumVolumes, &globalScoring);

// Create energy deposit in the detector
// This will call appropriate sensitive detector class
//      fHitMaker->make(G4FastHit(position, energy / fNbOfHits), aFastTrack);
//      generatedHits++;
    
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04EMShowerModel::Print() const
{
  G4cout << "Par04EMShowerModel: " << G4endl;
}

void Par04EMShowerModel::Initialize()
{

// Map VecGeom volume IDs to Geant4 material-cuts couples.
 
  // Fill world and calorimeter.
  MCIndex[0] = MCIndex[1] = WorldMC;
  for (int i = 2; i < NumVolumes; i += (1 + NbOfAbsorbers)) {
    MCIndex[i]     = WorldMC;
    MCIndex[i + 1] = GapMC;
    MCIndex[i + 2] = AbsorberMC;
  }

  CreateVecGeomWorld();  
} 
