// SPDX-FileCopyrightText: 2021 CERN
// SPDX-License-Identifier: Apache-2.0

#include "AdeptIntegration.h"
#include "AdeptIntegration.cuh"

#include <G4RunManager.hh>
#include <G4Threading.hh>

#include "VecGeom/management/GeoManager.h"

void AdeptIntegration::AddTrack(int pdg, double energy, double x, double y, double z, double dirx, double diry,
                                double dirz)
{
  fBuffer.toDevice.emplace_back(pdg, energy, x, y, z, dirx, diry, dirz);
  if (pdg == 11)
    fBuffer.nelectrons++;
  else if (pdg == -11)
    fBuffer.npositrons++;
  else if (pdg == 22)
    fBuffer.ngammas++;
}

void AdeptIntegration::Initialize(bool common_data)
{
  if (fInit) return;
  assert(fMaxBatch > 0 && "AdeptIntegration::Initialize - Maximum batch size not set.");

  if (common_data) {
    std::cout << "=== AdeptIntegration: initializing geometry and physics\n";
    // Initialize geometry on device
    assert(vecgeom::GeoManager::Instance().IsClosed() && "VecGeom geometry not closed.");
    const vecgeom::cxx::VPlacedVolume *world = vecgeom::GeoManager::Instance().GetWorld();
    if (!InitializeGeometry(world))
      throw std::runtime_error("AdeptIntegration::Initialize cannot initialize geometry on GPU");

    // Initialize G4HepEm
    if (!InitializePhysics()) throw std::runtime_error("AdeptIntegration::Initialize cannot initialize physics on GPU");

    // Initialize user segmentation on device
    UserMCIndex::GetInstance().InitializeOnGPU();
    return;
  }

  std::cout << "=== AdeptIntegration: initializing transport engine for thread: " << G4Threading::G4GetThreadId()
            << std::endl;

  // Initialize user scoring data
  fUserData.InitializeOnGPU();

  // Initialize the transport engine for the current thread
  InitializeGPU();

  fInit = true;
}

void AdeptIntegration::Cleanup()
{
  if (!fInit) return;
  // if (!fCleanup.test_and_set()) {
  //  std::cout << "=== Cleaning up GPU ...\n";
  AdeptIntegration::FreeGPU();
  //}
}

void AdeptIntegration::Shower(int event)
{
  AdeptIntegration::ShowerGPU(event, fBuffer);
  fBuffer.Clear();
}
