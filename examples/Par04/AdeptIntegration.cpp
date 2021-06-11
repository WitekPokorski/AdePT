// SPDX-FileCopyrightText: 2021 CERN
// SPDX-License-Identifier: Apache-2.0

#include "AdeptIntegration.h"
#include "AdeptIntegration.cuh"

#include "VecGeom/management/GeoManager.h"

void AdeptIntegration::AddTrack(int pdg, double energy, double x, double y, double z, double dirx, double diry, double dirz)
{
  fBuffer.toDevice.emplace_back(pdg, energy, x, y, z, dirx, diry, dirz);
  if (pdg == 11)
    fBuffer.nelectrons++;
  else if (pdg == -11)
    fBuffer.npositrons++;
  else if (pdg == 22)
    fBuffer.ngammas++;
}

void AdeptIntegration::Initialize()
{
  if (fInit) return;
  assert(fMaxBatch > 0 && "AdeptImtegration::Initialize - Maximum batch size not set.");

  assert(vecgeom::GeoManager::Instance().IsClosed() && "VecGeom geometry not closed.");
  const vecgeom::cxx::VPlacedVolume *world = vecgeom::GeoManager::Instance().GetWorld();

  AdeptIntegration::InitializeGPU(world, fMaxBatch);
  fInit = true;
}

void AdeptIntegration::Cleanup()
{
  AdeptIntegration::FreeGPU();
}

void AdeptIntegration::Shower(int event)
{
  AdeptIntegration::ShowerGPU(event, fBuffer);
}
