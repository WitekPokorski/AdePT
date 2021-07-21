// SPDX-FileCopyrightText: 2021 CERN
// SPDX-License-Identifier: Apache-2.0

#ifndef ADEPT_TRACK_CUH
#define ADEPT_TRACK_CUH

#include <AdePT/MParray.h>
#include <CopCore/SystemOfUnits.h>
#include <CopCore/Ranluxpp.h>

#include <VecGeom/base/Vector3D.h>
#include <VecGeom/navigation/NavStateIndex.h>

// A data structure to represent a particle track. The particle type is implicit
// by the queue and not stored in memory.
struct Track {
  RanluxppDouble rngState;
  double energy;
  double numIALeft[3];

  vecgeom::Vector3D<double> pos;
  vecgeom::Vector3D<double> dir;
  vecgeom::NavStateIndex currentState;
  vecgeom::NavStateIndex nextState;

  __host__ __device__ double Uniform() { return rngState.Rndm(); }

  __host__ __device__ void SwapStates()
  {
    auto state         = this->currentState;
    this->currentState = this->nextState;
    this->nextState    = state;
  }

  __host__ __device__ void InitAsSecondary(const Track &parent)
  {
    // The caller is responsible to branch a new RNG state and to set the energy.
    this->numIALeft[0] = -1.0;
    this->numIALeft[1] = -1.0;
    this->numIALeft[2] = -1.0;

    // A secondary inherits the position of its parent; the caller is responsible
    // to update the directions.
    this->pos          = parent.pos;
    this->currentState = parent.currentState;
    this->nextState    = parent.nextState;
  }
};
#endif
