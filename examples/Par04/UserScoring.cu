// SPDX-FileCopyrightText: 2021 CERN
// SPDX-License-Identifier: Apache-2.0

#include "UserScoring.h"
#include "UserScoring.cuh"

#include <CopCore/Global.h>
#include <CopCore/PhysicalConstants.h>

#include "Track.cuh" // not nice - we expose the track model here, interface of DepositEnergy to be changed

#include <iostream>
#include <iomanip>
#include <stdio.h>

__constant__ __device__ int *MCIndex = nullptr;
//__constant__ __device__ double BzFieldValue = 0;

// Deposit energy of particles still in flight.
__global__ void DepositEnergyKernel(Track *allTracks, const adept::MParray *queue, GlobalScoring *globalScoring,
                                    ScoringPerVolume *scoringPerVolume)
{
  int queueSize = queue->size();
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < queueSize; i += blockDim.x * gridDim.x) {
    const int slot      = (*queue)[i];
    Track &currentTrack = allTracks[slot];
    auto volume         = currentTrack.currentState.Top();
    if (volume == nullptr) {
      // The particle left the world, why wasn't it killed before?!
      continue;
    }
    int volumeID = volume->id();

    double energy = currentTrack.energy;
    atomicAdd(&globalScoring->energyDeposit, energy);
    atomicAdd(&scoringPerVolume->energyDeposit[volumeID], energy);
  }
}

void UserMCIndex::InitializeOnGPU()
{
  // Transfer MC indices.
  COPCORE_CUDA_CHECK(cudaMalloc(&fMCIndex_dev, sizeof(int) * NumVolumes));
  COPCORE_CUDA_CHECK(cudaMemcpy(fMCIndex_dev, fMCIndex, sizeof(int) * NumVolumes, cudaMemcpyHostToDevice));
  COPCORE_CUDA_CHECK(cudaMemcpyToSymbol(MCIndex, &fMCIndex_dev, sizeof(int *)));
}

void UserMCIndex::FreeGPU()
{
  COPCORE_CUDA_CHECK(cudaFree(fMCIndex_dev));
}

void UserData::InitializeOnGPU()
{
  // Allocate memory to score charged track length and energy deposit per volume.
  COPCORE_CUDA_CHECK(cudaMalloc(&fChargedTrackLength_dev, sizeof(double) * NumVolumes));
  COPCORE_CUDA_CHECK(cudaMemset(fChargedTrackLength_dev, 0, sizeof(double) * NumVolumes));
  COPCORE_CUDA_CHECK(cudaMalloc(&fEnergyDeposit_dev, sizeof(double) * NumVolumes));
  COPCORE_CUDA_CHECK(cudaMemset(fEnergyDeposit_dev, 0, sizeof(double) * NumVolumes));

  // Allocate and initialize scoring and statistics.
  COPCORE_CUDA_CHECK(cudaMalloc(&fGlobalScoring_dev, sizeof(GlobalScoring)));
  COPCORE_CUDA_CHECK(cudaMemset(fGlobalScoring_dev, 0, sizeof(GlobalScoring)));

  ScoringPerVolume scoringPerVolume_devPtrs;
  scoringPerVolume_devPtrs.chargedTrackLength = fChargedTrackLength_dev;
  scoringPerVolume_devPtrs.energyDeposit      = fEnergyDeposit_dev;
  COPCORE_CUDA_CHECK(cudaMalloc(&fScoringPerVolume_dev, sizeof(ScoringPerVolume)));
  COPCORE_CUDA_CHECK(
      cudaMemcpy(fScoringPerVolume_dev, &scoringPerVolume_devPtrs, sizeof(ScoringPerVolume), cudaMemcpyHostToDevice));
}

void UserData::FreeGPU()
{
  // Free resources.
  COPCORE_CUDA_CHECK(cudaFree(fChargedTrackLength_dev));
  COPCORE_CUDA_CHECK(cudaFree(fEnergyDeposit_dev));

  COPCORE_CUDA_CHECK(cudaFree(fGlobalScoring_dev));
  COPCORE_CUDA_CHECK(cudaFree(fScoringPerVolume_dev));
}

void UserData::CopyHitsToHost()
{
  // Transfer back scoring.
  COPCORE_CUDA_CHECK(cudaMemcpy(&fGlobalScoring, fGlobalScoring_dev, sizeof(GlobalScoring), cudaMemcpyDeviceToHost));

  // Transfer back the scoring per volume (charged track length and energy deposit).
  COPCORE_CUDA_CHECK(cudaMemcpy(fScoringPerVolume.chargedTrackLength, fChargedTrackLength_dev,
                                sizeof(double) * NumVolumes, cudaMemcpyDeviceToHost));
  COPCORE_CUDA_CHECK(cudaMemcpy(fScoringPerVolume.energyDeposit, fEnergyDeposit_dev, sizeof(double) * NumVolumes,
                                cudaMemcpyDeviceToHost));
  fGlobalScoring.Print();
}

void UserData::DepositEnergy(Track *allTracks, const adept::MParray *queue, int numBlocks, int numThreads,
                             cudaStream_t &stream)
{
  DepositEnergyKernel<<<numBlocks, numThreads, 0, stream>>>(allTracks, queue, fGlobalScoring_dev,
                                                            fScoringPerVolume_dev);
}
