// SPDX-FileCopyrightText: 2021 CERN
// SPDX-License-Identifier: Apache-2.0

#ifndef SCORING_H
#define SCORING_H

#include <CopCore/Global.h>
#include <CopCore/PhysicalConstants.h>

static const char *WorldMaterial    = "G4_Galactic";
static const char *GapMaterial      = "G4_Pb";
static const char *AbsorberMaterial = "G4_lAr";

enum MaterialCutCouples {
  WorldMC = 0,
  GapMC,
  AbsorberMC,
};

constexpr double ProductionCut = 0.7 * copcore::units::mm;

constexpr double CalorSizeYZ = 40 * copcore::units::cm;
constexpr int NbOfLayers     = 50;
constexpr int NbOfAbsorbers  = 2;
constexpr int NumVolumes     = 1 + 1 + NbOfLayers * (1 + NbOfAbsorbers);

constexpr double GapThickness      = 2.3 * copcore::units::mm;
constexpr double AbsorberThickness = 5.7 * copcore::units::mm;

constexpr double LayerThickness = GapThickness + AbsorberThickness;
constexpr double CalorThickness = NbOfLayers * LayerThickness;

constexpr double WorldSizeX  = 1.2 * CalorThickness;
constexpr double WorldSizeYZ = 1.2 * CalorSizeYZ;

// Data structures for scoring. The accessors must make sure to use atomic operations if needed.
struct GlobalScoring {
  double energyDeposit;
  // Not int to avoid overflows for more than 100,000 events; unsigned long long
  // is the only other data type available for atomicAdd().
  unsigned long long chargedSteps;
  unsigned long long neutralSteps;
  unsigned long long hits;
  unsigned long long numGammas;
  unsigned long long numElectrons;
  unsigned long long numPositrons;

  void Print()
  {
    printf("Global scoring: stpChg=%llu stpNeu=%llu hits=%llu numGam=%llu numEle=%llu numPos=%llu\n", chargedSteps,
           neutralSteps, hits, numGammas, numElectrons, numPositrons);
  }
};

struct ScoringPerVolume {
  double *energyDeposit;
  double *chargedTrackLength;
};

struct UserMCIndex {
  int fMCIndex[NumVolumes];
  int *fMCIndex_dev{nullptr};

  static UserMCIndex &GetInstance()
  {
    static UserMCIndex userMCIndex;
    return userMCIndex;
  }

  UserMCIndex()
  {
    // Map VecGeom volume IDs to Geant4 material-cuts couples.
    // Fill world and calorimeter.
    fMCIndex[0] = fMCIndex[1] = WorldMC;
    for (int i = 2; i < NumVolumes; i += (1 + NbOfAbsorbers)) {
      fMCIndex[i]     = WorldMC;
      fMCIndex[i + 1] = GapMC;
      fMCIndex[i + 2] = AbsorberMC;
    }
  }

  ~UserMCIndex() { FreeGPU(); }

  void InitializeOnGPU();
  void FreeGPU();
};

struct Track;
namespace adept {
class MParray;
}

struct UserData {
  double *fEnergyDeposit_dev{nullptr};
  double *fChargedTrackLength_dev{nullptr};

  double fChargedTrackLength[NumVolumes];
  double fEnergyDeposit[NumVolumes];
  ScoringPerVolume fScoringPerVolume;
  ScoringPerVolume *fScoringPerVolume_dev{nullptr};
  GlobalScoring fGlobalScoring;
  GlobalScoring *fGlobalScoring_dev{nullptr};

  UserData()
  {
    fScoringPerVolume.chargedTrackLength = fChargedTrackLength;
    fScoringPerVolume.energyDeposit      = fEnergyDeposit;
  }

  void InitializeOnGPU();
  void CopyHitsToHost();
  void FreeGPU();

  // We have to define here also the scoring kernel.
#ifdef COPCORE_CUDA_COMPILER
  void DepositEnergy(Track *allTracks, const adept::MParray *queue, int numBlocks, int numThreads,
                     cudaStream_t &stream);
#endif
};

#endif
