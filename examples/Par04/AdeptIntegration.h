// SPDX-FileCopyrightText: 2021 CERN
// SPDX-License-Identifier: Apache-2.0

/// The The Geant4 AdePT integration service. This provides the interfaces for:
/// - initializing the geometry and physics on the AdePT size
/// - filling the buffer with tracks to be transported on the GPU
/// - Calling the Shower method transporting a buffer on the GPU

#ifndef ADEPT_INTEGRATION_H
#define ADEPT_INTEGRATION_H

#include <vector>
#include <atomic>
#include <VecGeom/base/Config.h>
#ifdef VECGEOM_ENABLE_CUDA
#include <VecGeom/management/CudaManager.h> // forward declares vecgeom::cxx::VPlacedVolume
#endif

#include <G4HepEmData.hh>
#include <G4HepEmParameters.hh>

#include "UserScoring.h" // This should be removed from here!

class AdeptIntegration {
public:
  static constexpr int kMaxThreads = 256;

  /// @brief Track data exchanged between Geant4 and AdePT
  struct TrackData {
    double position[3];
    double direction[3];
    double energy;
    int pdg;

    TrackData(int pdg_id, double ene, double x, double y, double z, double dirx, double diry, double dirz)
        : position{x, y, z}, direction{dirx, diry, dirz}, energy{ene}, pdg{pdg_id}
    {
    }
  };

  /// @brief Buffer holding input tracks to be transported on GPU and output tracks to be
  /// re-injected in the Geant4 stack
  struct TrackBuffer {
    std::vector<TrackData> toDevice;   ///< Tracks to be transported on the device
    std::vector<TrackData> fromDevice; ///< Tracks coming from device to be transported on the CPU
    int nelectrons{0};
    int npositrons{0};
    int ngammas{0};

    void Clear()
    {
      toDevice.clear();
      fromDevice.clear();
      nelectrons = npositrons = ngammas = 0;
    }
  };

  struct G4HepEmState {
    G4HepEmData data;
    G4HepEmParameters parameters;
  };

  /// @brief Returns the service instance

private:
  bool fInit{false};                ///< Service initialized flag
  std::atomic_flag fCleanup{false}; ///< Make sure we cleanup only once
  int fNthreads{0};                 ///< Number of cpu threads
  int fMaxBatch{0};                 ///< Max batch size for allocating GPU memory
  TrackBuffer fBuffer;              ///< Vector of buffers of tracks to/from device (per thread)
  TrackData *toDevice_dev{nullptr}; ///< Track buffer on device
  void *fGPUstate{nullptr};         ///< CUDA state placeholder

  void InitializeUserSegmentation() { UserMCIndex::GetInstance().InitializeOnGPU(); }
  void InitializeUserData() { fUserData.InitializeOnGPU(); }
  bool InitializeGeometry(const vecgeom::cxx::VPlacedVolume *world);
  bool InitializePhysics();
  void InitializeGPU();
  void ShowerGPU(int event, TrackBuffer const &buffer);
  void FreeGPU();

public:
  UserData fUserData;                 ///< User data (to be removed)
  static G4HepEmState *g4hepem_state; ///< The HepEm state singleton

public:
  /// @brief Adds a track to the buffer
  void AddTrack(int pdg, double energy, double x, double y, double z, double dirx, double diry, double dirz);
  /// @brief Set maximum batch size
  void SetMaxBatch(int npart) { fMaxBatch = npart; }
  /// @brief Initialize service and copy geometry & physics data on device
  void Initialize(bool common_data = false);
  /// @brief Final cleanup
  void Cleanup();
  /// @brief Interface for transporting a buffer of tracks in AdePT.
  void Shower(int event);
};

#endif
