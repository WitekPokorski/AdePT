// SPDX-FileCopyrightText: 2020 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <curand.h>
#include <curand_kernel.h>

#include <AdePT/BlockData.h>
#include <AdePT/MParray.h>
#include <AdePT/Atomic.h>

#include "process.h"
#include "process_list.h"
#include "pair_production.h"
#include "energy_loss.h"

#include "track.h"

// some simple scoring
struct Scoring {
  adept::Atomic_t<int> secondaries;
  adept::Atomic_t<float> totalEnergyLoss;

  VECCORE_ATT_HOST_DEVICE
  Scoring() {}

  VECCORE_ATT_HOST_DEVICE
  static Scoring *MakeInstanceAt(void *addr)
  {
    Scoring *obj = new (addr) Scoring();
    return obj;
  }
};

// portable kernel functions have to reside in a backend-dependent inline namespace to avoid symbol duplication
// when building executables/libraries running same functions on different backends
inline namespace COPCORE_IMPL {

// kernel select processes based on interaction lenght and put particles in the appropriate queues
__device__ void DefinePhysicalStepLength(int i, adept::BlockData<track> *block, process_list* proclist, adept::MParray **queues, curandState_t *states)
{
//  int n = block->GetNused() + block->GetNholes();

//  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i += blockDim.x * gridDim.x) {
    
    // skip particles that are already dead
    if ((*block)[i].status == dead) return;

    (proclist)->GetPhysicsInteractionLength(i, block, states); // return value (if step limited by physics or geometry) not used for the moment
    // now, I know which process wins, so I add the particle to the appropriate queue
    queues[(*block)[i].current_process]->push_back(i);
//  }
}

COPCORE_CALLABLE_FUNC(DefinePhysicalStepLength)

// kernel to call Along Step function for particles in the queues
__device__ void CallAlongStepProcess(int i, adept::BlockData<track> *block, process* proc, adept::MParray *queue, 
                                        Scoring *scor, curandState_t *states)
{
  int particle_index;
  // get particles index from the queue
  particle_index = (*queue)[i];
  // and call the process for it
  proc->GenerateInteraction(particle_index, block, states);

  // a simple version of scoring 
  scor->totalEnergyLoss.fetch_add((*block)[particle_index].energy_loss);
  scor->secondaries.fetch_add((*block)[particle_index].number_of_secondaries);

  // if particles returns with 'dead' status, release the element from the block
  if ((*block)[particle_index].status == dead) block->ReleaseElement(particle_index);
}

COPCORE_CALLABLE_FUNC(CallAlongStepProcess)

// kernel function to initialize the random states
__device__ void init(int, curandState_t *states)
{
  /* we have to initialize the state */
  curand_init(0, 0, 0, states);
}

COPCORE_CALLABLE_FUNC(init)

// kernel to create the processes and process list
__device__ void create_processes(int, process_list **proclist, process **processes)
{
  // instantiate the existing processes
  *(processes) = new energy_loss();
  *(processes+1) = new pair_production();

  // add them to process_list (process manager)
  *proclist = new process_list(processes, 2);

//printf("processes %p \n", (*proclist));
//printf("processes %d \n", (*proclist)->list_size);

}

COPCORE_CALLABLE_FUNC(create_processes)

} // End namespace COPCORE_IMPL
