// SPDX-FileCopyrightText: 2020 CERN
// SPDX-License-Identifier: Apache-2.0

#include <curand.h>
#include <curand_kernel.h>
#include <iostream>
#include <stdio.h>

#include <AdePT/BlockData.h>

#include "process.h"
#include "process_list.h"
#include "pair_production.h"
#include "energy_loss.h"

#include "track.h"

#include <AdePT/MParray.h>

__global__ void DefinePhysicalStepLength(adept::BlockData<track> *block, process_list** proclist, adept::MParray **queues, curandState_t *states)
{
  int n = block->GetNused() + block->GetNholes();

  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i += blockDim.x * gridDim.x) {
    
    if ((*block)[i].status == dead) continue;

    (*proclist)->GetPhysicsInteractionLength(i, block, states);
    // now, I know which process wins, so I add the particle to the appropriate queue
    queues[(*block)[i].current_process]->push_back(i);
  }
}


__global__ void CallAlongStepProcesses(adept::BlockData<track> *block, process_list** proclist, adept::MParray **queues, curandState_t *states)
{
  int particle_index;

  for (int process_id=0 ; process_id < (*proclist)->list_size; process_id++) 
    {
      int queue_size = queues[process_id]->size();

      for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < queue_size; i += blockDim.x * gridDim.x) 
        {
          particle_index = (*(queues[process_id]))[i];

          ((*proclist)->list)[process_id]->GenerateInteraction(particle_index, block, states);

          if ((*block)[particle_index].status == dead) block->ReleaseElement(particle_index);
        }
    }
}

/* this GPU kernel function is used to initialize the random states */
__global__ void init(curandState_t *states)
{
  /* we have to initialize the state */
  curand_init(0, 0, 0, states);
}

// kernel to create the processes and process list
__global__ void create_processes(process_list **proclist, process **processes)
{
  *(processes) = new energy_loss();
  *(processes+1) = new pair_production();

  *proclist = new process_list(processes, 2);
}

//

int main()
{

  curandState_t *state;
  cudaMalloc((void **)&state, sizeof(curandState_t));
  init<<<1, 1>>>(state);
  cudaDeviceSynchronize();

  // kernel to create the processes to be run on the device

  process_list **proclist;
  process **processes;
  cudaMalloc((void**)&proclist, sizeof(process_list*));
  cudaMalloc((void**)&processes, 2*sizeof(process*));

  create_processes<<<1,1>>>(proclist, processes);
  cudaDeviceSynchronize();

  // Capacity of the different containers
  constexpr int capacity = 1 << 20;

  constexpr int numberOfProcesses = 3;
  char *buffer[numberOfProcesses];

  adept::MParray **queues = nullptr;
  cudaMallocManaged(&queues, numberOfProcesses * sizeof(adept::MParray *));

  size_t buffersize = adept::MParray::SizeOfInstance(capacity);

  for (int i = 0; i < numberOfProcesses; i++) {
    buffer[i] = nullptr;
    cudaMallocManaged(&buffer[i], buffersize);

    queues[i] = adept::MParray::MakeInstanceAt(capacity, buffer[i]);
  }

  
  // Allocate a block of tracks with capacity larger than the total number of spawned threads
  size_t blocksize = adept::BlockData<track>::SizeOfInstance(capacity);
  char *buffer2    = nullptr;
  cudaMallocManaged(&buffer2, blocksize);
  auto block = adept::BlockData<track>::MakeInstanceAt(capacity, buffer2);

  // initializing one track in the block
  auto track    = block->NextElement();
  track->energy = 100.0f;
  track->index = 1;

  // initializing second track in the block
  //  auto track2    = block->NextElement();
  //  track2->energy = 30.0f;
  //  track2->index = 2;

  //
  constexpr dim3 nthreads(32);
  constexpr dim3 maxBlocks(10);
  dim3 numBlocks;

  while (block->GetNused()>0) 
  {
    numBlocks.x = (block->GetNused() + block->GetNholes() + nthreads.x - 1) / nthreads.x;
    numBlocks.x = std::min(numBlocks.x, maxBlocks.x);

    // call the kernel to do check the step lenght and select process
    DefinePhysicalStepLength<<<numBlocks, nthreads>>>(block, proclist, queues, state);
    
    // call the kernel for Along Step Processes
    CallAlongStepProcesses<<<numBlocks, nthreads>>>(block, proclist, queues, state);

    cudaDeviceSynchronize();

    for (int i = 0; i < numberOfProcesses; i++) queues[i]->clear();
  
    cudaDeviceSynchronize();

    std::cout << "Nused: " << block->GetNused() << std::endl;

  }
}
