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

using Queue_t = adept::mpmc_bounded_queue<int>;



__global__ void DefinePhysicalStepLength(int n, adept::BlockData<track> *block, process_list** proclist, Queue_t *queues[], curandState_t *states)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i += blockDim.x * gridDim.x) {
    // check if you are not outside the used block
    if (i > block->GetNused() + block->GetNholes()) return;

    if ((*block)[i].status == dead) continue;

    (*proclist)->GetPhysicsInteractionLength(i, block, states);

    // now, I know which process wins, so I add the particle to the appropriate queue
    queues[(*block)[i].current_process]->enqueue(i);
  }
}


__global__ void CallAlongStepProcesses(int n, adept::BlockData<track> *block, process_list** proclist, Queue_t *queues[], curandState_t *states)
{
  int particle_index;

  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i += blockDim.x * gridDim.x) {
    for (int process_id=0; process_id < (*proclist)->list_size; process_id++)
    {
      if (!queues[process_id]->dequeue(particle_index)) continue;

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

  using Queue_t = adept::mpmc_bounded_queue<int>;

  constexpr int numberOfProcesses = 3;
  char *buffer[numberOfProcesses];

  Queue_t **queues = nullptr;
  cudaMallocManaged(&queues, numberOfProcesses * sizeof(Queue_t *));

  size_t buffersize = Queue_t::SizeOfInstance(capacity);

  for (int i = 0; i < numberOfProcesses; i++) {
    buffer[i] = nullptr;
    cudaMallocManaged(&buffer[i], buffersize);

    queues[i] = Queue_t::MakeInstanceAt(capacity, buffer[i]);
  }

  
  // Allocate a block of tracks with capacity larger than the total number of spawned threads
  // Note that if we want to allocate several consecutive block in a buffer, we have to use
  // Block_t::SizeOfAlignAware rather than SizeOfInstance to get the space needed per block

  using Block_t    = adept::BlockData<track>;
  size_t blocksize = Block_t::SizeOfInstance(capacity);
  char *buffer2    = nullptr;
  cudaMallocManaged(&buffer2, blocksize);
  auto block = Block_t::MakeInstanceAt(capacity, buffer2);

  // initializing one track in the block
  auto track    = block->NextElement();
  track->energy = 100.0f;
  track->index = 1;

  // initializing second track in the block
//   auto track2    = block->NextElement();
//  track2->energy = 30.0f;
//  track2->index = 2;
  //
  constexpr dim3 nthreads(32);
  constexpr dim3 maxBlocks(10);
  dim3 numBlocks;

  while (block->GetNused()>0) 
  {

    int n = block->GetNused() + block->GetNholes();

    numBlocks.x = (block->GetNused() + block->GetNholes() + nthreads.x - 1) / nthreads.x;
    numBlocks.x = std::min(numBlocks.x, maxBlocks.x);

    // call the kernel to do check the step lenght and select process
    DefinePhysicalStepLength<<<numBlocks, nthreads>>>(n, block, proclist, queues, state);
    
    // call the kernel for Along Step Processes
    CallAlongStepProcesses<<<numBlocks, nthreads>>>(n, block, proclist, queues, state);

    cudaDeviceSynchronize();

    std::cout << "Nused: " << block->GetNused() << std::endl;

  }
}
