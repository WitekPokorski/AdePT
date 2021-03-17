// SPDX-FileCopyrightText: 2020 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <curand.h>
#include <curand_kernel.h>
#include <iostream>
#include <iomanip>
#include <stdio.h>

#include "sim_kernels.h"

///______________________________________________________________________________________
template <copcore::BackendType backend>
int runSimulation()
{

  using Launcher_t     = copcore::Launcher<backend>;
  using StreamStruct   = copcore::StreamType<backend>;
  using Stream_t       = typename StreamStruct::value_type;

  // Boilerplate to get the pointers to the device functions to be used
  COPCORE_CALLABLE_DECLARE(stepLengthFunc, DefinePhysicalStepLength);
  COPCORE_CALLABLE_DECLARE(alongStepFunc, CallAlongStepProcess);
  COPCORE_CALLABLE_DECLARE(initFunc, init);
  COPCORE_CALLABLE_DECLARE(processesFunc, create_processes);

  std::cout << "Executing simulation on " << copcore::BackendName(backend) << "\n";
  //  const char *result[2] = {"FAILED", "OK"};

  //  bool testOK  = true;
  bool success = true;

  // Create a stream to work with. On the CPU backend, this will be equivalent with: int stream = 0;
  Stream_t stream;
  StreamStruct::CreateStream(stream);

  // call the kernel to initialize RND engine
  curandState_t *state;
  cudaMalloc((void **)&state, sizeof(curandState_t));

  Launcher_t init_launch(stream);

  init_launch.Run(initFunc, 1, {0,0}, state);
  init_launch.WaitStream();


  // call the kernel to create the processes to be run on the device

process_list **proclist_dev;
  process **processes;
  cudaMalloc((void **)&proclist_dev, sizeof(process_list *));
  cudaMalloc((void **)&processes, 2 * sizeof(process *));

  Launcher_t create_processes_launch(stream);
  create_processes_launch.Run(processesFunc, 1, {0,0}, proclist_dev, processes);
 create_processes_launch.WaitStream();

  cudaDeviceSynchronize();
  process_list *proclist;
  cudaMemcpy(&proclist, proclist_dev, sizeof(process_list *), cudaMemcpyDeviceToHost);
  //cudaFree(proclist_dev);

 cudaDeviceSynchronize();
/*
 copcore::Allocator<process*, backend> procAlloc;  
  process **processes = procAlloc.allocate(2);

 copcore::Allocator<process_list*, backend> proclistAlloc;  
process_list **proclist = proclistAlloc.allocate(1); 
*/

std::cout << "a" << std::endl;

std::cout << "created processes " << proclist->list << std::endl;

  // Capacity of the different containers
  constexpr int capacity = 1 << 20;

  // setting the number of existing processes
  constexpr int numberOfProcesses = 2;
  char *buffer1[numberOfProcesses];

  // reserving queues for each of the processes
  adept::MParray **queues = nullptr;
  cudaMallocManaged(&queues, numberOfProcesses * sizeof(adept::MParray *));
  size_t buffersize = adept::MParray::SizeOfInstance(capacity);

  for (int i = 0; i < numberOfProcesses; i++) {
    buffer1[i] = nullptr;
    cudaMallocManaged(&buffer1[i], buffersize);
    queues[i] = adept::MParray::MakeInstanceAt(capacity, buffer1[i]);
  }

  // Allocate the content of Scoring in a buffer
  char *buffer_scor = nullptr;
  cudaMallocManaged(&buffer_scor, sizeof(Scoring));
  Scoring *scor = Scoring::MakeInstanceAt(buffer_scor);
  // Initialize scoring
  scor->secondaries     = 0;
  scor->totalEnergyLoss = 0;
  
  // Allocate a block of tracks with capacity larger than the total number of spawned threads
  size_t blocksize = adept::BlockData<track>::SizeOfInstance(capacity);
  char *buffer2    = nullptr;
  cudaMallocManaged(&buffer2, blocksize);
  auto block = adept::BlockData<track>::MakeInstanceAt(capacity, buffer2);

  // initializing one track in the block
  auto track    = block->NextElement();
  track->energy = 100.0f;
  track->energy_loss = 0.0f;
  //  track->index = 1; // this is not use for the moment, but it should be a unique track index

  // initializing second track in the block
  auto track2    = block->NextElement();
  track2->energy = 30.0f;
  track2->energy_loss = 0.0f;
  //  track2->index = 2; // this is not use for the moment, but it should be a unique track index

  // simple version of scoring
  float* energy_deposition = nullptr;
  cudaMalloc((void **)&energy_deposition, sizeof(float));

  constexpr dim3 nthreads(32);
  constexpr dim3 maxBlocks(10);
  dim3 numBlocks;


cudaDeviceSynchronize();

  Launcher_t stepLengthFunc_launch(stream);
  Launcher_t alongStepFunc_launch(stream);

  while (block->GetNused()>0) 
  {
//    numBlocks.x = (block->GetNused() + block->GetNholes() + nthreads.x - 1) / nthreads.x;
//    numBlocks.x = std::min(numBlocks.x, maxBlocks.x);

    // call the kernel to do check the step lenght and select process
    stepLengthFunc_launch.Run(stepLengthFunc, block->GetNused() + block->GetNholes(), {0,0}, block, proclist, queues, state);
    stepLengthFunc_launch.WaitStream();

std::cout << "number of processes " << std::endl;
    // call the kernel for Along Step Processes for each process
    for (int process_id=0 ; process_id < (proclist)->list_size; process_id++)  {
      auto process = ((proclist)->list)[process_id];
      auto queue = queues[process_id];
      // Launch the kernel per process
      alongStepFunc_launch.Run(alongStepFunc, queue->size(), {0, 0}, block, process, queue, scor, state);
  }
  alongStepFunc_launch.WaitStream();

  // clear all the queues before next step
  for (int i = 0; i < numberOfProcesses; i++) queues[i]->clear();
  cudaDeviceSynchronize();

  std::cout << "Number of tracks in flight: " ;
  std::cout << std::setw(8) << block->GetNused() ;
  
  std::cout << " total energy depostion: " << std::setw(10) << scor->totalEnergyLoss.load() 
  << " total number of secondaries: " << scor->secondaries.load() << std::endl;
  }
  return 1;
}

