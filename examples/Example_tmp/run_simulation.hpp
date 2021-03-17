// SPDX-FileCopyrightText: 2020 CERN
// SPDX-License-Identifier: Apache-2.0

#pragma once

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
  COPCORE_CALLABLE_DECLARE(processesFunc, create_processes);

  std::cout << "Executing simulation on " << copcore::BackendName(backend) << "\n";
  
  bool success = true;

  // Create a stream to work with. On the CPU backend, this will be equivalent with: int stream = 0;
  Stream_t stream;
  StreamStruct::CreateStream(stream);

  // call the kernel to create the processes to be run on the device

  copcore::Allocator<process*, backend> procAlloc;  
  process **processes = procAlloc.allocate(2);

  copcore::Allocator<process_list*, backend> proclistAlloc;  
  process_list **proclist = proclistAlloc.allocate(1); 

  *proclist = new process_list(processes, 2);

 cudaDeviceSynchronize();
 
  Launcher_t create_processes_launch(stream);
  create_processes_launch.Run(processesFunc, 1, {0,0}/*, proclist*/, processes);
  create_processes_launch.WaitStream();

  std::cout << "proclist " << (*proclist) << std::endl;
  std::cout << "created processes: " << (*proclist)->list_size << std::endl;

  return 1;
}

