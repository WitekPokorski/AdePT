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


// portable kernel functions have to reside in a backend-dependent inline namespace to avoid symbol duplication
// when building executables/libraries running same functions on different backends
inline namespace COPCORE_IMPL {

// kernel to create the processes and process list
__device__ void create_processes(int, /*process_list **proclist, **/ process **processes)
{
  // instantiate the existing processes
  *(processes) = new energy_loss();
  *(processes+1) = new pair_production();

  // add them to process_list (process manager)
//  *proclist = new process_list(processes, 2);

}

COPCORE_CALLABLE_FUNC(create_processes)

} // End namespace COPCORE_IMPL
