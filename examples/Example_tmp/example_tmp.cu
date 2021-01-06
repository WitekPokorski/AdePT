// SPDX-FileCopyrightText: 2020 CERN
// SPDX-License-Identifier: Apache-2.0


#include "run_simulation.hpp"

///______________________________________________________________________________________
int main(void)
{
  int result;
  result = runSimulation<copcore::BackendType::CUDA>();
  return result;
}
