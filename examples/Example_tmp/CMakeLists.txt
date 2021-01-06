# SPDX-FileCopyrightText: 2020 CERN
# SPDX-License-Identifier: Apache-2.0

# Example tmp of particle processing with GPU
add_executable(example_tmp example_tmp.cu)
target_include_directories(example_tmp PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/base/inc>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/physics/processes/inc>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/tracking/inc>
  $<INSTALL_INTERFACE:base>
)
target_link_libraries(example_tmp PRIVATE CUDA::curand VecCore::VecCore CopCore::CopCore)
