cmake_minimum_required(VERSION 3.10.2)
enable_language(C ASM)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -march=native")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Winit-self")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wno-array-parameter")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -funsigned-char -Wall -Wextra -Wpedantic")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wunused -Wcomment -Wchar-subscripts -Wuninitialized -Wshadow")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wwrite-strings -Wno-deprecated-declarations -Wno-unknown-pragmas -Wformat-security")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wcast-qual -Wunused-result -fPIC -fopenmp")

add_library(mp_matrix mp_matrix.c)
target_link_libraries(mp_matrix mpfr mpc m)

add_library(fft_mpc fft_mpc.c)
target_link_libraries(fft_mpc mp_matrix)

add_executable(fft_mpc_test fft_mpc.c)
target_compile_definitions(fft_mpc_test PRIVATE FFT_MPC_TEST)
target_link_libraries(fft_mpc_test mp_matrix)

add_library(dfr_2d_codes dfr_2d_codes.c)
target_link_libraries(dfr_2d_codes mp_matrix)

add_executable(minal2d_dfr_computation joint_distribution_2d.c)
target_compile_definitions(minal2d_dfr_computation PUBLIC COMPUTE_JOINT_DISTRIBUTIONS_AND_DFR)
target_link_libraries(minal2d_dfr_computation fft_mpc mp_matrix dfr_2d_codes)

add_library(prob_dist prob_dist.c)
target_link_libraries(prob_dist mpfr m)

add_executable(minal4d_dfr_computation minal4d_dfr_computation.c)
target_compile_definitions(minal4d_dfr_computation PUBLIC SAVE_PARTIAL_DISTRIBUTION_RESULTS=1)
target_link_libraries(minal4d_dfr_computation prob_dist fft_mpc mp_matrix)

add_executable(minal2d_dfr_computation_with_pk_compression11 joint_distribution_2d.c)
target_compile_definitions(minal2d_dfr_computation_with_pk_compression11 PUBLIC COMPUTE_JOINT_DISTRIBUTIONS_AND_DFR SIMULATE_FOR_PK_COMPRESSION_11BITS)
target_link_libraries(minal2d_dfr_computation_with_pk_compression11 fft_mpc mp_matrix dfr_2d_codes)

add_executable(minal4d_dfr_computation_with_pk_compression11 minal4d_dfr_computation.c)
target_compile_definitions(minal4d_dfr_computation_with_pk_compression11 PUBLIC SAVE_PARTIAL_DISTRIBUTION_RESULTS=1 SIMULATE_FOR_PK_COMPRESSION_11BITS)
target_link_libraries(minal4d_dfr_computation_with_pk_compression11 prob_dist fft_mpc mp_matrix)
