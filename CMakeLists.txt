cmake_minimum_required(VERSION 3.10.2)
enable_language(C ASM)

project(kyber_minal)

add_subdirectory(minal_src/dfr_computation)
add_subdirectory(minal_src/kem_implementation)
add_subdirectory(minal_src/minal_codes)
