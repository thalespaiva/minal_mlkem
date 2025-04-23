set(SOURCES kem.c indcpa.c polyvec.c poly.c fq.S shuffle.S ntt.S invntt.S
            basemul.S consts.c rejsample.c cbd.c verify.c)
set(SOURCESKECCAK ${SOURCES} fips202.c fips202x4.c symmetric-shake.c
                  keccak4x/KeccakP-1600-times4-SIMD256.c)

if(MSVC)
  add_compile_options(/nologo /O2 /W4 /wd4146 /wd4244)
else()
  add_compile_options(-Wall -Wextra -Wpedantic -Wmissing-prototypes -Wredundant-decls)
  add_compile_options(-Wshadow -Wpointer-arith -mavx2 -mbmi2 -mpopcnt -maes)
  add_compile_options(-march=native -mtune=native -O3 -fomit-frame-pointer)
endif()

set(kyber_sets 512 768 1024)
set(kyber_ks 2 3 4)
set(code_dims 1 2 4 8)

list(LENGTH kyber_sets len)
math(EXPR len "${len} - 1")

foreach(idx RANGE ${len})
  list(GET kyber_sets ${idx} set)
  list(GET kyber_ks ${idx} k)
  foreach(code_dim ${code_dims})
    message("A ${k} B ${set} C ${code_dim}")
    add_executable(test_kyber${set}_${code_dim}d ${SOURCESKECCAK} randombytes.c test_kyber.c)
    target_compile_definitions(test_kyber${set}_${code_dim}d PUBLIC KYBER_K=${k} KYBER_N_ENCODING_DIMENSIONS=${code_dim})

    add_executable(test_speed${set}_${code_dim}d ${SOURCESKECCAK} cpucycles.c speed_print.c test_speed.c kex.c)
    target_compile_definitions(test_speed${set}_${code_dim}d PUBLIC KYBER_K=${k} KYBER_N_ENCODING_DIMENSIONS=${code_dim})

    add_executable(test_kyber${set}_${code_dim}d_null_eta2 ${SOURCESKECCAK} randombytes.c test_kyber.c)
    target_compile_definitions(test_kyber${set}_${code_dim}d_null_eta2 PUBLIC KYBER_K=${k}
                               KYBER_N_ENCODING_DIMENSIONS=${code_dim} USE_NULL_ETA2=1)

    add_executable(test_speed${set}_${code_dim}d_null_eta2 ${SOURCESKECCAK} cpucycles.c speed_print.c test_speed.c kex.c)
    target_compile_definitions(test_speed${set}_${code_dim}d_null_eta2 PUBLIC KYBER_K=${k}
                               KYBER_N_ENCODING_DIMENSIONS=${code_dim} USE_NULL_ETA2=1)
  endforeach()
endforeach()

