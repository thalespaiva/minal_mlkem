set(KYBER_SRCS kex.c kem.c indcpa.c polyvec.c poly.c ntt.c cbd.c reduce.c verify.c)
set(KYBER_FIPS202_SRCS ${KYBER_SRCS} symmetric-shake.c)
set(KYBER_NINETIES_SRCS ${KYBER_SRCS} symmetric-aes.c)
set(FIPS202_SRCS fips202.c)
set(AES256CTR_SRCS aes256ctr.c)
set(SHA2_SRCS sha256.c sha512.c)
set(TEST_KYBER_SRCS test_kyber.c randombytes.c)
set(TEST_KEX_SRCS test_kex.c randombytes.c)
set(TEST_VECTORS_SRCS test_vectors.c)
set(TEST_SPEED_SRCS test_speed.c speed_print.c cpucycles.c randombytes.c)

if(MSVC)
  add_compile_options(/nologo /O2 /W4 /wd4146 /wd4244)
else()
  add_compile_options(-Wall -Wextra -Wpedantic -Werror)
  add_compile_options(-Wmissing-prototypes -Wredundant-decls -Wshadow -Wpointer-arith)
  add_compile_options(-O3 -fomit-frame-pointer)
  add_compile_options(-march=native -mtune=native -O3 -fwrapv)
endif()

add_library(fips202_ref ${FIPS202_SRCS})

set(kyber_sets 512 768 1024)
set(kyber_ks 2 3 4)
set(code_dims 1 2 4 8)

list(LENGTH kyber_sets len)
math(EXPR len "${len} - 1")

foreach(idx RANGE ${len})
  list(GET kyber_sets ${idx} set)
  list(GET kyber_ks ${idx} k)
  foreach(code_dim ${code_dims})
    # message("A ${k} B ${set} C ${code_dim}")
    add_executable(test_kyber_ref${set}_${code_dim}d ${TEST_KYBER_SRCS} ${KYBER_FIPS202_SRCS})
    target_compile_definitions(test_kyber_ref${set}_${code_dim}d PUBLIC KYBER_K=${k} KYBER_N_ENCODING_DIMENSIONS=${code_dim})
    target_link_libraries(test_kyber_ref${set}_${code_dim}d fips202_ref)

    add_executable(test_speed_ref${set}_${code_dim}d ${TEST_SPEED_SRCS} ${KYBER_FIPS202_SRCS} ${FIPS202_SRCS})
    target_compile_definitions(test_speed_ref${set}_${code_dim}d PUBLIC KYBER_K=${k} KYBER_N_ENCODING_DIMENSIONS=${code_dim})
    target_link_libraries(test_speed_ref${set}_${code_dim}d fips202_ref)


    add_executable(test_kyber_ref${set}_${code_dim}d_null_eta2 ${TEST_KYBER_SRCS} ${KYBER_FIPS202_SRCS} ${FIPS202_SRCS})
    target_compile_definitions(test_kyber_ref${set}_${code_dim}d_null_eta2 PUBLIC KYBER_K=${k}
                               KYBER_N_ENCODING_DIMENSIONS=${code_dim} USE_NULL_ETA2=1)
    target_link_libraries(test_kyber_ref${set}_${code_dim}d_null_eta2 fips202_ref)

    add_executable(test_speed_ref${set}_${code_dim}d_null_eta2 ${TEST_SPEED_SRCS} ${KYBER_FIPS202_SRCS} ${FIPS202_SRCS})
    target_compile_definitions(test_speed_ref${set}_${code_dim}d_null_eta2 PUBLIC KYBER_K=${k}
                               KYBER_N_ENCODING_DIMENSIONS=${code_dim} USE_NULL_ETA2=1)
    target_link_libraries(test_speed_ref${set}_${code_dim}d_null_eta2 fips202_ref)

  endforeach()
endforeach()

