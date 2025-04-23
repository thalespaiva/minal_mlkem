# Tailorable codes for lattice-based KEMs with applications to compact ML-KEM instantiations

> Paiva, Thales B., Marcos A. Simplicio Jr, Syed Mahbub Hafiz, Bahattin Yildiz, Eduardo L. Cominetti, and Henrique S. Ogawa. "Tailorable codes for lattice-based KEMs with applications to compact ML-KEM instantiations." IACR Transactions on Cryptographic Hardware and Embedded Systems (2022): To appear.

The paper can also be downloaded from eprint at:
https://eprint.iacr.org/2024/1243.pdf

**Abstract.** This repository contains the code and data associated with our paper.
In this README, we present a guide to explore the repository, together with direct
instructions on how to reproduce all the results in the paper.
This repository will be made publicly available after the reviewing process is done.

**This file is also provided in HTML format as README.html, which is more convenient for reading.**

## Important notices

**1. External code:**
In this repository, we use external codes from two sources provided by Kyber's Team:

* [Kyber's reference implementation](https://github.com/pq-crystals/kyber):
    This was used as a basis for implementing our different encoding solution and
    also for comparing the performance between the different methods.

* [Kyber's reference scripts for security estimation](https://github.com/pq-crystals/security-estimates): We used their code to verify that it is safe to use more aggressive
ciphertext compression parameters (`du` and `dv`) with respect to the LWE instances protecting
the message.

**2. On the exact reproduction of data**: To exactly reproduce the results in this paper, you
have to use the exact number of threads we used to generate it using:
```
$ export OMP_NUM_THREADS=24
```
The reason is that, otherwise, we cannot have guarantees on the order of the summations and how they are split into
partial summations. This is a problem because the associative property is not valid for floating point
operations, even when with the arbitrary precision operations provided by MPFR and MPC.
For an illustration of the problem, see this example in Python (without arbitrary precision).
```Python
In [9]: (0.1 + 0.2) + 0.3 + 0.6
Out[9]: 1.2000000000000002

In [10]: 0.1 + 0.2 + (0.3 + 0.6)
Out[10]: 1.2
```
However, notice that we use 260 bits of precision for the operations, which should be more than enough
to accurately compute the decryption failure rates (DFRs) of interest (that typically are
above 2^{-180}). Therefore, we don't expect the (multithreaded) reproduction **with a number of threads
different than 24** to be absolutely exact, but the differences should be negligible.

**3. Running scripts.** All scripts should be called from the project's root.

**4. TLDR.** We provide a script `reproduce_me.sh` for full reproduction. One caveat: it creates ~1GB of files
and takes about 20 hours to run on a server with 2x Intel Xeon Gold 5118 CPU @ 2.30GHz
(24 cores, 48 threads), with 64GB of RAM, (using `export OMP_NUM_THREADS=24`).



# Table of contents

<!-- MarkdownTOC -->

- Overview of this package
- Reproducing the results
    - System setup and dependencies
    - The main script for reproducing all results
- Short guide to the code
    - Compilation
    - Source code
        - Implementation of 2D codes over Kyber's reference/AVX2 code
        - Implementation of 4D codes over Kyber's reference/AVX2 code
        - Computation of the joint noise distribution `delta_m[i, i + n/2]`
- Description of the results
    - Performance results in `results/performance`
    - The DFR results
        - Case 2D:
        - Case 4D:
- Reproducing the figures from the paper

<!-- /MarkdownTOC -->


# Overview of this package

* `README.md` or `README.html`: This file
* `reproduce_me.sh`: Shell script used to reproduce all our data
* `analysis_scripts/generate_graphs.py`: Python script to generate all plots used in the paper
* `analysis_scripts/minal.py`: Python implementation of Minal codes and auxiliary functions for their analysis
* `experiments/`: Our experiments' setup files and running scripts
* `results/`: The experiments' results
* `minal_src/dfr_computation/`: C code with the DFR analysis of Kyber using 2D and 4D codes
* `minal_src/kem_implementation/`: C code with the constant-time implementation of Kyber using 2D, 4D and 8D.
* `minal_src/minal_codes/`: C code with the constant-time the general decoding algorithm in dimensions from 2 to 14


# Reproducing the results

First we list the necessary dependencies, then we introduce the script for full reproduction.

## System setup and dependencies

**System setup**

* The main system used for developing experiments was
a Linux PC with an Intel i7-8700 CPU @ 3.20GHz (6 cores, 12 threads) and 32GB of RAM.

* For actually running the experiments, we used a Linux server with 2 Intel Xeon Gold 5118 CPU @ 2.30GHz
(24 cores, 48 threads), with 64GB of RAM.

All tools required for reproduction are somewhat standard in scientific computing
and should be available in the official repositories of major Linux distributions.

**C libraries**

* MPFR (for arbitrary precision real numbers)
* MPC (for arbitrary precision complex numbers)
* OpenMP
* OpenSSL (required by Kyber)

**Tools**

* CMake
* gcc or clang
* Python3 + `matplotlib` + `seaborn` + `pandas` + `numpy` + `scipy`
* Optional: `Pipenv` for easier installation of the same versions we used


**Versions**

* Linux PC
```
$ cat /proc/cpuinfo | grep "model name" | head -n 1
model name  : Intel(R) Core(TM) i7-8700 CPU @ 3.20GHz
$ gcc --version | head -n 1
gcc (GCC) 13.2.1 20240417
$ cmake --version | head -n 1
cmake version 3.29.2
$ python3 --version
Python 3.12.3
$ pipenv --version
pipenv, version 2023.12.1
```

* Linux Server
```
$ cat /proc/cpuinfo | grep "model name" | head -n 1
model name  : Intel(R) Xeon(R) Gold 5118 CPU @ 2.30GHz
$ gcc --version | head -n 1
gcc (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0
$ cmake --version | head -n 1
cmake version 3.14.3
$ python3 --version
Python 3.10.12
$ # Pipenv not used in the server
```

## The main script for reproducing all results

The results can be fully reproduced using the `reproduce_me.sh` script.
Notice, however, that this will take a while.
The results took approximately 20 hours to be generated on a high-end server
with 2x Intel Xeon Gold 5118 CPU @ 2.30GHz, using 24 threads of the 48 total threads.

More specifically, we used the following definition before running the scripts.
```
$ export OMP_NUM_THREADS=24
```

The main script should be run in the root of the project with a directory path as the argument `reproduced_results`.
To avoid overwriting, the directory `reproduced_results` should not yet be created when calling the script.
```
$ ./reproduce_me.sh reproduced_results
```

All programs are called with `tee` so the command above will write a lot to `stdout`,
which should give you some information in case a problem occurs.


# Short guide to the code

We now elaborate on how one interested in compiling and running the targets individually for their
own experimentation can compile and explore the code.

## Compilation

Compilation is easy using `cmake`, `make` and `gcc`.
Please make sure you have covered the dependencies mentioned in the previous section.

```
$ cmake -B build
$ cd build
$ make
```

This will take some time (about 3 minutes) and generate multiple targets.
Most of them are related to to the 1D/2D/4D implementation of Kyber (ref and AVX2).
The targets can be seen below.



**General implementation of the Minal decoding with timings**
```
build/minal_src/minal_codes/minal_codes
```

Running this program results in a simple CSV showing the number of dimensions and the decoding time per bit
for 2D to 14D codes.
**Remember that this is a general, non-optimized implementation.**
However, we implemented it following constant-time programming practices hoping that it can make it easier for
designers of lattice-based schemes to test our codes in their constructions.
```
$ ./build/minal_src/minal_codes/minal_codes
n_dimensions,decode_median_cycles_per_bit,decode_average_cycles_per_bit
2,134.00,174.50
3,108.67,109.33
4,128.75,131.25
5,223.20,232.60
6,419.17,375.50
7,314.57,386.00
8,427.62,457.00
9,643.22,655.11
10,1177.00,1185.40
11,2232.18,2243.82
12,4355.17,4391.83
13,8797.31,8965.23
14,17230.36,17457.43
```


**DFR computation of Minal 2D/4D codes when applied to Kyber**

The executables below use accept CSV files (as found in `experiments/setup/exp_kyber*d_dfr*.csv`).
The last two consider that the pk is compressed with 11 bits (instead of the 12 bits used by Kyber).
```
build/minal_src/dfr_computation/minal4d_dfr_computation
build/minal_src/dfr_computation/minal2d_dfr_computation
build/minal_src/dfr_computation/minal2d_dfr_computation_with_pk_compression11
build/minal_src/dfr_computation/minal4d_dfr_computation_with_pk_compression11
```

**Performance tests of Kyber (original and our 2D/4D implementation)**
```
build/minal_src/kem_implementation/ref/test_speed_ref1024_4d
build/minal_src/kem_implementation/ref/test_speed_ref512_2d
build/minal_src/kem_implementation/ref/test_speed_ref768_2d_null_eta2
build/minal_src/kem_implementation/ref/test_speed_ref768_1d
build/minal_src/kem_implementation/ref/test_speed_ref512_4d_null_eta2
build/minal_src/kem_implementation/ref/test_speed_ref512_1d
build/minal_src/kem_implementation/ref/test_speed_ref1024_1d_null_eta2
build/minal_src/kem_implementation/ref/test_speed_ref1024_2d_null_eta2
build/minal_src/kem_implementation/ref/test_speed_ref1024_1d
build/minal_src/kem_implementation/ref/test_speed_ref768_1d_null_eta2
build/minal_src/kem_implementation/ref/test_speed_ref768_4d
build/minal_src/kem_implementation/ref/test_speed_ref1024_4d_null_eta2
build/minal_src/kem_implementation/ref/test_speed_ref512_1d_null_eta2
build/minal_src/kem_implementation/ref/test_speed_ref512_4d
build/minal_src/kem_implementation/ref/test_speed_ref512_2d_null_eta2
build/minal_src/kem_implementation/ref/test_speed_ref768_2d
build/minal_src/kem_implementation/ref/test_speed_ref768_4d_null_eta2
build/minal_src/kem_implementation/ref/test_speed_ref1024_2d
build/minal_src/kem_implementation/avx2/test_speed768_2d
build/minal_src/kem_implementation/avx2/test_speed512_1d
build/minal_src/kem_implementation/avx2/test_speed1024_1d
build/minal_src/kem_implementation/avx2/test_speed512_2d
build/minal_src/kem_implementation/avx2/test_speed512_2d_null_eta2
build/minal_src/kem_implementation/avx2/test_speed1024_4d_null_eta2
build/minal_src/kem_implementation/avx2/test_speed512_1d_null_eta2
build/minal_src/kem_implementation/avx2/test_speed768_2d_null_eta2
build/minal_src/kem_implementation/avx2/test_speed768_4d
build/minal_src/kem_implementation/avx2/test_speed512_4d_null_eta2
build/minal_src/kem_implementation/avx2/test_speed1024_2d
build/minal_src/kem_implementation/avx2/test_speed512_4d
build/minal_src/kem_implementation/avx2/test_speed768_4d_null_eta2
build/minal_src/kem_implementation/avx2/test_speed768_1d
build/minal_src/kem_implementation/avx2/test_speed1024_2d_null_eta2
build/minal_src/kem_implementation/avx2/test_speed1024_4d
build/minal_src/kem_implementation/avx2/test_speed768_1d_null_eta2
build/minal_src/kem_implementation/avx2/test_speed1024_1d_null_eta2
```

**Sanity-check tests for Kyber (original and our 2D/4D implementation)**
```
build/minal_src/kem_implementation/ref/test_kyber_ref512_2d_null_eta2
build/minal_src/kem_implementation/ref/test_kyber_ref512_4d_null_eta2
build/minal_src/kem_implementation/ref/test_kyber_ref768_4d
build/minal_src/kem_implementation/ref/test_kyber_ref1024_2d_null_eta2
build/minal_src/kem_implementation/ref/test_kyber_ref768_2d
build/minal_src/kem_implementation/ref/test_kyber_ref768_4d_null_eta2
build/minal_src/kem_implementation/ref/test_kyber_ref768_1d_null_eta2
build/minal_src/kem_implementation/ref/test_kyber_ref1024_1d
build/minal_src/kem_implementation/ref/test_kyber_ref768_2d_null_eta2
build/minal_src/kem_implementation/ref/test_kyber_ref1024_4d_null_eta2
build/minal_src/kem_implementation/ref/test_kyber_ref1024_4d
build/minal_src/kem_implementation/ref/test_kyber_ref768_1d
build/minal_src/kem_implementation/ref/test_kyber_ref512_4d
build/minal_src/kem_implementation/ref/test_kyber_ref1024_2d
build/minal_src/kem_implementation/ref/test_kyber_ref1024_1d_null_eta2
build/minal_src/kem_implementation/ref/test_kyber_ref512_1d
build/minal_src/kem_implementation/ref/test_kyber_ref512_2d
build/minal_src/kem_implementation/ref/test_kyber_ref512_1d_null_eta2
build/minal_src/kem_implementation/avx2/test_kyber512_2d
build/minal_src/kem_implementation/avx2/test_kyber512_4d_null_eta2
build/minal_src/kem_implementation/avx2/test_kyber512_2d_null_eta2
build/minal_src/kem_implementation/avx2/test_kyber768_1d
build/minal_src/kem_implementation/avx2/test_kyber1024_2d_null_eta2
build/minal_src/kem_implementation/avx2/test_kyber512_1d_null_eta2
build/minal_src/kem_implementation/avx2/test_kyber1024_4d_null_eta2
build/minal_src/kem_implementation/avx2/test_kyber768_4d_null_eta2
build/minal_src/kem_implementation/avx2/test_kyber768_2d
build/minal_src/kem_implementation/avx2/test_kyber1024_2d
build/minal_src/kem_implementation/avx2/test_kyber512_1d
build/minal_src/kem_implementation/avx2/test_kyber1024_4d
build/minal_src/kem_implementation/avx2/test_kyber1024_1d_null_eta2
build/minal_src/kem_implementation/avx2/test_kyber768_1d_null_eta2
build/minal_src/kem_implementation/avx2/test_kyber768_4d
build/minal_src/kem_implementation/avx2/test_kyber512_4d
build/minal_src/kem_implementation/avx2/test_kyber1024_1d
build/minal_src/kem_implementation/avx2/test_kyber768_2d_null_eta2
```

## Source code

We now briefly discuss two important files of our source code that are good entry points
for exploration.

### Implementation of 2D codes over Kyber's reference/AVX2 code


We use conditional compilation to compile Kyber with our 2D/4D code implementation.

* **For REF implementation: the only additions we make are in file `minal_src/kem_implementation/ref/poly.c`.**
* **For the AVX2 implementation: the only additions we make are in file `minal_src/kem_implementation/avx2/poly.c`.**

First, we use `experiments/scripts/generate_decoding_lines.py` to generate the auxiliary decoding
lines for the optimal 2D codes found for the **original** parameters for Levels 1, 3 and 5.
We add them to the beginning of `poly.c` file.

```C
#if defined(KYBER_N_ENCODING_DIMENSIONS) && (KYBER_N_ENCODING_DIMENSIONS == 2)

///////////// GENERATED WITH experiments/scripts/generate_decoding_lines.py /////////////////////////////
# if KYBER_K == 2

# define CODE_ALPHA (KYBER_Q/2)
# define CODE_BETA 422

static int32_t L[5][3] = {
    {-211, 832, 1438953},
    {844, 3330, -2950309},
    {-832, 211, 736745},
    {-1, 1, -1243},
    {3330, 844, 139789},
};

# define ABOVE_L1(x, y) lower_than_mask(L[0][0]*x + L[0][2], L[0][1]*y)
# define ABOVE_L2(x, y) lower_than_mask(L[1][0]*x + L[1][2], L[1][1]*y)
# define ABOVE_L3(x, y) lower_than_mask(L[2][0]*x + L[2][2], L[2][1]*y)
# define ABOVE_L4(x, y) lower_than_mask(L[3][0]*x + L[3][2], L[3][1]*y)
# define ABOVE_L5(x, y) lower_than_mask(L[4][0]*x + L[4][2], L[4][1]*y)

# elif KYBER_K == 3

# define CODE_ALPHA (KYBER_Q/2)
# define CODE_BETA 422

static int32_t L[5][3] = {
    {-211, 832, 1438953},
    {844, 3330, -2950309},
    {-832, 211, 736745},
    {-1, 1, -1243},
    {3330, 844, 139789},
};

# define ABOVE_L1(x, y) lower_than_mask(L[0][0]*x + L[0][2], L[0][1]*y)
# define ABOVE_L2(x, y) lower_than_mask(L[1][0]*x + L[1][2], L[1][1]*y)
# define ABOVE_L3(x, y) lower_than_mask(L[2][0]*x + L[2][2], L[2][1]*y)
# define ABOVE_L4(x, y) lower_than_mask(L[3][0]*x + L[3][2], L[3][1]*y)
# define ABOVE_L5(x, y) lower_than_mask(L[4][0]*x + L[4][2], L[4][1]*y)

# elif KYBER_K == 4

# define CODE_ALPHA (KYBER_Q/2)
# define CODE_BETA 436

static int32_t L[5][3] = {
    {-109, 416, 732626},
    {872, 3330, -2962321},
    {-416, 109, 369874},
    {-1, 1, -1229},
    {3330, 872, 58561},
};

# define ABOVE_L1(x, y) lower_than_mask(L[0][0]*x + L[0][2], L[0][1]*y)
# define ABOVE_L2(x, y) lower_than_mask(L[1][0]*x + L[1][2], L[1][1]*y)
# define ABOVE_L3(x, y) lower_than_mask(L[2][0]*x + L[2][2], L[2][1]*y)
# define ABOVE_L4(x, y) lower_than_mask(L[3][0]*x + L[3][2], L[3][1]*y)
# define ABOVE_L5(x, y) lower_than_mask(L[4][0]*x + L[4][2], L[4][1]*y)

# endif
///////////////////////////////////////////////////////////////////////////////////////////////////

# endif
```


Then we add the actual code for encoding and decoding. We also define the additional function `lower_than_mask`
used as a constant-time comparison between its inputs.

```C
# elif defined(KYBER_N_ENCODING_DIMENSIONS) && (KYBER_N_ENCODING_DIMENSIONS == 2)

// Returns 0xffffffff if (v1 < v2) and 0x00000000 otherwise
static __inline__ uint32_t lower_than_mask(const uint32_t v1, const uint32_t v2) {
  return -((v1 - v2) >> 31);
}

void poly_frommsg(poly *r, const uint8_t msg[KYBER_INDCPA_MSGBYTES])
{
  unsigned int i,j;

# if (KYBER_INDCPA_MSGBYTES != KYBER_N/8)
# error "KYBER_INDCPA_MSGBYTES must be equal to KYBER_N/33 bytes!"
# endif

  size_t base = 0;
  for(i=0;i<KYBER_N/8;i++) {
    for(j=0;j<8;j+=2) {
      uint16_t bit0_mask = -((msg[i] >> j) & 1);
      uint16_t bit1_mask = -((msg[i] >> (j + 1)) & 1);
      r->coeffs[base] = (CODE_ALPHA & bit1_mask) + (CODE_BETA & bit0_mask);
      r->coeffs[base + KYBER_N/2] = (CODE_BETA & bit1_mask) + (CODE_ALPHA & bit0_mask);

      base++;
    }
  }
}

static __inline__ int decode_msg_pair(int32_t x, int32_t y) {
  uint32_t reflect_mask = lower_than_mask(x, y);

  int32_t x_prime = (x & ~reflect_mask) | (y & reflect_mask);
  int32_t y_prime = (y & ~reflect_mask) | (x & reflect_mask);

  uint8_t above_l1 = ABOVE_L1(x_prime, y_prime);
  uint8_t above_l2 = ABOVE_L2(x_prime, y_prime);
  uint8_t above_l3 = ABOVE_L3(x_prime, y_prime);
  uint8_t above_l4 = ABOVE_L4(x_prime, y_prime);
  uint8_t above_l5 = ABOVE_L5(x_prime, y_prime);

  // Not needed, but conceptually: uint8_t c00 = (~above_l3 & above_l2 & above_l4);
  uint8_t c01 = ~above_l2 & ~above_l5 & ~above_l3;
  uint8_t c10 = above_l2 & above_l3 & ~above_l1;
  uint8_t c11 = above_l1 | (above_l3 & ~above_l2) | (above_l5 & ~above_l4);

  c01 &= (1 ^ reflect_mask);
  c10 &= (2 ^ reflect_mask);

  uint8_t bits = c01 | c10 | c11;

  return bits & 3;
}

/*************************************************
* Name:        poly_tomsg
*
* Description: Convert polynomial to 32-byte message
*
* Arguments:   - uint8_t *msg: pointer to output message
*              - const poly *a: pointer to input polynomial
**************************************************/
void poly_tomsg(uint8_t msg[KYBER_INDCPA_MSGBYTES], const poly *a)
{
  unsigned int i,j;
  uint32_t base = 0;
  for(i=0;i<KYBER_N/8;i++) {
    msg[i] = 0;
    for(j=0;j<8;j+=2) {
      int x = a->coeffs[base];
      int y = a->coeffs[base + KYBER_N/2];
      base++;

      msg[i] |= (decode_msg_pair(x, y) << j);
    }
  }
}

```

### Implementation of 4D codes over Kyber's reference/AVX2 code

Again, only the `poly.c` file was changed (both for REF and AVX2 implementations).
For the 4D encoding/decoding, we added the following:

```C

#elif defined(KYBER_N_ENCODING_DIMENSIONS) && (KYBER_N_ENCODING_DIMENSIONS == 4)

// FishCode of dimension 4 and generator matrix:
//     [   0    0 1664  752]
//     [ 752    0    0 1664]
//     [1664  752    0    0]
//     [   0 1664  752    0]

#define CODE_ALPHA 1664
#define CODE_BETA  752

static int16_t CODE_VALS[4] = {0, CODE_BETA, CODE_ALPHA, CODE_BETA + CODE_ALPHA - KYBER_Q};

static int8_t CODEWORDS[16][4] = {
   {0, 0, 0, 0},
   {1, 2, 0, 0},
   {2, 0, 0, 1},
   {3, 2, 0, 1},
   {0, 0, 1, 2},
   {1, 2, 1, 2},
   {2, 0, 1, 3},
   {3, 2, 1, 3},
   {0, 1, 2, 0},
   {1, 3, 2, 0},
   {2, 1, 2, 1},
   {3, 3, 2, 1},
   {0, 1, 3, 2},
   {1, 3, 3, 2},
   {2, 1, 3, 3},
   {3, 3, 3, 3},
};

// Returns 0xffffffff if (v1 < v2) and 0x00000000 otherwise
static __inline__ uint32_t lower_than_mask(const uint32_t v1, const uint32_t v2) {
  return -((v1 - v2) >> 31);
}

void poly_frommsg(poly *r, const uint8_t msg[KYBER_INDCPA_MSGBYTES])
{
  unsigned int i,j;

#if (KYBER_INDCPA_MSGBYTES != KYBER_N/8)
#error "KYBER_INDCPA_MSGBYTES must be equal to KYBER_N/8 bytes!"
#endif

  size_t base = 0;
  for(i=0;i<KYBER_N/8;i++) {
    for(j=0;j<8;j+=4) {
      uint16_t b3 = -((msg[i] >> (j + 0)) & 1);
      uint16_t b2 = -((msg[i] >> (j + 1)) & 1);
      uint16_t b1 = -((msg[i] >> (j + 2)) & 1);
      uint16_t b0 = -((msg[i] >> (j + 3)) & 1);

      r->coeffs[base]               = (CODE_ALPHA & b2) + (CODE_BETA & b3);
      r->coeffs[base + KYBER_N/4]   = (CODE_ALPHA & b3) + (CODE_BETA & b0);
      r->coeffs[base + 2*KYBER_N/4] = (CODE_ALPHA & b0) + (CODE_BETA & b1);
      r->coeffs[base + 3*KYBER_N/4] = (CODE_ALPHA & b1) + (CODE_BETA & b2);

      base++;
    }
  }
}

// Returns `abs(centered_mod(value, KYBER_Q)` assuming `-KYBER_Q <= value <= KYBER_Q`
static __inline__ int32_t abs_center_mod_of_2q_centered_value(int32_t value) {
  // value = abs(value):
  uint32_t mask_sign = value >> 31;
  value ^= mask_sign;
  value += mask_sign & 1;
  value -= KYBER_Q & lower_than_mask(KYBER_Q/2, value);
  return value;
}

static __inline__ uint32_t get_distance_sqr_to_codeword(uint16_t idx, uint32_t distsqr_matrix[4][4]) {
  int32_t a0 = distsqr_matrix[0][CODEWORDS[idx][0]];
  int32_t a1 = distsqr_matrix[1][CODEWORDS[idx][1]];
  int32_t a2 = distsqr_matrix[2][CODEWORDS[idx][2]];
  int32_t a3 = distsqr_matrix[3][CODEWORDS[idx][3]];

  // Returns `distance_sqr | codeword_index`
  return (a0 + a1 + a2 + a3) << 4 | idx;
}

// Computes min(v1, v2) in constant time
static __inline__ uint32_t secure_min(int32_t v1, int32_t v2) {
  uint32_t mask_min_v1 = lower_than_mask(v1, v2);
  return (mask_min_v1 & v1) | (~mask_min_v1 & v2);
}

static __inline__ int decode_msg_4d(int16_t target[4]) {
  // This function assumes that:
  //    * -q/2 < target[i] <= q/2
  //    * -q/2 < CODEWORD[idx][i] <= q/2

  // Build matrix with square distances to target coordinates
  uint32_t distsqr_matrix[4][4] = {0};
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      distsqr_matrix[i][j] = abs_center_mod_of_2q_centered_value(target[i] - CODE_VALS[j]);
      distsqr_matrix[i][j] *= distsqr_matrix[i][j];
    }
  }

  uint32_t min_dist_codeword = get_distance_sqr_to_codeword(0, distsqr_matrix);
  for (size_t i = 1; i < 16; i++) {
    min_dist_codeword = secure_min(get_distance_sqr_to_codeword(i, distsqr_matrix), min_dist_codeword);
  }

  // Retrieves `codeword_index` from `distance_sqr | codeword_index`
  return min_dist_codeword & 0xF;
}

/*************************************************
* Name:        poly_tomsg
*
* Description: Convert polynomial to 32-byte message
*
* Arguments:   - uint8_t *msg: pointer to output message
*              - const poly *a: pointer to input polynomial
**************************************************/
void poly_tomsg(uint8_t msg[KYBER_INDCPA_MSGBYTES], const poly *a)
{
  unsigned int i;
  uint32_t base = 0;
  for(i=0;i<KYBER_N/8;i++) {
    msg[i] = 0;

    int16_t target[4] = {0};
    target[0] = a->coeffs[base];
    target[1] = a->coeffs[base + 1*KYBER_N/4];
    target[2] = a->coeffs[base + 2*KYBER_N/4];
    target[3] = a->coeffs[base + 3*KYBER_N/4];
    msg[i] |= (decode_msg_4d(target) << 0);

    target[0] = a->coeffs[base + 1];
    target[1] = a->coeffs[base + 1 + 1*KYBER_N/4];
    target[2] = a->coeffs[base + 1 + 2*KYBER_N/4];
    target[3] = a->coeffs[base + 1 + 3*KYBER_N/4];
    msg[i] |= (decode_msg_4d(target) << 4);

    base += 2;
  }
}


#else
# error "KYBER_N_ENCODING_DIMENSIONS must be either 1, 2, or 4"
#endif // !defined(KYBER_N_ENCODING_DIMENSIONS) || (KYBER_N_ENCODING_DIMENSIONS == 1)

```

### Computation of the joint noise distribution `delta_m[i, i + n/2]`

The computation of the joint noise distribution is done by function `init_delta_m_as_joint_distribution`
from `kyber2d_dfr_analysis/joint_distribution_2d.c`, which is a good starting point to understand our code.
It is essentially the function below, but here we show a simplification without considering the generalized
compression we studied in a previous version of the paper.

While it is a large function, we tried to be very explicit on how the 1D distributions are
initialized, and how hey are combined to make the intermediate distributions finally
obtain `delta_m[i, i + n/2]`. Functions `fft2` and `ifft2` denote the FFT and Inverse-FFT
in 2 dimensions.

```C

// In this simplification, we assume that there is no generalized compression.
// That is `params->KYBER_N_BLOCKS_FOR_DU0 = KYBER_K`
int init_delta_m_as_joint_distribution(mpc_matrix_t *dist_delta_m, kyber_parameters_t *params) {
    init_fft_mpc();

    prob_dist_1d_t dist_dv;
    prob_dist_1d_init_as_compress_decompress_error(&dist_dv, params->KYBER_DV, params->KYBER_Q);

    prob_dist_1d_t dist_eta1;
    prob_dist_1d_init_as_centered_binomial(&dist_eta1, params->KYBER_ETA1);

    prob_dist_1d_t dist_eta2;
    prob_dist_1d_init_as_centered_binomial(&dist_eta2, params->KYBER_ETA2);

    prob_dist_1d_t dist_du0;
    prob_dist_1d_init_as_compress_decompress_error(&dist_du0, params->KYBER_DU0, params->KYBER_Q);

    prob_dist_1d_t dist_du0_sum_eta2;
    prob_dist_1d_init_as_sum(&dist_du0_sum_eta2, &dist_du0, &dist_eta2);

    prob_dist_1d_t dist_dv_sum_eta2;
    prob_dist_1d_init_as_sum(&dist_dv_sum_eta2, &dist_dv, &dist_eta2);

    mpc_matrix_t dist_s_du0_sum_eta2_product;
    get_base_2x2_distribution_for_self_convolution(&dist_s_du0_sum_eta2_product, &dist_du0_sum_eta2, &dist_eta1);
    fft2(&dist_s_du0_sum_eta2_product);
    self_convolution_in_fft_domain(&dist_s_du0_sum_eta2_product, 128 * params->KYBER_N_BLOCKS_FOR_DU0);

    mpc_matrix_t dist_s_du_sum_eta2_product;
    init_as_copy(&dist_s_du_sum_eta2_product, &dist_s_du0_sum_eta2_product);
    mpc_matrix_free(&dist_s_du0_sum_eta2_product);

    mpc_matrix_t dist_e_r_product;
    get_base_2x2_distribution_for_self_convolution(&dist_e_r_product, &dist_eta1, &dist_eta1);
    fft2(&dist_e_r_product);
    self_convolution_in_fft_domain(&dist_e_r_product, 128 * params->KYBER_K);

    mpc_matrix_t dist_sum_dot_products;
    init_as_convolution_in_fft_domain(&dist_sum_dot_products, &dist_e_r_product, &dist_s_du_sum_eta2_product);

    mpc_matrix_free(&dist_e_r_product);
    mpc_matrix_free(&dist_s_du_sum_eta2_product);

    mpc_matrix_t dist_dv_sum_eta2_as_2d;
    init_as_2d_from_1d(&dist_dv_sum_eta2_as_2d, &dist_dv_sum_eta2);
    fft2(&dist_dv_sum_eta2_as_2d);

    init_as_convolution_in_fft_domain(dist_delta_m, &dist_sum_dot_products, &dist_dv_sum_eta2_as_2d);
    ifft2(dist_delta_m);

    // Now dist_delta_m contains the 2D error distribution of delta_m[i, i + n/2]

    prob_dist_1d_clear(&dist_du0);
    prob_dist_1d_clear(&dist_dv);
    prob_dist_1d_clear(&dist_eta1);
    prob_dist_1d_clear(&dist_eta2);
    prob_dist_1d_clear(&dist_du0_sum_eta2);
    prob_dist_1d_clear(&dist_dv_sum_eta2);

    mpc_matrix_free(&dist_sum_dot_products);
    mpc_matrix_free(&dist_dv_sum_eta2_as_2d);

    clear_fft_mpc();

    return 0;
}
```


# Description of the results

Our results can be found in the `results` directory, whose structure is given below (omitting the
`figures` directory).

```
$ tree results
results
├── kyber2d_dfr # DFR results using the efficient decoder based on approximate Voronoi cells (default)
│   ├── result_kyber2d_dfr_effect_of_dv_in_level5.csv
│   └── result_kyber2d_dfr_multiple_parameters.csv
├── kyber2d_dfr_mindist_decoding # DFR results using the minimum distance decoder (old setup, put here for reference)
│   ├── result_kyber2d_dfr_effect_of_dv_in_level5.csv
│   └── result_kyber2d_dfr_multiple_parameters.csv
├── kyber2d_joint_error_distribution
│   ├── delta_m_pair_q=3329_l=1_k=2_e1=3_e2=2_du0rep=2_du0=10_du1=0_dv=4.csv
│   ├── delta_m_pair_q=3329_l=3_k=3_e1=2_e2=2_du0rep=3_du0=10_du1=0_dv=4.csv
│   ├── delta_m_pair_q=3329_l=5_k=4_e1=2_e2=2_du0rep=4_du0=11_du1=0_dv=5.csv
│   ├── delta_m_pair_q=3329_l=.tar.xz
│   ├── kyber2d_joint_error_distribution.zip
│   ├── params.csv_q=3329_l=1_k=2_e1=3_e2=2_du0rep=2_du0=10_du1=0_dv=4.csv
│   ├── params.csv_q=3329_l=3_k=3_e1=2_e2=2_du0rep=3_du0=10_du1=0_dv=4.csv
│   └── params.csv_q=3329_l=5_k=4_e1=2_e2=2_du0rep=4_du0=11_du1=0_dv=5.csv
├── performance
│   ├── test_speed1024_2d_ref.dat
│   ├── test_speed1024_ref.dat
│   ├── test_speed512_2d_ref.dat
│   ├── test_speed512_ref.dat
│   ├── test_speed768_2d_ref.dat
│   └── test_speed768_ref.dat
└── toy_example_of_main_lemma
    └── toy_example_of_main_lemma.txt

```

In high-level, the results presented in each subdirectory can be described as follows.

* `kyber2d_dfr`: CSV files with selected Kyber parameters using different 2D codes and their
                 decryption failure rates (DFR). All results here consider the **efficient decoder**
                 based on approximate Voronoi cells described in the paper.
* `kyber2d_dfr_mindist_decoding`: Similar to `kyber2d_dfr` but using **
* `kyber2d_joint_error_distribution`: CSV files with 2D distributions of `delta_m[i, i + n/2]`
                                      for the 3 Kyber original parameter sets. These are only used
                                      for plotting Figure 1.
* `performance`: Performance results for comparing our implementation of 2D codes with Kyber's
                 reference implementation.
* `toy_example_of_main_lemma`: The output of our Python script showing a toy example of the lemma
                               allowing the 2D distribution to be computing with a series of
                               convolutions.


Most results are associated with experiments, whose setup files and scripts are
in the `experiments` directory, listed below.
```
experiments
├── scripts
│   ├── generate_decoding_lines.py
│   ├── run_kyber2d_dfr_effect_of_dv_in_level5.sh
│   ├── run_kyber_2d_dfr_joint_error_distribution.sh
│   ├── run_kyber2d_dfr_multiple_parameters.sh
│   └── run_speed_tests.sh
└── setup
    ├── exp_kyber2d_dfr_effect_of_dv_in_level5.csv
    ├── exp_kyber2d_dfr_multiple_parameters.csv
    ├── exp_kyber2d_joint_error_distribution.csv
    └── gen_kyber2d_dfr_multiple_parameters_setup.py
```

In the following, we explain how each experiment is defined and provide intuition on
what is the expected result.

**Important notice.** All reproduction steps are taken care when you call the `./reproduce_me.sh` main script.
Therefore, this section explains the internal scripts called from the main script.

**Script parameters.** All scripts except for the main one require two arguments:

* `<reproduction_directory>`: the directory where it will write the data (possibly creating a subdirectory in it)
* `<build_directory>`: the directory used as the build directory by `cmake`


## Performance results in `results/performance`

These are the easiest results to interpret and to generate.
Simply run
```
$ ./experiments/scripts/run_speed_tests.sh <reproduction_directory> <build_directory>
```

For us, the relevant parts of these files are the functions `poly_frommg` (encoding), `poly_tomsg` (decoding),
and `kyber_decaps` (full decapsulation).
For Kyber1024 (level5), these can be compared using the following command lines.

```
$ grep -A2 -E  "frommsg|tomsg|decaps" results/performance/avx2/*1024*d.txt
results/performance/avx2/test_speed1024_1d.txt:poly_tomsg:
results/performance/avx2/test_speed1024_1d.txt-median: 17 cycles/ticks
results/performance/avx2/test_speed1024_1d.txt-average: 16 cycles/ticks
--
results/performance/avx2/test_speed1024_1d.txt:poly_frommsg:
results/performance/avx2/test_speed1024_1d.txt-median: 23 cycles/ticks
results/performance/avx2/test_speed1024_1d.txt-average: 23 cycles/ticks
--
results/performance/avx2/test_speed1024_1d.txt:kyber_decaps:
results/performance/avx2/test_speed1024_1d.txt-median: 46104 cycles/ticks
results/performance/avx2/test_speed1024_1d.txt-average: 46490 cycles/ticks
--
results/performance/avx2/test_speed1024_2d.txt:poly_tomsg:
results/performance/avx2/test_speed1024_2d.txt-median: 404 cycles/ticks
results/performance/avx2/test_speed1024_2d.txt-average: 405 cycles/ticks
--
results/performance/avx2/test_speed1024_2d.txt:poly_frommsg:
results/performance/avx2/test_speed1024_2d.txt-median: 78 cycles/ticks
results/performance/avx2/test_speed1024_2d.txt-average: 78 cycles/ticks
--
results/performance/avx2/test_speed1024_2d.txt:kyber_decaps:
results/performance/avx2/test_speed1024_2d.txt-median: 46486 cycles/ticks
results/performance/avx2/test_speed1024_2d.txt-average: 46608 cycles/ticks
--
results/performance/avx2/test_speed1024_4d.txt:poly_tomsg:
results/performance/avx2/test_speed1024_4d.txt-median: 990 cycles/ticks
results/performance/avx2/test_speed1024_4d.txt-average: 1002 cycles/ticks
--
results/performance/avx2/test_speed1024_4d.txt:poly_frommsg:
results/performance/avx2/test_speed1024_4d.txt-median: 65 cycles/ticks
results/performance/avx2/test_speed1024_4d.txt-average: 78 cycles/ticks
--
results/performance/avx2/test_speed1024_4d.txt:kyber_decaps:
results/performance/avx2/test_speed1024_4d.txt-median: 46092 cycles/ticks
results/performance/avx2/test_speed1024_4d.txt-average: 46961 cycles/ticks

```


## The DFR results

The DFR results in `results/*dfr*` are obtained using:
```
$ ./experiments/scripts/run_exp_kyber2d.sh <reproduction_directory> <build_directory>
$ ./experiments/scripts/run_exp_kyber4d.sh <reproduction_directory> <build_directory>
```

### Case 2D:

The setup file of the main experiment is CSV `experiments/setup/exp_kyber2d_dfr`.
It defines different Kyber parameters for which we want to compute the 2D DFR of different codes.
```
KYBER_Q,KYBER_LEVEL,KYBER_K,KYBER_ETA1,KYBER_ETA2,KYBER_N_BLOCKS_FOR_DU0,KYBER_DU0,KYBER_DU1,KYBER_DV
3329,1,2,3,2,2,9,0,3
3329,1,2,3,2,2,9,0,4
3329,1,2,3,2,2,9,0,5
3329,1,2,3,2,2,10,0,3
3329,1,2,3,2,2,10,0,4
3329,3,3,2,2,3,9,0,3
3329,3,3,2,2,3,9,0,4
3329,3,3,2,2,3,9,0,5
3329,3,3,2,2,3,10,0,3
3329,3,3,2,2,3,10,0,4
3329,5,4,2,2,4,10,0,4
3329,5,4,2,2,4,10,0,5
3329,5,4,2,2,4,10,0,6
3329,5,4,2,2,4,11,0,4
3329,5,4,2,2,4,11,0,5
3329,1,2,3,0,2,9,0,5
3329,3,3,2,0,3,9,0,5
3329,3,3,2,1,3,10,0,3
```

The output is another CSV file in `results/kyber2d_dfr/result_kyber2d_dfr.csv`.
Notice that this is a large file because we evaluate the DFR for each possible `beta` parameter from 0 to 500.
```
KYBER_Q,KYBER_LEVEL,KYBER_K,KYBER_ETA1,KYBER_ETA2,KYBER_N_BLOCKS_FOR_DU0,KYBER_DU0,KYBER_DU1,KYBER_DV,CODE_ALPHA,CODE_BETA,DFR
3329,1,2,3,2,2,9,0,3,1664,0,1.0969336836229559799652030636118520979595057895107241911494817822181072339490170e-18
3329,1,2,3,2,2,9,0,3,1664,1,1.0971363028130969128667818258777957062814933486214388685626861625927579537404139e-18
3329,1,2,3,2,2,9,0,3,1664,2,1.0970478752223682701615992037155600973644802028889651646078901793070622833550620e-18
3329,1,2,3,2,2,9,0,3,1664,3,1.0970066460652067508910572566302821608013107114161830079950657655535862843832733e-18
3329,1,2,3,2,2,9,0,3,1664,4,1.0971331973525063043929074572682580577188056708791153016391053727131528442383535e-18
3329,1,2,3,2,2,9,0,3,1664,5,1.0973353392738909969494740072403879017610226024478546459453498365895142190896346e-18
3329,1,2,3,2,2,9,0,3,1664,6,1.0975652998794586038986772203644061864531717744534631454745129038568145068817234e-18
3329,1,2,3,2,2,9,0,3,1664,7,1.0978265847334240171373451042252188327898819869331151053590190083111761243067297e-18
3329,1,2,3,2,2,9,0,3,1664,8,1.0981344777247930828887138732155327219628471214304809033957031009755928496947822e-18
...
```

Notice that there are other setup files in `experiments/setup/exp_kyber2d*`, for computing the
effect of `dv` in the DFR and also to consider pk compression.

### Case 4D:

The setup file of the main experiment is CSV `experiments/setup/exp_kyber4d_dfr.csv`.
It defines different Kyber parameters for which we want to compute the 4D DFR of different codes.
```
KYBER_Q,KYBER_LEVEL,KYBER_K,KYBER_ETA1,KYBER_ETA2,KYBER_N_BLOCKS_FOR_DU0,KYBER_DU0,KYBER_DU1,KYBER_DV,CODE_BETA,CODE_PNORM
3329,1,2,3,2,2,9,0,3,688,2.3833461905097693
3329,1,2,3,2,2,9,0,4,731,2.1171550193952813
3329,1,2,3,2,2,9,0,5,751,2.0000003890434996
3329,1,2,3,0,2,9,0,5,745,2.0312567350190927
3329,1,2,3,2,2,10,0,3,660,2.5647353744238117
3329,1,2,3,2,2,10,0,4,717,2.197292261676462
3329,3,3,2,2,3,9,0,3,680,2.430772527953634
3329,3,3,2,2,3,9,0,4,725,2.151345091917993
3329,3,3,2,2,3,9,0,5,746,2.0263103836432723
3329,3,3,2,2,3,10,0,3,646,2.663492196293784
3329,3,3,2,2,3,10,0,4,709,2.249047620905394
3329,3,3,2,0,3,9,0,5,741,2.055417880003882
3329,3,3,2,1,3,10,0,3,633,2.7558103094443007
3329,5,4,2,2,4,10,0,4,715,2.213894798056261
3329,5,4,2,2,4,10,0,5,741,2.057696900078069
3329,5,4,2,2,4,10,0,6,751,2.0000003890434996
3329,5,4,2,2,4,11,0,4,708,2.252335697078637
3329,5,4,2,2,4,11,0,5,738,2.075737935656825
3329,1,2,3,3,2,10,0,4,722,2.169739739983144
```

Notice that this is different to the setup files for 2D because now there are the `CODE_BETA` and `CODE_PNORM`
columns. The `CODE_PNORM` column is not used by the DFR computation, which needs only `CODE_BETA`.

**Important:** This file is computed using the script `analysis_script/minal.py` with another input
file `experiments/setup/exp_kyber4d_dfr_without_betas.csv` that does not have the values of `CODE_BETA` as follows:

```
./analysis_scripts/minal.py > experiments/setup/exp_kyber4d_dfr.csv
```

Furthermore, the 4D DFR experiments requires we define the set of DFR relevant vectors in
`minal_src/dfr_computation/minal4d_defs/dfr_relevant_points_for_selected_betas.h`.
Naturally, this file is already precomputed for the values of `beta` that appear in the current experiment setup.
However, if you want to change the betas, you'll have to open call the following method from the `MinalCode` class
(defined in `analysis_scripts/minal.py`):
```Python
MinalCode.generate_4d_dfr_relevant_vectors_for_betas(filepath, kyber_params_csv='experiments/setup/exp_kyber4d_dfr.csv'):
```
Where `filepath = minal_src/dfr_computation/minal4d_defs/dfr_relevant_points_for_selected_betas.h`.



The output of the experiment is another CSV file in `results/kyber4d_dfr/result_kyber4d_dfr.csv`.
```
KYBER_Q,KYBER_LEVEL,KYBER_K,KYBER_ETA1,KYBER_ETA2,KYBER_N_BLOCKS_FOR_DU0,KYBER_DU0,KYBER_DU1,KYBER_DV,CODE_DIMENSION,CODE_ALPHA,CODE_BETA,DFR
3329,1,2,3,2,2,9,0,3,4,1664,688,1.4279506510673429876482040194918633985730564789759112873936091531170520531596806e-19
3329,1,2,3,2,2,9,0,4,4,1664,731,9.7873418752094013301438754282503659122968650543117638388865987214534685855249170e-28
3329,1,2,3,2,2,9,0,5,4,1664,751,5.2000714461584772231620704135909689115639911504352690177353419432259612862561077e-32
3329,1,2,3,0,2,9,0,5,4,1664,745,1.6360071698666228250869396813709177312917621323931065961117992994996678763833795e-38
3329,1,2,3,2,2,10,0,3,4,1664,660,3.4310041565196695999889103809953859314275057700027408084422422577938516680539371e-33
3329,1,2,3,2,2,10,0,4,4,1664,717,6.9549458696223591492135636006871917838528459105550893663546822664738929323384240e-47
3329,3,3,2,2,3,9,0,3,4,1664,680,1.6959376067939888093219612267332078340267347769756510080465600344103623244532892e-21
3329,3,3,2,2,3,9,0,4,4,1664,725,1.2665340648989775722833448794547715808154239317824817254250852705822281081569391e-30
3329,3,3,2,2,3,9,0,5,4,1664,746,1.6805197596579943260731528381570364174259794371018543792880189983562049514113176e-35
3329,3,3,2,2,3,10,0,3,4,1664,646,5.0588838697119121803589548875312417523584347898075088979250404231861828812669927e-39
3329,3,3,2,2,3,10,0,4,4,1664,709,7.2609736730825394396320058674186569486541114615913912556345853421252153571551832e-56
3329,3,3,2,0,3,9,0,5,4,1664,741,2.3141296465065289581571466822156211860334608525405051260170711673708017828265308e-43
3329,3,3,2,1,3,10,0,3,4,1664,633,2.3414934375386296906921358275621589929749549417807619512518483081020326052027107e-46
3329,5,4,2,2,4,10,0,4,4,1664,715,1.6342618222903739440755855212782258221257670652877012551020678265732222639820731e-43
3329,5,4,2,2,4,10,0,5,4,1664,741,1.1114050961338923075309020517591102637976159398421661887403607300371799290970856e-50
3329,5,4,2,2,4,10,0,6,4,1664,751,3.5422485814451757641398595499253616952750885587922096961586130579600057398419132e-54
3329,5,4,2,2,4,11,0,4,4,1664,708,1.0008500485548108749790834729011957800984969511838091924706070588024847253315466e-52
3329,5,4,2,2,4,11,0,5,4,1664,738,1.7545790099057973553165130448089826942862004809659631486595520259252652510052882e-61
3329,1,2,3,3,2,10,0,4,4,1664,722,1.4875435433593403537548938165119651290381802529642758701758077168982064552987185e-41
```

Notice that there are other setup files in `experiments/setup/exp_kyber4d*`, for computing the
effect of `dv` in the DFR and also to consider pk compression.


# Reproducing the figures from the paper

To generate the paper's figures, we use a single graphing tool: `analysis_scripts/generate_graphs.py`.
This script uses typical depends on common python libraries such as `matplotlib`, `seaborn`, `pandas`, `numpy` and
`scipy`.

If you have have these libraries installed, you may be ready to run the script. If you don't
or if you have problems running the script, we provide a `Pipfile` to help you get the correct
versions of the libraries. In this case, you'll need to install [`Pipenv`](https://pipenv.pypa.io/en/latest/).
Then, simply run:

```
pipenv shell
pipenv install
./analysis_scripts/generate_graphs.py figures
```

This will generate the paper's figures in the `figures` directory.
Below are the generated figures and the number of the figure in the paper.
```
$ tree results/figures
results/figures
├── codes
│   ├── candidate_code_pnorm2_446_1664.pdf
│   ├── candidate_code_pnorm3_339_1664.pdf
│   ├── dist_shape_level5.pdf
│   ├── kyber_code_2cols.pdf
│   ├── kyber_code_decoding_2cols.pdf
│   ├── kyber_code_decoding_symmetry_2cols.pdf
│   └── lattice_code_2cols.pdf
└── compression
    ├── balancing_du_and_dv.pdf
    └── effect_of_dv_in_level5.pdf

3 directories, 9 files

```

Note that the file is large, since the figures require lots of tuning.
However, the `main` function can be a useful guide to explore how each function relates to
the figures in the paper.

```Python
def main():

    if len(sys.argv) != 2:
        print(f'Usage {sys.argv[0]} figures_directory')
        sys.exit(1)

    figures_dir = sys.argv[1]
    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(os.path.join(figures_dir, 'codes'), exist_ok=True)
    os.makedirs(os.path.join(figures_dir, 'compression'), exist_ok=True)

    # Prepares matplotlib with lcns fonts and more visually appealing settings
    latexify_tches()

    # Figures 1a, 1b, and 1c
    CodePlots.draw_codes(figures_dir)

    # Figure 2
    DFRPlots.plot_effect_of_dimension_and_beta(figures_dir)

    # Figure 3
    DFRPlots.plot_effect_of_dv_dense(figures_dir)

    # Figure 4
    DistributionPlots.plot_shape_of_distributions_different_dv(figures_dir)

    # Figures 5a, 5b, and 5c
    CodePlots.draw_pair_of_good_codes_voronoi_pnorm(figures_dir)

    # Figures 6a and 6b
    CodePlots.draw_code_decoding_voronoi(figures_dir)

    # Figure 7
    DFRPlots.plot_compression_simple_factors_2d_4d(figures_dir)

```

# License

Apache 2.0

> Copyright 2025 LG Electronics, Inc. All Rights Reserved.
> SPDX-License-Identifier: Apache-2.0
