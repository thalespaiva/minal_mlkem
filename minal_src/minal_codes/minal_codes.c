// Copyright 2025 LG Electronics, Inc. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "cpucycles.h"

#define N_CODE_VALS 4

typedef struct minal_params_s {
    int q;
    int beta;
    int alpha;
    int n_dim;
    int16_t **CODEWORDS; // CODEWORDS[2^(n_dim)][n_dim]
    int16_t CODE_VALS[N_CODE_VALS];
} minal_params_t;


void print_vector(int16_t vector[], int n_dim);
void print_vector8(uint8_t vector[], int n_dim);


void fake_encode(int16_t codeword[], uint8_t msg_bits[], minal_params_t *minal) {

    memset(codeword, 0, minal->n_dim * sizeof(*codeword));

    if (minal->n_dim == 1) {
        codeword[0] = msg_bits[0] * 1;
        return;
    }

    for (int i = 0; i < minal->n_dim; i++) {
        int idx = (minal->n_dim - 2 + i) % minal->n_dim;
        codeword[i] = msg_bits[idx] * 2 + msg_bits[(idx + 1) % minal->n_dim] * 1;
    }
}

void minal_code_init(minal_params_t *minal, int q, int beta, int n_dim) {
    minal->q = q;
    minal->alpha = q/2;
    minal->beta = beta;
    minal->n_dim = n_dim;

    int n = 1 << minal->n_dim;

    minal->CODE_VALS[0] = 0;
    minal->CODE_VALS[1] = minal->beta;
    minal->CODE_VALS[2] = minal->alpha;
    minal->CODE_VALS[3] = minal->alpha + minal->beta - minal->q;

    minal->CODEWORDS = malloc(n * sizeof(*minal->CODEWORDS));

    for (int i = 0; i < n; i++) {
        minal->CODEWORDS[i] = malloc(n_dim * sizeof(*minal->CODEWORDS[i]));

        uint8_t msg_bits[minal->n_dim];
        for (int k = 0; k < minal->n_dim; k++) {
            msg_bits[k] = (i >> (minal->n_dim - 1 - k)) & 1;
        }
        memset(minal->CODEWORDS[i], 0, minal->n_dim * sizeof(*minal->CODEWORDS[i]));
        fake_encode(minal->CODEWORDS[i], msg_bits, minal);
        // printf("msg_bits = ");
        // print_vector8(msg_bits, minal->n_dim);
        // printf(" -> ");
        // print_vector(minal->CODEWORDS[i], minal->n_dim);
        // printf("\n");
    }

}

void minal_code_clear(minal_params_t *minal) {
    minal->CODE_VALS[0] = -1;
    minal->CODE_VALS[1] = -1;
    minal->CODE_VALS[2] = -1;
    minal->CODE_VALS[3] = -1;

    int n = 1 << minal->n_dim;
    for (int i = 0; i < n; i++) {
        free(minal->CODEWORDS[i]);
    }
    free(minal->CODEWORDS);

    minal->q = -1;
    minal->alpha = -1;
    minal->beta = -1;
    minal->n_dim = -1;

}

// Returns 0xffffffff if (v1 < v2) and 0x00000000 otherwise
static __inline__ uint32_t lower_than_mask(const uint32_t v1, const uint32_t v2) {
  return -((v1 - v2) >> 31);
}

static __inline__ uint64_t lower_than_mask64(const uint64_t v1, const uint64_t v2) {
  return -((v1 - v2) >> 63);
}

static __inline__ void minal_code_encode(int16_t codeword[], uint8_t msg_bits[], minal_params_t *minal) {
    int alpha = minal->q / 2;

    memset(codeword, 0, minal->n_dim * sizeof(*codeword));
    if (minal->n_dim == 1) {
        codeword[0] = msg_bits[0] * alpha;
        return;
    }

    for (int i = 0; i < minal->n_dim; i++) {
        int idx = (minal->n_dim - 2 + i) % minal->n_dim;
        codeword[i] = msg_bits[idx] * alpha + msg_bits[(idx + 1) % minal->n_dim] * minal->beta;
    }
}

// Returns `abs(centered_mod(value, KYBER_Q)` assuming `-KYBER_Q <= value <= KYBER_Q`
static __inline__ int32_t abs_center_mod_of_2q_centered_value(int32_t value, int q) {
  // value = abs(value):
  uint32_t mask_sign = value >> 31;
  value ^= mask_sign;
  value += mask_sign & 1;
  value -= q & lower_than_mask(q/2, value);
  return value;
}

uint64_t get_distance_sqr_to_codeword_large_dim(uint16_t idx, uint32_t distsqr_matrix[][N_CODE_VALS],
                                                        minal_params_t *minal) {

    uint64_t dist_sqr = 0;
    for (int i = 0; i < minal->n_dim; i++) {
        dist_sqr += distsqr_matrix[i][minal->CODEWORDS[idx][i]];
    }

    // Returns `distance_sqr | codeword_index`
    return (dist_sqr) << 16 | idx;
}

uint32_t get_distance_sqr_to_codeword_small_dim(uint16_t idx, uint32_t distsqr_matrix[][N_CODE_VALS],
                                                        minal_params_t *minal) {

    uint64_t dist_sqr = 0;
    for (int i = 0; i < minal->n_dim; i++) {
        dist_sqr += distsqr_matrix[i][minal->CODEWORDS[idx][i]];
    }

    // Returns `distance_sqr | codeword_index`
    return (dist_sqr) << 8 | idx;
}


// Computes min(v1, v2) in constant time
static __inline__ uint32_t secure_min(int32_t v1, int32_t v2) {
    uint32_t mask_min_v1 = lower_than_mask(v1, v2);
    return (mask_min_v1 & v1) | (~mask_min_v1 & v2);
}

static __inline__ uint64_t secure_min64(int64_t v1, int64_t v2) {
    uint64_t mask_min_v1 = lower_than_mask64(v1, v2);
    return (mask_min_v1 & v1) | (~mask_min_v1 & v2);
}



static __inline__ uint16_t minal_code_decode_large_dim(int16_t target[], minal_params_t *minal) {
    // Build matrix with square distances to target coordinates
    uint32_t distsqr_matrix[minal->n_dim][N_CODE_VALS];
    for (int i = 0; i < minal->n_dim; i++) {
        for (int j = 0; j < N_CODE_VALS; j++) {
            distsqr_matrix[i][j] = abs_center_mod_of_2q_centered_value(target[i] - minal->CODE_VALS[j], minal->q);
            distsqr_matrix[i][j] *= distsqr_matrix[i][j];
        }
    }

    uint64_t min_dist_codeword = get_distance_sqr_to_codeword_large_dim(0, distsqr_matrix, minal);
    for (int i = 1; i < (1 << minal->n_dim); i++) {
        min_dist_codeword = secure_min64(get_distance_sqr_to_codeword_large_dim(i, distsqr_matrix, minal), min_dist_codeword);
    }
    // Retrieves `codeword_index` from `distance_sqr | codeword_index`

    // return min_dist_codeword & 0xFF;
    return min_dist_codeword & 0xFFFF;
}

static __inline__ uint16_t minal_code_decode_small_dim(int16_t target[], minal_params_t *minal) {
    // Build matrix with square distances to target coordinates
    uint32_t distsqr_matrix[minal->n_dim][N_CODE_VALS];
    for (int i = 0; i < minal->n_dim; i++) {
        for (int j = 0; j < N_CODE_VALS; j++) {
            distsqr_matrix[i][j] = abs_center_mod_of_2q_centered_value(target[i] - minal->CODE_VALS[j], minal->q);
            distsqr_matrix[i][j] *= distsqr_matrix[i][j];
        }
    }

    uint32_t min_dist_codeword = get_distance_sqr_to_codeword_small_dim(0, distsqr_matrix, minal);
    for (int i = 1; i < (1 << minal->n_dim); i++) {
        min_dist_codeword = secure_min(get_distance_sqr_to_codeword_small_dim(i, distsqr_matrix, minal), min_dist_codeword);
    }
    // Retrieves `codeword_index` from `distance_sqr | codeword_index`

    // return min_dist_codeword & 0xFF;
    return min_dist_codeword & 0xFF;
}


static __inline__ uint16_t minal_code_decode(int16_t target[], minal_params_t *minal) {
    // This function assumes that:
    //    * -q/2 < target[i] <= q/2
    //    * -q/2 < CODEWORD[idx][i] <= q/2
    if (minal->n_dim <= 5) {
        return minal_code_decode_small_dim(target, minal);
    }
    return minal_code_decode_large_dim(target, minal);
}





void print_vector8(uint8_t vector[], int n_dim) {
    printf("[ ");
    for (int i = 0; i < n_dim; i++) {
        printf("%d ", vector[i]); fflush(stdout);
    }
    printf("]");
}

void print_vector(int16_t vector[], int n_dim) {
    printf("[ ");
    for (int i = 0; i < n_dim; i++) {
        printf("%d ", vector[i]);
    }
    printf("]");
}


#define NTESTS 10000

void test_encoding_decoding(minal_params_t *minal) {

    for (int i = 0; i < (1 << minal->n_dim); i++) {
        uint8_t msg_bits[minal->n_dim];

        for (int j = 0; j < minal->n_dim; j++) {
            msg_bits[j] = (i >> (minal->n_dim - 1 - j)) & 1;
        }
        int16_t codeword[minal->n_dim];
        minal_code_encode(codeword, msg_bits, minal);
        uint16_t dec = minal_code_decode(codeword, minal);
        if(dec != i) {
            printf("dec = %d, i = %d\n", dec, i);
            fflush(stdout);
        }
        assert(dec == i);
    }

}

void test_encoding_decoding_speed(minal_params_t *minal) {

    for (int i = 0; i < 2; i++) {
        uint8_t msg_bits[minal->n_dim];

        for (int j = 0; j < minal->n_dim; j++) {
            msg_bits[j] = (i >> (minal->n_dim - 1 - j)) & 1;
        }
        int16_t codeword[minal->n_dim];
        minal_code_encode(codeword, msg_bits, minal);
        uint16_t dec = minal_code_decode(codeword, minal);
        if(dec != i) {
            printf("dec = %d, i = %d\n", dec, i);
            fflush(stdout);
        }
        assert(dec == i);
    }

}


#define MIN_DIMENSION 2
// #define MAX_DIMENSION 16
#define MAX_DIMENSION 14

// Using `minal.py`:
// BETAS_USED_FOR_TESTING_FOR_DIMENSION = [-1, -1] + [MinalCode.find_best_beta(3329, d) for d in range(2, 9)]
static int BETAS_USED_FOR_TESTING_FOR_DIMENSION[MAX_DIMENSION + 1] = {-1, -1, 446, 636, 751, 832, 981, 981,981,981,981,981,981,981};


int main(void) {
    uint64_t t[NTESTS];


    printf("n_dimensions,decode_median_cycles_per_bit,decode_average_cycles_per_bit\n");

    for (int n_dim = MIN_DIMENSION; n_dim <= MAX_DIMENSION; n_dim++) {

        minal_params_t minal;
        minal_code_init(&minal, 3329, BETAS_USED_FOR_TESTING_FOR_DIMENSION[n_dim], n_dim);

        // printf("dim = %d\n", minal.n_dim);
        // printf("Testing dim = %d...", n_dim); fflush(stdout);
        test_encoding_decoding(&minal);
        // printf("Passed!\n"); fflush(stdout);

        uint8_t msg_bits[minal.n_dim];

        int16_t msg = rand() % (1 << n_dim);
        for (int j = 0; j < minal.n_dim; j++) {
            msg_bits[j] = (msg >> (minal.n_dim - 1 - j)) & 1;
        }
        int16_t codeword[minal.n_dim];
        minal_code_encode(codeword, msg_bits, &minal);

        uint16_t dec;
        for(int i = 0; i < 10000; i++) {
            t[i] = cpucycles();
            dec = minal_code_decode(codeword, &minal);
        }
        assert(dec == msg);
        // print_results("decode_msg: ", t, 10000);
        uint64_t median;
        uint64_t average;
        get_results(&average, &median, t, 10000);
        printf("%d,%.2lf,%.2lf\n", n_dim, (double) median/n_dim, (double) average/n_dim);

        minal_code_clear(&minal);
    }
    return 0;
}


