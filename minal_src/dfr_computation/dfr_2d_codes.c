// Copyright 2025 LG Electronics, Inc. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <mpfr.h>
#include <sys/resource.h>
#include <omp.h>
#include <stdint.h>

#include "mp_matrix.h"
#include "common.h"
#include "decoding_lines.h"

#define STR_MAX 1000

#define N_CODEWORDS 4

typedef struct point_s {
    int x[2];
} point_t;

void print_point(point_t *point) {
    printf("(%d, %d)\n", point->x[0], point->x[1]);
}

static __inline__ int get_distance_sqr(point_t *a, point_t *b) {
    int x0_diff = a->x[0] - b->x[0];
    int x1_diff = a->x[1] - b->x[1];
    return x0_diff * x0_diff + x1_diff * x1_diff;
}

static __inline__ uint32_t lower_than_mask(const uint32_t v1, const uint32_t v2) {
  return -((v1 - v2) >> 31);
}

static __inline__ int64_t get_distance_sqr_reduced_points(point_t *a, point_t *b, int kyber_q) {

    int64_t min_dist_sqr = kyber_q * kyber_q;

    #pragma GCC unroll 10
    for (int x0 = -kyber_q; x0 <= kyber_q; x0 += kyber_q) {
        #pragma GCC unroll 10
        for (int x1 = -kyber_q; x1 <= kyber_q; x1 += kyber_q) {
            point_t tmp_a = {{a->x[0] + x0, a->x[1] + x1}};
            int64_t dist_sqr = get_distance_sqr(&tmp_a, b);
            min_dist_sqr = dist_sqr < min_dist_sqr ? dist_sqr : min_dist_sqr;
        }
    }
    return min_dist_sqr;
}

int find_closest_codeword(point_t codewords[], point_t *target_brute, int kyber_q) {

    point_t target = {{MOD(target_brute->x[0], kyber_q), MOD(target_brute->x[1], kyber_q)}};

    int64_t min_dist_sqr = kyber_q * kyber_q;
    int closest_codeword_index = -1;
    #pragma GCC unroll 10
    for (int i = 0; i < N_CODEWORDS; i++) {
        int64_t dist_sqr = get_distance_sqr_reduced_points(&target, &codewords[i], kyber_q);
        if (dist_sqr < min_dist_sqr) {
            closest_codeword_index = i;
            min_dist_sqr = dist_sqr;
        }
    }

    return closest_codeword_index;
}


// This is a general decoding algorithm, that, although less inefficient, can be applied to multiple parameter
// sets without much trouble. In particular, we don't need decoding lines or to know other geometrical properties
// of the code.
void get_dfr_for_2d_code_with_minimum_distance_decoding(mpfr_t error_prob, mpfr_matrix_t *delta_m, int code_alpha, int code_beta, int kyber_q) {

    int nthreads = omp_get_max_threads();

    assert(0 <= code_alpha);
    assert(code_alpha < kyber_q);
    assert(0 <= code_beta);
    assert(code_beta < kyber_q);
    assert(code_beta + code_beta < kyber_q);

    point_t points[N_CODEWORDS] = {
        {{0, 0}},
        {{code_beta, code_alpha}},
        {{code_alpha, code_beta}},
        {{code_alpha + code_beta, code_alpha + code_beta}},
    };


    mpfr_t pair_error_prob_for_thread[nthreads];
    for (int t = 0; t < nthreads; t++) {
        mpfr_init2(pair_error_prob_for_thread[t], N_PRECISION_BITS);
        mpfr_set_d(pair_error_prob_for_thread[t], 0, MPFR_RNDD);
    }

    for (int j = 0; j < N_CODEWORDS; j++) {

        // Iterates over possible pairs of errors (e0, e1)
        #pragma omp parallel for schedule(static) // Using static scheduling the results are reproducible
        for (int a0 = 0; a0 < N_SUPPORT_DIST_2D; a0++) {
            int t = omp_get_thread_num();
            for (int a1 = 0; a1 < N_SUPPORT_DIST_2D; a1++) {
                int e0 = center_mod(a0, N_SUPPORT_DIST_2D);
                int e1 = center_mod(a1, N_SUPPORT_DIST_2D);

                point_t target = {{points[j].x[0] + e0, points[j].x[1] + e1}};
                int i = find_closest_codeword(points, &target, kyber_q);

                if (i != j) {
                    mpfr_add(pair_error_prob_for_thread[t],
                             pair_error_prob_for_thread[t],
                             delta_m->m[a0][a1],
                             MPFR_RNDD);
                }
            }
        }
    }
    // mpfr_t error_prob;
    // mpfr_init2(error_prob, N_PRECISION_BITS);
    mpfr_set_d(error_prob, 0, MPFR_RNDD);

    // Sums the partial results from each thread
    for (int t = 0; t < nthreads; t++) {
        mpfr_add(error_prob, error_prob, pair_error_prob_for_thread[t], MPFR_RNDD);
    }
    // Multiplies by probability of each codeword appearing (division by 4)
    mpfr_div_ui(error_prob, error_prob, N_CODEWORDS, MPFR_RNDD);
    // Union bound over pairs of coefficients (n/2)
    mpfr_mul_ui(error_prob, error_prob, 128, MPFR_RNDD);

    for (int t = 0; t < nthreads; t++) {
        mpfr_clear(pair_error_prob_for_thread[t]);
    }
}


static __inline__ int decode_msg_pair(int code_beta, int32_t x, int32_t y) {
  uint32_t reflect_mask = lower_than_mask(x, y);

  int32_t x_prime = (x & ~reflect_mask) | (y & reflect_mask);
  int32_t y_prime = (y & ~reflect_mask) | (x & reflect_mask);

  uint8_t above_l1 = ABOVE_L1(code_beta, x_prime, y_prime);
  uint8_t above_l2 = ABOVE_L2(code_beta, x_prime, y_prime);
  uint8_t above_l3 = ABOVE_L3(code_beta, x_prime, y_prime);
  uint8_t above_l4 = ABOVE_L4(code_beta, x_prime, y_prime);
  uint8_t above_l5 = ABOVE_L5(code_beta, x_prime, y_prime);

  // Not needed, but conceptually: uint8_t c00 = (~above_l3 & above_l2 & above_l4);
  uint8_t c01 = ~above_l2 & ~above_l5 & ~above_l3;
  uint8_t c10 = above_l2 & above_l3 & ~above_l1;
  uint8_t c11 = above_l1 | (above_l3 & ~above_l2) | (above_l5 & ~above_l4);

  c01 &= (1 ^ reflect_mask);
  c10 &= (2 ^ reflect_mask);

  uint8_t bits = c01 | c10 | c11;

  return bits & 3;
}

static __inline__ int decode_msg_pair_original_kyber(int x, int y) {
    int KYBER_Q = 3329;

    int t_x = x;
    t_x += ((int16_t)t_x >> 15) & KYBER_Q;
    t_x  = (((t_x << 1) + KYBER_Q/2)/KYBER_Q) & 1;

    int t_y = y;
    t_y += ((int16_t)t_y >> 15) & KYBER_Q;
    t_y  = (((t_y << 1) + KYBER_Q/2)/KYBER_Q) & 1;

    return (t_x << 1) | (t_y);
}


// This dfr computation uses our decoder, which is based on approximate Voronoi cells and exploit the
// geometrical properties of our codes.
// Notice that this is not generic: we need q = 3320 and alpha = 1664, because the decoding lines defined
// in `decoding_lines.h` are computed with these values.
void get_dfr_for_2d_code_with_real_decoder(mpfr_t error_prob, mpfr_matrix_t *delta_m, int code_alpha, int code_beta, int kyber_q) {
    // THIS FUNCTION IS ONLY IMPLEMENTED FOR Q = 3329!
    assert(kyber_q == 3329);
    assert(code_alpha == 1664);

    int nthreads = omp_get_max_threads();

    assert(0 <= code_alpha);
    assert(code_alpha < kyber_q);
    assert(0 <= code_beta);
    assert(code_beta < kyber_q);
    assert(code_beta + code_beta < kyber_q);

    point_t codewords[N_CODEWORDS] = {
        {{0, 0}},
        {{code_beta, code_alpha}},
        {{code_alpha, code_beta}},
        {{code_alpha + code_beta, code_alpha + code_beta}},
    };


    mpfr_t pair_error_prob_for_thread[nthreads];
    for (int t = 0; t < nthreads; t++) {
        mpfr_init2(pair_error_prob_for_thread[t], N_PRECISION_BITS);
        mpfr_set_d(pair_error_prob_for_thread[t], 0, MPFR_RNDD);
    }

    mpfr_t zero;
    mpfr_init2(zero, N_PRECISION_BITS);
    mpfr_set_d(zero, ZERO_DOUBLE, MPFR_RNDD);

    for (int j = 0; j < N_CODEWORDS; j++) {

        // Iterates over possible pairs of errors (e0, e1)
        #pragma omp parallel for schedule(static) // Using static scheduling the results are reproducible
        for (int a0 = 0; a0 < N_SUPPORT_DIST_2D; a0++) {
            int t = omp_get_thread_num();

            int e0 = center_mod(a0, N_SUPPORT_DIST_2D);
            int x = center_mod(codewords[j].x[0] + e0, kyber_q);
            #pragma GCC unroll 64
            for (int a1 = 0; a1 < N_SUPPORT_DIST_2D; a1++) {
                int cmp_real = mpfr_cmp_abs(delta_m->m[a0][a1], zero);
                if (cmp_real <= 0) {
                    continue;
                }

                int e1 = center_mod(a1, N_SUPPORT_DIST_2D);

                // Remember that `decode_msg_pair` takes the centralized representatives
                int y = center_mod(codewords[j].x[1] + e1, kyber_q);

                int decoded = -1;
                if (code_beta == 0) {
                    // Remember that the decoding lines used for `decode_msg_pair` are undefined for beta = 0.
                    // Therefore we use the original kyber decoding procedure in this case.
                    decoded = decode_msg_pair_original_kyber(x, y);
                }
                else {
                    decoded = decode_msg_pair(code_beta, x, y);
                }
                if (decoded != j) {
                    mpfr_add(pair_error_prob_for_thread[t],
                             pair_error_prob_for_thread[t],
                             delta_m->m[a0][a1],
                             MPFR_RNDD);
                }
            }
        }
    }
    mpfr_clear(zero);
    mpfr_set_d(error_prob, 0, MPFR_RNDD);

    // Sums the partial results from each thread
    for (int t = 0; t < nthreads; t++) {
        mpfr_add(error_prob, error_prob, pair_error_prob_for_thread[t], MPFR_RNDD);
    }
    // Multiplies by probability of each codeword appearing (division by 4)
    mpfr_div_ui(error_prob, error_prob, N_CODEWORDS, MPFR_RNDD);
    // Union bound over pairs of coefficients (n/2)
    mpfr_mul_ui(error_prob, error_prob, 128, MPFR_RNDD);

    for (int t = 0; t < nthreads; t++) {
        mpfr_clear(pair_error_prob_for_thread[t]);
    }
}

void find_best_code(mpfr_matrix_t *delta_m, int kyber_q) {

    for (int code_beta = 0; code_beta < 500; code_beta++) {

        mpfr_t error_prob;
        mpfr_init2(error_prob, N_PRECISION_BITS);
        mpfr_set_d(error_prob, 0, MPFR_RNDD);


        printf("%d,", code_beta);

        // You can uncomment the lines below to print the DFR computation using the minimum distance decoder.
        // get_dfr_for_2d_code_with_minimum_distance_decoding(error_prob, delta_m, kyber_q/2, code_beta, kyber_q);
        // mpfr_out_str(stdout, 10, 0, error_prob, MPFR_RNDD);
        // printf(",");

        get_dfr_for_2d_code_with_real_decoder(error_prob, delta_m, kyber_q/2, code_beta, kyber_q);
        mpfr_out_str(stdout, 10, 0, error_prob, MPFR_RNDD);
        printf("\n"); fflush(stdout);

        mpfr_clear(error_prob);
    }
}

#ifdef MAIN_PROGRAM

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s kyber_q delta_m.csv\n", argv[0]);
        exit(1);
    }
    int kyber_q = atoi(argv[1]);
    mpfr_matrix_t delta_m;
    printf("Reading delta_m file... "); fflush(stdout);
    mpfr_matrix_init(&delta_m, N_SUPPORT_DIST_2D, N_SUPPORT_DIST_2D);
    mpfr_matrix_init_from_csv(&delta_m, argv[2]);
    printf("Done.\n"); fflush(stdout);
    find_best_code(&delta_m, kyber_q);
    mpfr_matrix_free(&delta_m);
}

#endif
