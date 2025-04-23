// Copyright 2025 LG Electronics, Inc. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <mpfr.h>
#include <sys/resource.h>
#include <omp.h>

#include "mp_matrix.h"
#include "fft_mpc.h"
#include "dfr_2d_codes.h"
#include "prob_dist.h"



#ifndef SAVE_PARTIAL_DISTRIBUTION_RESULTS
#define SAVE_PARTIAL_DISTRIBUTION_RESULTS 0
#endif

#define MOD(a, b) ((((a) % (b)) + (b)) % (b))

#define MAX_STR 1000

// #define EXPECTED_SETUP_CSV_HEADER "KYBER_Q,KYBER_LEVEL,KYBER_K,KYBER_ETA1,KYBER_ETA2,KYBER_N_BLOCKS_FOR_DU0,KYBER_DU0,KYBER_DU1,KYBER_DV,CODE_DIMENSION\n"
#define EXPECTED_SETUP_CSV_HEADER "KYBER_Q,KYBER_LEVEL,KYBER_K,KYBER_ETA1,KYBER_ETA2,KYBER_N_BLOCKS_FOR_DU0,KYBER_DU0,KYBER_DU1,KYBER_DV,CODE_BETA,CODE_PNORM\n"

#define KYBER_LOG2_OF_N 8 // KYBER_N = 2**KYBER_LOG2_OF_N

// #include "fish_codes/all_codes.h"
// #include "minal4d_defs/dfr_relevant_points.h"
// #include "fish_codes/weight3_code.h"

#include "minal4d_defs/dfr_relevant_points_for_selected_betas.h"

// #define CODE_BETA 746
// #define CODE_BETA 740
// #define CODE_BETA 717
//
// #define CODE_IDX (CODE_BETA - STARTING_BETA)


#define MAX_DIFFERENT_CODEWORD_VALUES 10000

typedef struct precomputed_fft_dv_eta2_dot_codeword_s {
    int n_codeword_values;
    int codeword_values[MAX_DIFFERENT_CODEWORD_VALUES];
    int codeword_value_to_index_map;
    mpc_matrix_t *table;

} precomputed_fft_dv_eta2_dot_codeword_t;

typedef struct kyber_parameters_s {
    int KYBER_Q;
    int KYBER_LEVEL;
    int KYBER_K;
    int KYBER_ETA1;
    int KYBER_ETA2;
    // We use a generalized view of `KYBER_DU` as a vector of `KYBER_K` entries of the form:
    // [KYBER_DU0] * KYBER_N_BLOCKS_FOR_DU0 + [KYBER_DU1] * (KYBER_K - KYBER_N_BLOCKS_FOR_DU0)
    int KYBER_N_BLOCKS_FOR_DU0;
    int KYBER_DU0;
    int KYBER_DU1;
    int KYBER_DV;
    int CODE_BETA;
    int CODE_DIMENSION;

    // This table contains the precomputed probability distributions of the product
    // `(delta_v + eta_2) * value` in FFT domain for all values `value` that appear in
    // the DFR relevant vectors.
    precomputed_fft_dv_eta2_dot_codeword_t precomputed_fft_dv_eta2_dot_codeword;
} kyber_parameters_t;


void print_dfr_csv_for_2d_codes_header();
int find_value_in_list(int value, int list[], int list_size);

static int file_exists(const char filename[]) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        return 0;
    }
    fclose(file);
    return 1;
}

void print_mpc_matrix_to_file(const char filename[], mpc_matrix_t *dist) {
    assert(!file_exists(filename));

    FILE *f = fopen(filename, "w");
    if (!f) {
        fprintf(stderr, "Could not open output file %s\n", filename);
        exit(1);
    }
    mpc_matrix_print_to_file_ptr(f, dist);
    fclose(f);
}


// Returns 1 if updated, and 0 if last vector reached
int update_vector_to_next_one(prob_dist_1d_t *coeff_distribution, int current_vector[], int n_dimensions) {
    for (int i = n_dimensions - 1; i >= 0; i--) {
        int v = current_vector[i];
        int ret = prob_dist_1d_update_value_to_next_one(coeff_distribution, &v);
        if (ret == 1) {
            current_vector[i] = v;
            return 1;
        }
        prob_dist_1d_update_value_to_next_one(coeff_distribution, &v);
        current_vector[i] = v;
    }
    return 0;
}

int get_first_nonnegligible_vector(prob_dist_1d_t *coeff_distribution, int current_vector[], int n_dimensions) {
    for (int i = 0; i < n_dimensions; i++) {
        current_vector[i] = coeff_distribution->min_value - 1;
        prob_dist_1d_update_value_to_next_one(coeff_distribution, &current_vector[i]);
    }
    return 1;
}


// Extracted from: https://github.com/mkannwischer/polymul/blob/master/C/01schoolbook.c
static void polymul_schoolbook(int32_t* c, int32_t* a, size_t aN, int32_t* b, size_t bN){
    size_t i,j;
    uint32_t t;
    for(i=0; i<aN+bN-1; i++){
        c[i] = 0;
    }

    for(i=0;i<aN;i++){
        for(j=0;j<bN;j++){
            // since Q is very small, we could compute a 16 bits product only
            t = ((uint32_t)a[i] * b[j]);
            c[i+j] = (c[i+j] + t);
        }
    }
}

// Extracted from: https://github.com/mkannwischer/polymul/blob/master/C/01schoolbook.c
static void polymul_schoolbook_negacyclic(int32_t* c, int32_t* a, int32_t* b, size_t n, int32_t q) {
    int32_t t[2*n-1];
    size_t i;
    // perform double-sized schoolbook
    polymul_schoolbook(t, a, n, b, n);

    // reduce mod (X^n+1, q)
    for(i=0;i<n;i++){
        c[i] = t[i];
    }
    for(i=0;i<n-1;i++){
        c[i] = center_mod(c[i] - t[i+n], q);
    }
}


void print_vector(int vector[], int n) {
    fprintf(stderr, "[ ");
    for (int i = 0; i < n; i++) {
        fprintf(stderr, "%2d ", vector[i]);
    }
    fprintf(stderr, "]");
}

void vector_get_mpc_probability(mpc_t prob, prob_dist_1d_t *coeff_dist, int vector[], int n) {
    mpc_set_d(prob, 1, ROUND_DOWN);
    for (int i = 0; i < n; i++) {
        int idx = prob_dist_1d_get_idx_for_value(coeff_dist, vector[i]);
        mpc_mul_fr(prob, prob, coeff_dist->d[idx], MPFR_RNDD);
    }
}

void vector_get_mpfr_probability(mpfr_t prob, prob_dist_1d_t *coeff_dist, int vector[], int n) {
    mpfr_set_d(prob, 1, ROUND_DOWN);
    for (int i = 0; i < n; i++) {
        int idx = prob_dist_1d_get_idx_for_value(coeff_dist, vector[i]);
        mpfr_mul(prob, prob, coeff_dist->d[idx], MPFR_RNDD);
    }
}


void get_base_dot_product_distribution(mpc_matrix_t *base_dist,
                                                    prob_dist_1d_t *coeff_dist_a,
                                                    prob_dist_1d_t *coeff_dist_b,
                                                    int voronoi_relevant_vector[],
                                                    int n_dims,
                                                    int kyber_q) {
    assert(prob_dist_1d_is_symmetric(coeff_dist_b));

    mpc_matrix_init(base_dist, 1, MAX_N_FOR_FFT); // TODO: Have to be very careful with this bound

    mpfr_t tmp_prob;
    mpfr_init2(tmp_prob, N_PRECISION_BITS);

    int vector_a[n_dims];
    int ret_a = get_first_nonnegligible_vector(coeff_dist_a, vector_a, n_dims);

    mpc_t prob_a;
    mpc_init2(prob_a, N_PRECISION_BITS);
    mpc_t prob_b;
    mpc_init2(prob_b, N_PRECISION_BITS);
    mpc_t prob_ab;
    mpc_init2(prob_ab, N_PRECISION_BITS);
    while (ret_a == 1) {

        vector_get_mpc_probability(prob_a, coeff_dist_a, vector_a, n_dims);
        int vector_b[n_dims];
        int ret_b = get_first_nonnegligible_vector(coeff_dist_b, vector_b, n_dims);
        while (ret_b == 1) {
            vector_get_mpc_probability(prob_b, coeff_dist_b, vector_b, n_dims);

            mpc_mul(prob_ab, prob_a, prob_b, MPFR_RNDD);


            int vector_product[n_dims];
            polymul_schoolbook_negacyclic(vector_product, vector_a, vector_b, n_dims, kyber_q);

            int dot_product = 0;
            for (int i = 0; i < n_dims; i++)
                dot_product += vector_product[i] * voronoi_relevant_vector[i];

            // print_vector(vector_a, n_dims);
            // print_vector(vector_b, n_dims);
            // printf(" = ");
            // print_vector(vector_product, n_dims);
            // printf("| ");
            // mpc_out_str(stdout, 10, 0, prob_ab, MPFR_RNDD);
            // printf(" | dot_product = %d\n", dot_product);

            // TODO: assert(dot_product)
            int idx = MOD(dot_product, base_dist->n_cols);
            mpc_add(base_dist->m[0][idx], base_dist->m[0][idx], prob_ab, MPFR_RNDD);
            ret_b = update_vector_to_next_one(coeff_dist_b, vector_b, n_dims);
        }
        ret_a = update_vector_to_next_one(coeff_dist_a, vector_a, n_dims);
    }

    mpc_clear(prob_a);
    mpc_clear(prob_b);
    mpc_clear(prob_ab);
    mpfr_clear(tmp_prob);
}


void self_convolution_in_fft_domain(mpc_matrix_t *dist, unsigned long order) {
    for (int i = 0; i < dist->n_rows; i++) {
        #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < dist->n_cols; j++) {
            mpc_pow_ui(dist->m[i][j], dist->m[i][j], order, MPFR_RNDU);
        }
    }
}


void init_as_convolution_in_fft_domain(mpc_matrix_t *out, mpc_matrix_t *dist1, mpc_matrix_t *dist2) {
    assert(dist1->n_cols == dist2->n_cols);
    assert(dist1->n_rows == dist2->n_rows);
    mpc_matrix_init(out, dist1->n_rows, dist1->n_cols);

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < out->n_rows; i++) {
        for (int j = 0; j < out->n_cols; j++) {
            mpc_mul(out->m[i][j], dist1->m[i][j], dist2->m[i][j], MPFR_RNDU);
        }
    }
}

void init_as_copy(mpc_matrix_t *out, mpc_matrix_t *dist) {
    mpc_matrix_init(out, dist->n_rows, dist->n_cols);

    for (int i = 0; i < out->n_rows; i++) {
        for (int j = 0; j < out->n_cols; j++) {
            mpc_set(out->m[i][j], dist->m[i][j], MPFR_RNDU);
        }
    }
}

void print_params(const char filename[], kyber_parameters_t *params) {
    assert(!file_exists(filename));

    FILE *file = fopen(filename, "w");
    if (!file) {
        fprintf(stderr, "Could not open %s\n", filename);
        exit(1);
    }

    fprintf(file, "parameter,value\n");
    fprintf(file, "KYBER_Q,%d\n", params->KYBER_Q);
    fprintf(file, "KYBER_LOG2_OF_N,%d\n", KYBER_LOG2_OF_N);
    fprintf(file, "KYBER_LEVEL,%d\n", params->KYBER_LEVEL);
    fprintf(file, "KYBER_K,%d\n", params->KYBER_K);
    fprintf(file, "KYBER_ETA1,%d\n", params->KYBER_ETA1);
    fprintf(file, "KYBER_ETA2,%d\n", params->KYBER_ETA2);
    fprintf(file, "KYBER_N_BLOCKS_FOR_DU0,%d\n", params->KYBER_N_BLOCKS_FOR_DU0);
    fprintf(file, "KYBER_DU0,%d\n", params->KYBER_DU0);
    fprintf(file, "KYBER_DU1,%d\n", params->KYBER_DU1);
    fprintf(file, "KYBER_DV,%d\n", params->KYBER_DV);
    fprintf(file, "N_PRECISION_BITS,%d\n", N_PRECISION_BITS);
    fprintf(file, "ZERO_DOUBLE,%.30g\n", ZERO_DOUBLE);

    fclose(file);
}

#if SAVE_PARTIAL_DISTRIBUTION_RESULTS != 1
void DEBUG_print_params(__attribute__((unused)) const char output_dir[],
                        __attribute__((unused)) kyber_parameters_t *params,
                        __attribute__((unused)) const char base_name[]) {
    return;
}

void DEBUG_print_dist_1d_to_file(__attribute__((unused)) const char output_dir[],
                                 __attribute__((unused)) kyber_parameters_t *params,
                                 __attribute__((unused)) const char base_name[],
                                 __attribute__((unused)) prob_dist_1d_t *dist) {

    return;
}

void DEBUG_print_mpc_matrix_to_file(__attribute__((unused)) const char output_dir[],
                                    __attribute__((unused)) kyber_parameters_t *params,
                                    __attribute__((unused)) const char base_name[],
                                    __attribute__((unused)) mpc_matrix_t *dist) {

    return;
}

#else

void DEBUG_print_params(const char output_dir[],
                        kyber_parameters_t *params,
                        const char base_name[]) {
    char filename[MAX_STR];
    sprintf(filename,
            "%s/%s_q=%d_l=%d_k=%d_e1=%d_e2=%d_du0rep=%d_du0=%d_du1=%d_dv=%d.csv", output_dir, base_name,
            params->KYBER_Q, params->KYBER_LEVEL, params->KYBER_K, params->KYBER_ETA1, params->KYBER_ETA2,
            params->KYBER_N_BLOCKS_FOR_DU0, params->KYBER_DU0, params->KYBER_DU1, params->KYBER_DV);
    print_params(filename, params);
}

void DEBUG_print_dist_1d_to_file(const char output_dir[],
                                 kyber_parameters_t *params,
                                 const char base_name[],
                                 prob_dist_1d_t *dist) {

    char filename[MAX_STR];
    sprintf(filename,
            "%s/%s_q=%d_l=%d_k=%d_e1=%d_e2=%d_du0rep=%d_du0=%d_du1=%d_dv=%d.csv", output_dir, base_name,
            params->KYBER_Q, params->KYBER_LEVEL, params->KYBER_K, params->KYBER_ETA1, params->KYBER_ETA2,
            params->KYBER_N_BLOCKS_FOR_DU0, params->KYBER_DU0, params->KYBER_DU1, params->KYBER_DV);
    print_dist_1d_to_file(filename, dist);
}

void DEBUG_print_mpc_matrix_to_file(const char output_dir[],
                                    kyber_parameters_t *params,
                                    const char base_name[],
                                    mpc_matrix_t *dist) {
    char filename[MAX_STR];
    sprintf(filename,
            "%s/%s_q=%d_l=%d_k=%d_e1=%d_e2=%d_du0rep=%d_du0=%d_du1=%d_dv=%d.csv", output_dir, base_name,
            params->KYBER_Q, params->KYBER_LEVEL, params->KYBER_K, params->KYBER_ETA1, params->KYBER_ETA2,
            params->KYBER_N_BLOCKS_FOR_DU0, params->KYBER_DU0, params->KYBER_DU1, params->KYBER_DV);
    print_mpc_matrix_to_file(filename, dist);
}

#endif


void init_as_2d_from_1d(mpc_matrix_t *out, prob_dist_1d_t *dist) {
    mpc_matrix_init(out, N_SUPPORT_DIST_2D, N_SUPPORT_DIST_2D);

    mpfr_t tmp_prob;
    mpfr_init2(tmp_prob, N_PRECISION_BITS);
    for (int a0 = dist->min_value; a0 <= dist->max_value; a0++) {
        int idx_a0 = prob_dist_1d_get_idx_for_value(dist, a0);

        for (int a1 = dist->min_value; a1 <= dist->max_value; a1++) {
            int idx_a1 = prob_dist_1d_get_idx_for_value(dist, a1);

            mpfr_mul(tmp_prob, dist->d[idx_a0], dist->d[idx_a1], MPFR_RNDD);
            mpc_set_fr(out->m[MOD(a0, N_SUPPORT_DIST_2D)][MOD(a1, N_SUPPORT_DIST_2D)], tmp_prob, ROUND_DOWN);
        }
    }
    mpfr_clear(tmp_prob);
}


void init_as_dv_eta2_error_slow_for_testing(prob_dist_1d_t *out, prob_dist_1d_t *coeff_dist, int voronoi_relevant_vector[], int n_dims) {

    // mpc_matrix_init(out, 1, MAX_N_FOR_FFT);
    prob_dist_1d_init(out, -(1 << 20), (1 << 20));

    int vector_a[n_dims];
    int ret_a = get_first_nonnegligible_vector(coeff_dist, vector_a, n_dims);

    // mpc_t prob_a;
    // mpc_init2(prob_a, N_PRECISION_BITS);

    mpfr_t prob_a;
    mpfr_init2(prob_a, N_PRECISION_BITS);

    // int max_dot_product = 0;
    // int max_vec[n_dims];
    fprintf(stderr, "entering while\n"); fflush(stderr);
    while (ret_a == 1) {

        vector_get_mpfr_probability(prob_a, coeff_dist, vector_a, n_dims);

        int dot_product = 0;
        for (int i = 0; i < n_dims; i++)
            dot_product += vector_a[i] * voronoi_relevant_vector[i];

        int idx = prob_dist_1d_get_idx_for_value(out, dot_product);

        mpfr_add(out->d[idx], out->d[idx], prob_a, MPFR_RNDD);

        ret_a = update_vector_to_next_one(coeff_dist, vector_a, n_dims);
    }
    mpfr_clear(prob_a);
    fprintf(stderr, "finished while\n"); fflush(stderr);

}

void init_prob_dist_1d_from_mpc_matrix_1d(prob_dist_1d_t *dist, mpc_matrix_t *matrix);

void init_as_dv_eta2_error(prob_dist_1d_t *out, prob_dist_1d_t *coeff_dist, int voronoi_relevant_vector[], int n_dims) {

    mpc_matrix_t out_as_matrix;
    mpc_matrix_init(&out_as_matrix, 1, MAX_N_FOR_FFT);

    fprintf(stderr, "Sec A\n"); fflush(stderr);

    for (int i = 0; i < n_dims; i++) {
        mpc_matrix_t tmp;
        mpc_matrix_init(&tmp, 1, MAX_N_FOR_FFT);

        for (int v = coeff_dist->min_value; v <= coeff_dist->max_value; v++) {
            int idx = prob_dist_1d_get_idx_for_value(coeff_dist, v);

            int value = v * voronoi_relevant_vector[i];
            mpc_add_fr(tmp.m[0][MOD(value, tmp.n_cols)], tmp.m[0][MOD(value, tmp.n_cols)],
                       coeff_dist->d[idx], MPFR_RNDD);
        }
        fft(tmp.m[0], tmp.n_cols);

        if (i == 0) {
            for (int j = 0; j < tmp.n_cols; j++) {
                mpc_set(out_as_matrix.m[0][j], tmp.m[0][j], MPFR_RNDD);
            }
        }
        else {
            for (int j = 0; j < tmp.n_cols; j++) {
                mpc_mul(out_as_matrix.m[0][j], out_as_matrix.m[0][j], tmp.m[0][j], MPFR_RNDD);
            }
        }
        mpc_matrix_free(&tmp);
    }
    fprintf(stderr, "Sec B\n"); fflush(stderr);

    ifft(out_as_matrix.m[0], out_as_matrix.n_cols);
    fprintf(stderr, "Sec C\n"); fflush(stderr);

    init_prob_dist_1d_from_mpc_matrix_1d(out, &out_as_matrix);
    fprintf(stderr, "Sec D\n"); fflush(stderr);

    mpc_matrix_free(&out_as_matrix);
}


void init_as_dv_eta2_error_fast(prob_dist_1d_t *out, kyber_parameters_t *params, int voronoi_relevant_vector[], int n_dims) {

    mpc_matrix_t out_as_matrix;
    mpc_matrix_init(&out_as_matrix, 1, MAX_N_FOR_FFT);

    precomputed_fft_dv_eta2_dot_codeword_t *fft_table = &(params->precomputed_fft_dv_eta2_dot_codeword);


    fprintf(stderr, "Sec A\n"); fflush(stderr);
    for (int i = 0; i < n_dims; i++) {
        int table_idx = find_value_in_list(voronoi_relevant_vector[i], fft_table->codeword_values, fft_table->n_codeword_values);
        assert(table_idx != -1);
        mpc_matrix_t *tmp = &(fft_table->table[table_idx]);
        if (i == 0) {
            // #pragma omp parallel for
            for (int j = 0; j < tmp->n_cols; j++) {
                mpc_set(out_as_matrix.m[0][j], tmp->m[0][j], MPFR_RNDD);
            }
        }
        else {
            // #pragma omp parallel for
            for (int j = 0; j < tmp->n_cols; j++) {
                mpc_mul(out_as_matrix.m[0][j], out_as_matrix.m[0][j], tmp->m[0][j], MPFR_RNDD);
            }
        }
    }
    fprintf(stderr, "Sec B\n"); fflush(stderr);

    ifft(out_as_matrix.m[0], out_as_matrix.n_cols);
    fprintf(stderr, "Sec C\n"); fflush(stderr);

    init_prob_dist_1d_from_mpc_matrix_1d(out, &out_as_matrix);
    fprintf(stderr, "Sec D\n"); fflush(stderr);

    mpc_matrix_free(&out_as_matrix);
}

void init_as_dv_eta2_error_fast_fft_domain(mpc_matrix_t *out_as_matrix, kyber_parameters_t *params, int voronoi_relevant_vector[], int n_dims) {

    precomputed_fft_dv_eta2_dot_codeword_t *fft_table = &(params->precomputed_fft_dv_eta2_dot_codeword);

    fprintf(stderr, "Sec A\n"); fflush(stderr);

    #pragma omp parallel for
    for (int j = 0; j < out_as_matrix->n_cols; j++) {
        mpc_set_d(out_as_matrix->m[0][j], 1, MPFR_RNDD);
    }

    for (int i = 0; i < n_dims; i++) {
        if (voronoi_relevant_vector[i] == 0) continue;

        int table_idx = find_value_in_list(abs(voronoi_relevant_vector[i]), fft_table->codeword_values, fft_table->n_codeword_values);
        int sign = voronoi_relevant_vector[i] >= 0 ? 1 : -1;
        assert(table_idx != -1);

        mpc_matrix_t *tmp = &(fft_table->table[table_idx]);
        #pragma omp parallel for
        for (int j = 0; j < tmp->n_cols; j++) {
        int idx = MOD(sign * j, tmp->n_cols);
            mpc_mul(out_as_matrix->m[0][idx], out_as_matrix->m[0][idx], tmp->m[0][idx], MPFR_RNDD);
        }
    }
}

void init_as_dv_eta2_error_fft_domain(mpc_matrix_t *out_as_matrix, prob_dist_1d_t *coeff_dist, int voronoi_relevant_vector[], int n_dims) {

    fprintf(stderr, "Sec A\n"); fflush(stderr);

    #pragma omp parallel for
    for (int j = 0; j < out_as_matrix->n_cols; j++) {
        mpc_set_d(out_as_matrix->m[0][j], 1, MPFR_RNDD);
    }

    for (int i = 0; i < n_dims; i++) {
        if (voronoi_relevant_vector[i] == 0) continue;

        mpc_matrix_t tmp;
        mpc_matrix_init(&tmp, 1, MAX_N_FOR_FFT);

        #pragma omp parallel for
        for (int v = coeff_dist->min_value; v <= coeff_dist->max_value; v++) {
            int idx = prob_dist_1d_get_idx_for_value(coeff_dist, v);

            int value = v * voronoi_relevant_vector[i];
            mpc_add_fr(tmp.m[0][MOD(value, tmp.n_cols)], tmp.m[0][MOD(value, tmp.n_cols)],
                       coeff_dist->d[idx], MPFR_RNDD);
        }
        fft(tmp.m[0], tmp.n_cols);

        #pragma omp parallel for
        for (int j = 0; j < tmp.n_cols; j++) {
            mpc_mul(out_as_matrix->m[0][j], out_as_matrix->m[0][j], tmp.m[0][j], MPFR_RNDD);
        }
        mpc_matrix_free(&tmp);
    }
    fprintf(stderr, "Sec B\n"); fflush(stderr);
}

void init_prob_dist_1d_from_mpc_matrix_1d(prob_dist_1d_t *dist, mpc_matrix_t *matrix) {
    assert(matrix->n_rows == 1);

    prob_dist_1d_init(dist, -matrix->n_cols/2 - 1, matrix->n_cols/2 + 1);

    for (int i = 0; i < matrix->n_cols; i++) {
        int value = center_mod(i, matrix->n_cols);
        int idx = prob_dist_1d_get_idx_for_value(dist, value);
        mpc_real(dist->d[idx], matrix->m[0][i], ROUND_DOWN);
    }
}

double norm_squared(int vector[], int n_dims) {
    double norm_sqr = 0;
    for (int i = 0; i < n_dims; i++) {
        norm_sqr += vector[i] * vector[i];
    }
    return norm_sqr;
}



void mpc_matrix_get_probability_of_greater_than_value(mpfr_t out, mpc_matrix_t *matrix, double tgt_value) {
    mpfr_set_d(out, 0, MPFR_RNDD);

    mpfr_t tmp_real;
    mpfr_init2(tmp_real, N_PRECISION_BITS);

    assert(matrix->n_rows == 1);
    for (int i = 0; i < matrix->n_cols; i++) {
        int value = center_mod(i, matrix->n_cols);
        if (value >= tgt_value) {
            mpc_real(tmp_real, matrix->m[0][i], ROUND_DOWN);
            mpfr_add(out, out, tmp_real, MPFR_RNDD);
        }
    }

    mpfr_clear(tmp_real);
}

int compute_error_probability_for_a_voronoi_relevant_vector(mpfr_t upper_bound_for_one_vector,
                                                            kyber_parameters_t *params,
                                                            int voronoi_relevant_vector[]) {
    assert(params->KYBER_N_BLOCKS_FOR_DU0 == params->KYBER_K);
    // DEBUG_print_params(output_dir, params, "params.csv");

    prob_dist_1d_t dist_dv;
    prob_dist_1d_init_as_compress_decompress_error(&dist_dv, params->KYBER_DV, params->KYBER_Q);

    prob_dist_1d_t dist_eta1;
    prob_dist_1d_init_as_centered_binomial(&dist_eta1, params->KYBER_ETA1);

    prob_dist_1d_t dist_eta2;
    prob_dist_1d_init_as_centered_binomial(&dist_eta2, params->KYBER_ETA2);

    prob_dist_1d_t dist_du0;
    prob_dist_1d_init_as_compress_decompress_error(&dist_du0, params->KYBER_DU0, params->KYBER_Q);
    // print_dist_1d(stdout, &dist_du0);

    prob_dist_1d_t dist_du0_sum_eta2;
    prob_dist_1d_init_as_sum(&dist_du0_sum_eta2, &dist_du0, &dist_eta2);

    prob_dist_1d_t dist_dv_sum_eta2;
    prob_dist_1d_init_as_sum(&dist_dv_sum_eta2, &dist_dv, &dist_eta2);

    // Base distribution `dist_s_du_sum_eta2_product`
    mpc_matrix_t dist_s_du_sum_eta2_product;
    get_base_dot_product_distribution(&dist_s_du_sum_eta2_product, &dist_du0_sum_eta2, &dist_eta1, voronoi_relevant_vector,
                                      params->CODE_DIMENSION, params->KYBER_Q);
    fft(dist_s_du_sum_eta2_product.m[0], dist_s_du_sum_eta2_product.n_cols);

    int conv_order = params->KYBER_K * 256 / (params->CODE_DIMENSION);
    self_convolution_in_fft_domain(&dist_s_du_sum_eta2_product, conv_order);

// Compile with -DSIMULATE_FOR_PK_COMPRESSION_11BITS to compute the DFR when 11 bits are used for the
// public key (instead of the regular 12 bits)
#ifdef SIMULATE_FOR_PK_COMPRESSION_11BITS
    prob_dist_1d_t dist_dt;
    prob_dist_1d_init_as_compress_decompress_error(&dist_dt, 11, params->KYBER_Q);
    prob_dist_1d_t dist_dt_sum_eta1;
    prob_dist_1d_init_as_sum(&dist_dt_sum_eta1, &dist_dt, &dist_eta1);

    // Base distribution `dist_e_r_product`
    mpc_matrix_t dist_e_r_product;
    get_base_dot_product_distribution(&dist_e_r_product, &dist_dt_sum_eta1, &dist_eta1, voronoi_relevant_vector,
                                      params->CODE_DIMENSION, params->KYBER_Q);

    prob_dist_1d_clear(&dist_dt);
    prob_dist_1d_clear(&dist_dt_sum_eta1);
#else
    // Base distribution `dist_e_r_product`
    mpc_matrix_t dist_e_r_product;
    get_base_dot_product_distribution(&dist_e_r_product, &dist_eta1, &dist_eta1, voronoi_relevant_vector,
                                      params->CODE_DIMENSION, params->KYBER_Q);
#endif

    // char filename[MAX_STR];
    // sprintf(filename, "tmp/dist_e_r_product%d.txt", I_CODEWORD);
    // print_mpc_matrix_to_file(filename, &dist_e_r_product);

    // return 1;

    fft(dist_e_r_product.m[0], dist_e_r_product.n_cols);
    self_convolution_in_fft_domain(&dist_e_r_product, conv_order);

    fprintf(stderr, "doing dist_dv_sum_eta2_dot_prod\n"); fflush(stderr);
    mpc_matrix_t dist_dv_sum_eta2_dot_prod;
    mpc_matrix_init(&dist_dv_sum_eta2_dot_prod, 1, MAX_N_FOR_FFT);

#if USE_PRECOMPUTED_FFT_DV_ETA2_DOT_CODEWORD == 0
    init_as_dv_eta2_error_fft_domain(&dist_dv_sum_eta2_dot_prod, &dist_dv_sum_eta2, voronoi_relevant_vector, params->CODE_DIMENSION);
#else
    init_as_dv_eta2_error_fast_fft_domain(&dist_dv_sum_eta2_dot_prod, params, voronoi_relevant_vector, params->CODE_DIMENSION);
#endif

    mpc_matrix_t dist_sum_dot_products;
    init_as_convolution_in_fft_domain(&dist_sum_dot_products, &dist_e_r_product, &dist_s_du_sum_eta2_product);

    mpc_matrix_t dist_dot_product_with_dfr_relevant_point;
    init_as_convolution_in_fft_domain(&dist_dot_product_with_dfr_relevant_point, &dist_sum_dot_products, &dist_dv_sum_eta2_dot_prod);

    ifft(dist_dot_product_with_dfr_relevant_point.m[0], dist_dot_product_with_dfr_relevant_point.n_cols);

    fprintf(stderr, "done\n"); fflush(stderr);

    double half_codeword_norm_squared = norm_squared(voronoi_relevant_vector, params->CODE_DIMENSION) / 2;

    // char filename[MAX_STR];
    // sprintf(filename, "tmp/dist_dot_product_with_dfr_relevant_point%d.txt", I_CODEWORD);
    // print_mpc_matrix_to_file(filename, &dist_dot_product_with_dfr_relevant_point);

    // return 1;
    // print_mpc_matrix_to_file("tmp/xxx_discard.txt", &dist_dot_product_with_dfr_relevant_point);
    mpc_matrix_get_probability_of_greater_than_value(upper_bound_for_one_vector,
                                                     &dist_dot_product_with_dfr_relevant_point,
                                                     half_codeword_norm_squared);



    print_vector(voronoi_relevant_vector, params->CODE_DIMENSION);
    fprintf(stderr, ": ");
    fprintf(stderr, "Failure probability with respect to this point is: ");
    mpfr_out_str(stderr, 10, 0, upper_bound_for_one_vector, ROUND_DOWN);
    fprintf(stderr, "\n");
    fflush(stderr);

    prob_dist_1d_clear(&dist_du0);
    prob_dist_1d_clear(&dist_dv);
    prob_dist_1d_clear(&dist_eta1);
    prob_dist_1d_clear(&dist_eta2);
    prob_dist_1d_clear(&dist_du0_sum_eta2);
    prob_dist_1d_clear(&dist_dv_sum_eta2);

    mpc_matrix_free(&dist_dv_sum_eta2_dot_prod);
    mpc_matrix_free(&dist_sum_dot_products);
    mpc_matrix_free(&dist_e_r_product);
    mpc_matrix_free(&dist_s_du_sum_eta2_product);
    mpc_matrix_free(&dist_dot_product_with_dfr_relevant_point);

    return 0;
}

int compute_error_probability_for_a_voronoi_relevant_vector_chernoff(kyber_parameters_t *params,
                                                            int voronoi_relevant_vector[]) {
    assert(params->KYBER_N_BLOCKS_FOR_DU0 == params->KYBER_K);
    // DEBUG_print_params(output_dir, params, "params.csv");

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

    // Base distribution `dist_s_du_sum_eta2_product`
    mpc_matrix_t dist_s_du_sum_eta2_product;
    get_base_dot_product_distribution(&dist_s_du_sum_eta2_product, &dist_du0_sum_eta2, &dist_eta1, voronoi_relevant_vector,
                                      params->CODE_DIMENSION, params->KYBER_Q);
    fft(dist_s_du_sum_eta2_product.m[0], dist_s_du_sum_eta2_product.n_cols);

    int conv_order = 1;
    self_convolution_in_fft_domain(&dist_s_du_sum_eta2_product, conv_order);

    // Base distribution `dist_e_r_product`
    mpc_matrix_t dist_e_r_product;
    get_base_dot_product_distribution(&dist_e_r_product, &dist_eta1, &dist_eta1, voronoi_relevant_vector,
                                      params->CODE_DIMENSION, params->KYBER_Q);

    fft(dist_e_r_product.m[0], dist_e_r_product.n_cols);
    self_convolution_in_fft_domain(&dist_e_r_product, conv_order);


    mpc_matrix_t dist_sum_dot_products;
    init_as_convolution_in_fft_domain(&dist_sum_dot_products, &dist_e_r_product, &dist_s_du_sum_eta2_product);
    ifft(dist_sum_dot_products.m[0], dist_sum_dot_products.n_cols);

    prob_dist_1d_t dist_sum_dot_products_1d;
    init_prob_dist_1d_from_mpc_matrix_1d(&dist_sum_dot_products_1d, &dist_sum_dot_products);

    // print_dist_1d(stdout, &dist_sum_dot_products_1d);

    fprintf(stderr, "doing dist_dv_sum_eta2_dot_prod\n"); fflush(stderr);
    prob_dist_1d_t dist_dv_sum_eta2_dot_prod;
    // init_as_dv_eta2_error(&dist_dv_sum_eta2_dot_prod, &dist_dv_sum_eta2, voronoi_relevant_vector, params->CODE_DIMENSION);
    init_as_dv_eta2_error_fast(&dist_dv_sum_eta2_dot_prod, params, voronoi_relevant_vector, params->CODE_DIMENSION);

    fprintf(stderr, "done\n"); fflush(stderr);

    // Precompute Chernoff aux values
    double half_codeword_norm_squared = 0;
    int n_dims = params->CODE_DIMENSION;
    int self_convolution_order = params->KYBER_K * 256 / (n_dims * conv_order);

    for (int i = 0; i < n_dims; i++) {
        half_codeword_norm_squared += voronoi_relevant_vector[i] * voronoi_relevant_vector[i];
    }
    half_codeword_norm_squared /= 2;

    mpfr_t tmp_chernoff_bound2;
    mpfr_init2(tmp_chernoff_bound2, N_PRECISION_BITS);

    fprintf(stderr, "A start\n");
    fflush(stderr);

    double best_t = prob_dist_1d_minimize_t_for_chernoff(half_codeword_norm_squared, self_convolution_order, &dist_sum_dot_products_1d);
    // double best_t = 0.000113;
    // double best_t = 0.00020;
    fprintf(stderr, "A min\n");
    fflush(stderr);

    // double best_t = 0.000177;
    mpfr_t expectation_of_exp_tX;
    mpfr_init2(expectation_of_exp_tX, N_PRECISION_BITS);
    prob_dist_1d_get_expectation_of_exp_tX(expectation_of_exp_tX, &dist_sum_dot_products_1d, best_t);

    fprintf(stderr, "A part\n");
    fflush(stderr);

    prob_dist_1d_chernoff_with_precomputed_expectation_of_exp_tX(tmp_chernoff_bound2, half_codeword_norm_squared, self_convolution_order, best_t,
                                                                 expectation_of_exp_tX);
    fprintf(stderr, "A Done\n");
    fflush(stderr);

    fprintf(stderr, "Codeword: ");
    for (int i = 0; i < n_dims; i++) {
        fprintf(stderr, "%d, ", voronoi_relevant_vector[i]);
    }
    fprintf(stderr, "\n");
    fprintf(stderr, "  best_t = %lf\n", best_t);

    fprintf(stderr, "  Pr(sum_{i=1}^{self_convolution_order} X_i >= %lf) <= ", half_codeword_norm_squared);
    mpfr_out_str(stdout, 10, 0, tmp_chernoff_bound2, ROUND_DOWN);
    fprintf(stderr, "\n");
    fflush(stderr);

    // Begin chernoff bound computation
    mpfr_t upper_bound_for_one_vector;
    mpfr_init2(upper_bound_for_one_vector, N_PRECISION_BITS);
    mpfr_set_d(upper_bound_for_one_vector, 0, MPFR_RNDD);

    mpfr_t tmp_chernoff_bound;
    mpfr_init2(tmp_chernoff_bound, N_PRECISION_BITS);


    // printf("Computing bound\n");
    // int max_v = dist_dv_sum_eta2_dot_prod.max_value;
    // for (int v = dist_dv_sum_eta2_dot_prod.max_value; v >= dist_dv_sum_eta2_dot_prod.min_value; v--) {

    //     int idx = prob_dist_1d_get_idx_for_value(&dist_dv_sum_eta2_dot_prod, v);
    //     if (mpfr_cmp_d(dist_dv_sum_eta2_dot_prod.d[idx], ZERO_DOUBLE) <= 0)
    //         continue;

    //     max_v = v;
    //     break;

    // }
    // printf("max_v = %d\n", max_v);
    best_t = prob_dist_1d_minimize_t_for_chernoff(half_codeword_norm_squared, self_convolution_order, &dist_sum_dot_products_1d);
    prob_dist_1d_get_expectation_of_exp_tX(expectation_of_exp_tX, &dist_sum_dot_products_1d, best_t);
    for (int v = dist_dv_sum_eta2_dot_prod.min_value; v <= dist_dv_sum_eta2_dot_prod.max_value; v++) {

        int idx = prob_dist_1d_get_idx_for_value(&dist_dv_sum_eta2_dot_prod, v);
        if (mpfr_cmp_d(dist_dv_sum_eta2_dot_prod.d[idx], ZERO_DOUBLE) <= 0)
            continue;


        // if (dist_dv_sum_eta2_dot_prod.max_value - dist_dv_sum_eta2_dot_prod.min_value < 10000) {
        //     best_t = prob_dist_1d_minimize_t_for_chernoff(half_codeword_norm_squared - v, self_convolution_order, &dist_sum_dot_products_1d);
        //     prob_dist_1d_get_expectation_of_exp_tX(expectation_of_exp_tX, &dist_sum_dot_products_1d, best_t);
        // }

        prob_dist_1d_chernoff_with_precomputed_expectation_of_exp_tX(tmp_chernoff_bound, half_codeword_norm_squared - v,
                                                                     self_convolution_order, best_t,
                                                                     expectation_of_exp_tX);

        mpfr_mul(tmp_chernoff_bound, tmp_chernoff_bound, dist_dv_sum_eta2_dot_prod.d[idx], MPFR_RNDD);
        mpfr_add(upper_bound_for_one_vector, upper_bound_for_one_vector, tmp_chernoff_bound, MPFR_RNDD);
    }

    fprintf(stderr, "Chernoff bound is: ");
    mpfr_out_str(stdout, 10, 0, upper_bound_for_one_vector, ROUND_DOWN);
    fprintf(stderr, "\n");
    fflush(stderr);

    mpfr_clear(tmp_chernoff_bound2);
    mpfr_clear(expectation_of_exp_tX);
    mpfr_clear(upper_bound_for_one_vector);
    mpfr_clear(tmp_chernoff_bound);

    prob_dist_1d_clear(&dist_du0);
    prob_dist_1d_clear(&dist_dv);
    prob_dist_1d_clear(&dist_eta1);
    prob_dist_1d_clear(&dist_eta2);
    prob_dist_1d_clear(&dist_du0_sum_eta2);
    prob_dist_1d_clear(&dist_dv_sum_eta2);
    prob_dist_1d_clear(&dist_sum_dot_products_1d);
    prob_dist_1d_clear(&dist_dv_sum_eta2_dot_prod);

    mpc_matrix_free(&dist_sum_dot_products);
    mpc_matrix_free(&dist_e_r_product);
    mpc_matrix_free(&dist_s_du_sum_eta2_product);

    return 0;
}

int find_value_in_list(int value, int list[], int list_size) {
    for (int i = 0; i < list_size; i++) {
        if (list[i] == value)
            return i;
    }
    return -1;
}

void mpc_matrix_init_as_dv_eta2_dot_codeword(mpc_matrix_t *out, kyber_parameters_t *params, int codeword_value) {

    prob_dist_1d_t dist_dv;
    prob_dist_1d_init_as_compress_decompress_error(&dist_dv, params->KYBER_DV, params->KYBER_Q);
    prob_dist_1d_t dist_eta2;
    prob_dist_1d_init_as_centered_binomial(&dist_eta2, params->KYBER_ETA2);

    prob_dist_1d_t dist_dv_sum_eta2;
    prob_dist_1d_init_as_sum(&dist_dv_sum_eta2, &dist_dv, &dist_eta2);

    mpc_matrix_init(out, 1, MAX_N_FOR_FFT);

    for (int v = dist_dv_sum_eta2.min_value; v <= dist_dv_sum_eta2.max_value; v++) {
        int idx = prob_dist_1d_get_idx_for_value(&dist_dv_sum_eta2, v);

        int value = v * codeword_value;
        mpc_add_fr(out->m[0][MOD(value, out->n_cols)], out->m[0][MOD(value, out->n_cols)],
                   dist_dv_sum_eta2.d[idx], MPFR_RNDD);
    }
    fft(out->m[0], out->n_cols);

    prob_dist_1d_clear(&dist_dv);
    prob_dist_1d_clear(&dist_eta2);
    prob_dist_1d_clear(&dist_dv_sum_eta2);

}


void clear_precomputed_fft_dv_eta2_dot_codeword(kyber_parameters_t *params) {
    precomputed_fft_dv_eta2_dot_codeword_t *fft_table = &(params->precomputed_fft_dv_eta2_dot_codeword);

    for (int i = 0; i < fft_table->n_codeword_values; i++) {
        mpc_matrix_free(&fft_table->table[i]);
    }
    free(fft_table->table);

    fft_table->n_codeword_values = 0;
}


int init_kyber_parameters(kyber_parameters_t *params, FILE *experiment_setup) {
    char csv_line[MAX_STR] = {0};
    fgets(csv_line, MAX_STR, experiment_setup);

    if (strlen(csv_line) == 0 || strlen(csv_line) == 1) {
        // Empty line
        return 0;
    }
    // We discard (%*f) the input corresponding to the CODE_PNORM column
    int n_read = sscanf(csv_line, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%*f\n",
                        &params->KYBER_Q, &params->KYBER_LEVEL, &params->KYBER_K, &params->KYBER_ETA1,
                        &params->KYBER_ETA2, &params->KYBER_N_BLOCKS_FOR_DU0, &params->KYBER_DU0, &params->KYBER_DU1,
                        &params->KYBER_DV, &params->CODE_BETA);
    assert (n_read == 10);

    // KYBER_Q,KYBER_LEVEL,KYBER_K,KYBER_ETA1,KYBER_ETA2,KYBER_N_BLOCKS_FOR_DU0,KYBER_DU0,KYBER_DU1,KYBER_DV
    fprintf(stderr, "    KYBER_Q = %d\n", params->KYBER_Q);
    fprintf(stderr, "    KYBER_LEVEL = %d\n", params->KYBER_LEVEL);
    fprintf(stderr, "    KYBER_K = %d\n", params->KYBER_K);
    fprintf(stderr, "    KYBER_ETA1 = %d\n", params->KYBER_ETA1);
    fprintf(stderr, "    KYBER_ETA2 = %d\n", params->KYBER_ETA2);
    fprintf(stderr, "    KYBER_N_BLOCKS_FOR_DU0 = %d\n", params->KYBER_N_BLOCKS_FOR_DU0);
    fprintf(stderr, "    KYBER_DU0 = %d\n", params->KYBER_DU0);
    fprintf(stderr, "    KYBER_DU1 = %d\n", params->KYBER_DU1);
    fprintf(stderr, "    KYBER_DV = %d\n", params->KYBER_DV);
    fprintf(stderr, "    CODE_BETA = %d\n", params->CODE_BETA);
    params->CODE_DIMENSION = 4;

    assert(params->KYBER_N_BLOCKS_FOR_DU0 <= params->KYBER_K);
    assert(params->KYBER_N_BLOCKS_FOR_DU0 >= 0);

    return 1;
}


void clear_kyber_parameters(kyber_parameters_t *params) {
    params->KYBER_Q = 0;
    params->KYBER_LEVEL = 0;
    params->KYBER_K = 0;
    params->KYBER_ETA1 = 0;
    params->KYBER_ETA2 = 0;
    params->KYBER_N_BLOCKS_FOR_DU0 = 0;
    params->KYBER_DU0 = 0;
    params->KYBER_DU1 = 0;
    params->KYBER_DV = 0;
    params->CODE_DIMENSION = 0;
    params->CODE_BETA = 0;
}

void compute_4d_dfr_upper_bound(kyber_parameters_t *params) {
    mpfr_t total_dfr;
    mpfr_init2(total_dfr, N_PRECISION_BITS);
    mpfr_set_d(total_dfr, 0, MPFR_RNDD);

    mpfr_t upper_bound_for_one_vector;
    mpfr_init2(upper_bound_for_one_vector, N_PRECISION_BITS);
    mpfr_set_d(upper_bound_for_one_vector, 0, MPFR_RNDD);

    for (int i = 0; i < N_NEIGHBORS; i++) {

        int beta_index = find_value_in_list(params->CODE_BETA, INDEX_TO_BETA, N_BETAS);

        // Remember that we use padding to always have len(NEIGHBORS_WEIGHTS[beta_index]) == 46.
        // So we mark the weight as 0 when NEIGHBORS_WEIGHTS[beta_index][i] is used only for paadding, and
        // is not a real DFR relevant vector.
        if (NEIGHBORS_WEIGHTS[beta_index][i] == 0) continue;

        fprintf(stderr, "beta = %d\n", params->CODE_BETA);
        fflush(stderr);
        fprintf(stderr, "Beta = %d (index = %d)\n", params->CODE_BETA, beta_index);
        assert(beta_index >= 0);
        compute_error_probability_for_a_voronoi_relevant_vector(upper_bound_for_one_vector,
                                                                params, NEIGHBORS[beta_index][i]);
        mpfr_mul_si(upper_bound_for_one_vector, upper_bound_for_one_vector, NEIGHBORS_WEIGHTS[beta_index][i], MPFR_RNDD);
        mpfr_add(total_dfr, total_dfr, upper_bound_for_one_vector, MPFR_RNDD);

        fprintf(stderr, "Current DFR is: ");
        mpfr_out_str(stderr, 10, 0, total_dfr, ROUND_DOWN);
        fprintf(stderr, "\n");

    }
    // Each codeword is appears with probability 1/16
    mpfr_div_si(total_dfr, total_dfr, 16, MPFR_RNDD);

    // Union bound over all possible (256 / 4) 4D codewords used in Kyber
    mpfr_mul_si(total_dfr, total_dfr, 256 / 4, MPFR_RNDD);

    fprintf(stderr, "\n");
    fprintf(stderr, "Final DFR is: ");
    mpfr_out_str(stderr, 10, 0, total_dfr, ROUND_DOWN);
    fprintf(stderr, "\n\n");
    fprintf(stderr, "=========================\n");

    printf("%d,", params->KYBER_Q);
    printf("%d,", params->KYBER_LEVEL);
    printf("%d,", params->KYBER_K);
    printf("%d,", params->KYBER_ETA1);
    printf("%d,", params->KYBER_ETA2);
    printf("%d,", params->KYBER_N_BLOCKS_FOR_DU0);
    printf("%d,", params->KYBER_DU0);
    printf("%d,", params->KYBER_DU1);
    printf("%d,", params->KYBER_DV);
    printf("%d,", params->CODE_DIMENSION);
    printf("%d,", params->KYBER_Q/2);
    printf("%d,", params->CODE_BETA);
    mpfr_out_str(stdout, 10, 0, total_dfr, MPFR_RNDD);
    printf("\n");
    fflush(stdout);

    mpfr_clear(upper_bound_for_one_vector);
    mpfr_clear(total_dfr);
}

void run_experiment_compute_joint_distributions_and_save(const char experiment_setup_filename[]) {
    FILE *experiment_setup = fopen(experiment_setup_filename, "r");
    if (!experiment_setup) {
        fprintf(stderr, "Could not open experiment_setup file %s\n", experiment_setup_filename);
        exit(1);
    }

    char first_line[MAX_STR];
    fgets(first_line, MAX_STR, experiment_setup);
    assert(strncmp(first_line, EXPECTED_SETUP_CSV_HEADER, MAX_STR) == 0);

    print_dfr_csv_for_2d_codes_header();

    while (!feof(experiment_setup)) {
        fprintf(stderr, "==========================================================\n");
        fprintf(stderr, "Computing Pr(delta m[i, i + n/2]) for parameters:\n");
        kyber_parameters_t params = {0};
        int ret = init_kyber_parameters(&params, experiment_setup);
        if (ret == 0) {
            fprintf(stderr, "Found empty line in setup file, skipping it...\n");
            continue;
        }
        compute_4d_dfr_upper_bound(&params);
        clear_kyber_parameters(&params);
    }
    fclose(experiment_setup);
}


void print_dfr_csv_for_2d_codes_header() {
    printf("KYBER_Q,");
    printf("KYBER_LEVEL,");
    printf("KYBER_K,");
    printf("KYBER_ETA1,");
    printf("KYBER_ETA2,");
    printf("KYBER_N_BLOCKS_FOR_DU0,");
    printf("KYBER_DU0,");
    printf("KYBER_DU1,");
    printf("KYBER_DV,");
    printf("CODE_DIMENSION,");
    printf("CODE_ALPHA,");
    printf("CODE_BETA,");
    printf("DFR");
    printf("\n");
}

int main(int argc, char const *argv[]) {

    if (argc != 2) {
        fprintf(stderr, "Usage: %s experiment_setup.csv\n", argv[0]);
        exit(1);
    }
    if (SAVE_PARTIAL_DISTRIBUTION_RESULTS == 1)
        fprintf(stderr, "[*] Printing partial distributions to files.\n");

    init_fft_mpc();
    run_experiment_compute_joint_distributions_and_save(argv[1]);
    clear_fft_mpc();
}


