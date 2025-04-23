// Copyright 2025 LG Electronics, Inc. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

// This is the general (slow) decoder that is used for computing the best possible decoding.
void get_dfr_for_2d_code_with_minimum_distance_decoding(mpfr_t error_prob,
                                                        mpfr_matrix_t *delta_m,
                                                        int code_alpha,
                                                        int code_beta,
                                                        int kyber_q);

// This estimates the DFR using the real (fast and isochronous) decoder.
void get_dfr_for_2d_code_with_real_decoder(mpfr_t error_prob,
                                           mpfr_matrix_t *delta_m,
                                           int code_alpha,
                                           int code_beta,
                                           int kyber_q);
