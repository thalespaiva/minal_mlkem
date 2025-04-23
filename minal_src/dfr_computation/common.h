// Copyright 2025 LG Electronics, Inc. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <mpc.h>

#define MOD(a, b) ((((a) % (b)) + (b)) % (b))
#define N_PRECISION_BITS 260
// #define N_PRECISION_BITS 512
#define ROUND_DOWN MPC_RNDDN

#define ZERO_DOUBLE (1e-78)
#define N_SUPPORT_DIST_2D 2*2048

static __inline__ int center_mod(int x, int m) {
    x = MOD(x, m);
    if (x <= m - x) {
        return x;
    }
    return - (m - x);
}

#define MAX_CODE_BETA 500
