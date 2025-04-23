# Copyright 2025 LG Electronics, Inc. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

import pandas as pd


def generate_performance_table(ref_perf='./results/performance/ref_performance.csv',
                               avx2_perf='./results/performance/avx2_performance.csv',
                               m4_perf='./results/performance/m4_performance.csv',
                               a53_perf='./results/performance/a53_performance.csv',
                               function='kyber_decaps'):

    df_ref = pd.read_csv(ref_perf)
    df_avx2 = pd.read_csv(avx2_perf)
    df_m4 = pd.read_csv(m4_perf)
    df_a53 = pd.read_csv(a53_perf)

    levels_interest = [1, 3, 5, 5]
    dimension_interest = [4, 4, 2, 4]
    null_eta2_interest = [True, True, False, False]


    for impl, df in zip(['ref', 'avx2', 'm4', 'a53'], [df_ref, df_avx2, df_m4, df_a53]):
        for l, d, null_eta2 in zip(levels_interest, dimension_interest, null_eta2_interest):

            baseline_decaps_df = df[    (df.level == l) &
                                        (df.code_dimension == 1) &
                                        (df.null_eta2 == False) &
                                        (df.function == function)
                                    ]
            assert(len(baseline_decaps_df) == 1)
            baseline_decaps = baseline_decaps_df.iloc[0].cycles_avg

            decaps_df = df[    (df.level == l) &
                               (df.code_dimension == d) &
                               (df.null_eta2 == null_eta2) &
                               (df.function == function)
                        ]
            assert(len(decaps_df) == 1)
            decaps = decaps_df.iloc[0].cycles_avg

            print(f'{impl, l, d, null_eta2=}, speedup: {(baseline_decaps / decaps):.2f}')


def generate_pk_compression_performance_table(
    avx2_perf='./results/performance/avx2_performance.csv',
    m4_perf='./results/performance/m4_performance.csv',
    a53_perf='./results/performance/a53_performance.csv',
    function='kyber_encaps'):

    df_avx2 = pd.read_csv(avx2_perf)
    df_m4 = pd.read_csv(m4_perf)
    df_a53 = pd.read_csv(a53_perf)

    levels_interest = [1, 3, 5]
    dimension_interest = [4, 4, 4]
    null_eta2_interest = [False, False, False]

    ks = {1: 2, 3: 3, 5: 4}

    # We use the polyvec_decompress from level 5 (which decompresses 11 bit values)
    # to estimate the extra decompression time for u for each parameter set.
    # Dictionary poly_decompress_11bit_times normalizes this value for one polynomial,
    # so we have to divide it by k (k = 4 for level 5), and then multiply by
    # the appropriate k later. (See below `if (function == 'kyber_encaps')`)
    poly_decompress_11bit_times = {
        # AVX2 raw row:
        # level,function,cycles_avg,cycles_median,code_dimension,null_eta2
        # 5,polyvec_decompress,170,170,1,True
        'avx2': 170 / 4,

        # M4 raw row:
        # function,cycles_avg,cycles_stddev
        # polyvec_decompress_1024,13525,16
        'm4': 13525 / 4,

        # A53 raw row:
        # level,function,cycles_avg,cycles_stddev
        # 5,polyvec_decompress,6132,393
        'a53': 6132 / 4,
    }

    for impl, df in zip(['avx2', 'm4', 'a53'], [df_avx2, df_m4, df_a53]):

        for l, d, null_eta2 in zip(levels_interest, dimension_interest, null_eta2_interest):


            baseline_func_df = df[    (df.level == l) &
                                        (df.code_dimension == 1) &
                                        (df.null_eta2 == False) &
                                        (df.function == function)
                                    ]
            assert(len(baseline_func_df) == 1)
            baseline_func = baseline_func_df.iloc[0].cycles_avg

            func_df = df[    (df.level == l) &
                               (df.code_dimension == d) &
                               (df.null_eta2 == null_eta2) &
                               (df.function == function)
                        ]
            assert(len(func_df) == 1)
            func = func_df.iloc[0].cycles_avg
            if (function == 'kyber_encaps'):
                func += poly_decompress_11bit_times[impl] * ks[l]

            print(f'{impl, l, d, null_eta2=}, speedup: {(baseline_func / func):.2f}')



def generate_code_performance_table(ref_perf='./results/performance/ref_performance.csv',
                                    avx2_perf='./results/performance/avx2_performance.csv',
                                    m4_perf='./results/performance/m4_performance.csv',
                                    a53_perf='./results/performance/a53_performance.csv',
                                    function='kyber_decaps'):

    df_ref = pd.read_csv(ref_perf)
    df_avx2 = pd.read_csv(avx2_perf)
    df_m4 = pd.read_csv(m4_perf)
    df_a53 = pd.read_csv(a53_perf)

    for function in ['poly_frommsg', 'poly_tomsg']:
        for impl, df in zip(['avx2', 'm4', 'a53'], [df_avx2, df_m4, df_a53]):

            for idx, g in df[df.function == function].groupby('code_dimension'):
                print(function, impl, idx, round(g.cycles_avg.mean()))
