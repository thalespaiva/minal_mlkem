#!/usr/bin/python3
# Copyright 2025 LG Electronics, Inc. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0


import numpy as np

KYBER_Q = 3329

def get_equidistant_line(p1, p2):
    a, b = p1
    c, d = p2
    out = (2*(a-c), 2*(b-d), c**2 + d**2 - a**2 - b**2)
    if out[1] < 0:
        out = [-o for o in out]

    div = np.gcd.reduce(out)

    return [o//div for o in out]


def decode_x_y(xy, lines):
    true_x, true_y = xy

    reflect = False
    if (true_x > true_y):
        x, y = true_x, true_y

    else:
        x, y = true_y, true_x
        reflect = True

    # print(f'{true_x, true_y=}')
    # print(f'{x, y=}')
    above_l1 = (lines[1][1]*y > - lines[1][0]*x - lines[1][2])
    above_l2 = (lines[2][1]*y > - lines[2][0]*x - lines[2][2])
    above_l3 = (lines[3][1]*y > - lines[3][0]*x - lines[3][2])
    above_l4 = (lines[4][1]*y > - lines[4][0]*x - lines[4][2])
    above_l5 = (lines[5][1]*y > - lines[5][0]*x - lines[5][2])

    # print(f'{above_l1=}, {above_l2=}, {above_l3=}, {above_l4=}, {above_l5=}')
    c00 = (not above_l3 and above_l2 and above_l4);
    c01 = (not above_l2 and not above_l5 and not above_l3);
    c10 = (above_l2 and above_l3 and not above_l1);
    c11 = (above_l1) or (above_l3 and not above_l2) or (above_l5 and not above_l4);

    assert [c00, c01, c10, c11].count(True) == 1

    if c00:
        return (0, 0)
    if c11:
        return (1, 1)

    if not reflect:
        if (c01):
            return (0, 1)
        if (c10):
            return (1, 0)
    else:
        if (c01):
            return (1, 0)
        if (c10):
            return (0, 1)

    return None


def decode_x_y_positive_representative(xy, lines):
    true_x, true_y = xy

    reflect = False
    if (true_x > true_y):
        x, y = true_x, true_y

    else:
        x, y = true_y, true_x
        reflect = True

    # print(f'{true_x, true_y=}')
    # print(f'{x, y=}')
    above_l1 = (lines[1][1]*y > - lines[1][0]*x - lines[1][2])
    above_l2 = (lines[2][1]*y > - lines[2][0]*x - lines[2][2])
    above_l3 = (lines[3][1]*y > - lines[3][0]*x - lines[3][2])
    above_l4 = (lines[4][1]*y > - lines[4][0]*x - lines[4][2])
    above_l5 = (lines[5][1]*y > - lines[5][0]*x - lines[5][2])

    c00 = (not above_l5) or (not above_l2 and not above_l3) or (above_l4 and above_l1);
    c01 = (not above_l1 and not above_l3 and above_l2);
    c10 = (above_l3 and not above_l2 and above_l5);
    c11 = (above_l2 and above_l3 and not above_l4);

    assert [c00, c01, c10, c11].count(True) == 1

    if c00:
        return (0, 0)
    if c11:
        return (1, 1)

    if not reflect:
        if (c01):
            return (0, 1)
        if (c10):
            return (1, 0)
    else:
        if (c01):
            return (1, 0)
        if (c10):
            return (0, 1)

    return None

import sys

def ensure_balls_and_lines_are_correctly_related(lines, centered_balls):
    ball_representatives = [
        [(1, 1), (0, 1), (1, 1)],
        [(1, 0), (0, 0), (1, 0)],
        [(1, 1), (0, 1), (1, 1)],
    ]


    for i in range(3):
        for j in range(3):
            x, y = centered_balls[i][j]
            expected = ball_representatives[i][j]
            assert(expected == decode_x_y((x, y), lines))


    limits = [
        [(-(KYBER_Q//2 + 1), KYBER_Q//2), (0,  KYBER_Q//2), ( KYBER_Q//2,  KYBER_Q//2)],
        [(-(KYBER_Q//2 + 1), 0), (0, 0), ( KYBER_Q//2, 0)],
        [(-(KYBER_Q//2 + 1), -(KYBER_Q//2 + 1)), (0, -(KYBER_Q//2 + 1)), ( KYBER_Q//2, -(KYBER_Q//2 + 1))],
    ]
    for i in range(3):
        for j in range(3):
            x, y = limits[i][j]
            expected = ball_representatives[i][j]
            assert(expected == decode_x_y((x, y), lines))


def ensure_balls_and_lines_are_correctly_related_positive_representatives(lines, centered_balls):
    ball_representatives = [
        [(0, 0), (1, 0), (0, 0)],
        [(0, 1), (1, 1), (0, 1)],
        [(0, 0), (1, 0), (0, 0)],
    ]


    for i in range(3):
        for j in range(3):
            x, y = centered_balls[i][j]
            expected = ball_representatives[i][j]
            assert(expected == decode_x_y_positive_representative((x, y), lines))

    limits = [
        [(0, KYBER_Q), (KYBER_Q//2, KYBER_Q), (KYBER_Q, KYBER_Q)],
        [(0, KYBER_Q//2), (KYBER_Q//2, KYBER_Q//2), (KYBER_Q, KYBER_Q//2)],
        [(0, 0), (KYBER_Q//2, 0), (KYBER_Q, 0)],
    ]
    for i in range(3):
        for j in range(3):
            x, y = limits[i][j]
            expected = ball_representatives[i][j]
            assert(expected == decode_x_y_positive_representative((x, y), lines))


# This function generates the decoding lines when the points (x, y) to be decoded are in [-q/2, q/2] x [-q/2, q/2]
def get_decoding_lines(alpha, beta):

    base = np.array([
        np.array([alpha, beta]),
        np.array([beta, alpha]),
    ])

    zero = np.array((0, 0))

    centered_balls = [
        [base[0] + base[1] - np.array((KYBER_Q, 0)), base[1], base[0] + base[1]],
        [base[0] - np.array((KYBER_Q, 0)), [0, 0], base[0]],
    ]
    centered_balls.append([b - np.array((0, KYBER_Q)) for b in centered_balls[0]])
    # print(centered_balls)
    # LINES FORMAT:
    # l[0]*x + l[1] * y + l[2] = 0
    lines = []
    lines.append([-1, 1, 0])
    lines.append(get_equidistant_line(centered_balls[1][2], centered_balls[0][2]))
    lines.append(get_equidistant_line(centered_balls[1][1], centered_balls[2][1]))
    lines.append(get_equidistant_line(centered_balls[1][1], centered_balls[1][2]))
    lines.append(get_equidistant_line(centered_balls[2][0], centered_balls[1][1]))
    lines.append(get_equidistant_line(centered_balls[2][0], centered_balls[2][1]))

    # for l in lines:
    #     print(f'{l[1]}*y = {-l[0]}*x + ({-l[2]})')

    ensure_balls_and_lines_are_correctly_related(lines, centered_balls)

    return lines


# This function generates the decoding lines when the points (x, y) to be decoded are in [0, q] x [0, q]
def get_decoding_lines_positive_representatives(alpha, beta):

    base = np.array([
        np.array([alpha, beta]),
        np.array([beta, alpha]),
    ])

    zero = np.array((0, 0))

    centered_balls = [
        [base[0] + base[1] - np.array((KYBER_Q, 0)), base[1], base[0] + base[1]],
        [base[0] - np.array((KYBER_Q, 0)), [0, 0], base[0]],
    ]
    centered_balls.append([b - np.array((0, KYBER_Q)) for b in centered_balls[0]])


    centered_balls = [
        [zero + (0, KYBER_Q), base[0] + (0, KYBER_Q), zero + (KYBER_Q, KYBER_Q)],
        [base[1], base[0] + base[1], base[1] + (KYBER_Q, 0)],
        [zero, base[0], zero + (KYBER_Q, 0)],
    ]

    # print(centered_balls)
    # LINES FORMAT:
    # l[0]*x + l[1] * y + l[2] = 0
    lines = []
    lines.append([-1, 1, 0])
    lines.append(get_equidistant_line(centered_balls[1][2], centered_balls[0][2]))
    lines.append(get_equidistant_line(centered_balls[1][1], centered_balls[2][1]))
    lines.append(get_equidistant_line(centered_balls[1][1], centered_balls[1][2]))
    # lines.append(get_equidistant_line(centered_balls[2][0], centered_balls[1][1]))
    lines.append(get_equidistant_line(centered_balls[0][2], centered_balls[1][1]))
    lines.append(get_equidistant_line(centered_balls[2][0], centered_balls[2][1]))

    # for l in lines:
    #     print(f'{l[1]}*y = {-l[0]}*x + ({-l[2]})')

    ensure_balls_and_lines_are_correctly_related_positive_representatives(lines, centered_balls)

    return lines




def main_decoding_lines_for_performance_test_ref_implementation():
    print("//////////////////////////////////////////////////////////////////////////////////////////////////////////")
    print("//")
    print("// DECODING LINES FOR TESTING THE REFERENCE IMPLEMENTATION")
    print("//")
    print("// THIS SHOULD BE INCLUDED IN: kyber2d_kem_implementation/ref/poly.c ")
    print()
    first = True
    for k, beta in zip(range(2, 5), [422, 422, 436]):
        lines = get_decoding_lines(alpha=KYBER_Q//2, beta=beta)

        if first:
            print(f'#if KYBER_K == {k}')
            first = False
        else:
            print(f'#elif KYBER_K == {k}')
        print()
        print('#define CODE_ALPHA (KYBER_Q/2)')
        print(f'#define CODE_BETA {beta}')
        print()
        print('static int32_t L[5][3] = {')
        for i in range(1, 6):
            print(f'    {{{-lines[i][0]}, {lines[i][1]}, {-lines[i][2]}}},')
        print('};')
        print()

        for i in range(1, 6):
            print(f'#define ABOVE_L{i}(x, y) lower_than_mask(L[{i - 1}][0]*x + L[{i - 1}][2], L[{i - 1}][1]*y)')
        print()

    print('#endif')

def main_decoding_lines_for_performance_test_avx2_implementation():
    print("//////////////////////////////////////////////////////////////////////////////////////////////////////////")
    print("// ")
    print("// DECODING LINES FOR TESTING THE AVX2 IMPLEMENTATION")
    print("// ")
    print("// THIS SHOULD BE INCLUDED IN: kyber2d_kem_implementation/avx2/poly.c ")
    print()
    first = True
    for k, beta in zip(range(2, 5), [422, 422, 436]):
        lines = get_decoding_lines_positive_representatives(alpha=KYBER_Q//2, beta=beta)

        if first:
            print(f'#if KYBER_K == {k}')
            first = False
        else:
            print(f'#elif KYBER_K == {k}')
        print()
        print('#define CODE_ALPHA (KYBER_Q/2)')
        print(f'#define CODE_BETA {beta}')
        print()
        print('static int32_t L[5][3] = {')
        for i in range(1, 6):
            print(f'    {{{-lines[i][0]}, {lines[i][1]}, {-lines[i][2]}}},')
        print('};')
        print()

        for i in range(1, 6):
            print(f'#define ABOVE_L{i}(x, y) lower_than_mask(L[{i - 1}][0]*x + L[{i - 1}][2], L[{i - 1}][1]*y)')
        print()

    print('#endif')



def main_decoding_lines_for_dfr_computation():

    print("//////////////////////////////////////////////////////////////////////////////////////////////////////////")
    print("// THIS DEFINES THE FILE: kyber2d_dfr_analysis/decoding_lines.h ")
    print()
    print('#pragma once')
    print()

    print('static int32_t L[600][5][3] = {')

    for beta in range(600):
        if beta == 0:
            print('  // THE DECODING LINES ARE UNDEFINED FOR BETA = 0 (L[0])')
            print('  {')
            print('    {0,0,0},\n' * 5, end='')
            print('  },')
            continue
        print('  {')
        lines = get_decoding_lines(alpha=KYBER_Q//2, beta=beta)

        for i in range(1, 6):
            print(f'    {{{-lines[i][0]}, {lines[i][1]}, {-lines[i][2]}}},')
        print('  },')
    print('};')
    print()

    for i in range(1, 6):
        print(f'#define ABOVE_L{i}(beta, x, y) lower_than_mask(L[beta][{i - 1}][0]*x + L[beta][{i - 1}][2], L[beta][{i - 1}][1]*y)')
    print()


if __name__ == '__main__':
    print('\n' * 3)

    print('// BEGIN DECODING_LINES_FOR_DFR_COMPUTATION')
    main_decoding_lines_for_dfr_computation()
    print('// END DECODING_LINES_FOR_DFR_COMPUTATION')

    print('\n' * 3)

    print('// BEGIN DECODING_LINES_FOR_PERFORMANCE_TEST_REF_IMPLEMENTATION')
    main_decoding_lines_for_performance_test_ref_implementation()
    print('// END DECODING_LINES_FOR_PERFORMANCE_TEST_REF_IMPLEMENTATION')

    print('\n' * 3)

    print('// BEGIN DECODING_LINES_FOR_PERFORMANCE_TEST_AVX2_IMPLEMENTATION')
    main_decoding_lines_for_performance_test_avx2_implementation()
    print('// END DECODING_LINES_FOR_PERFORMANCE_TEST_AVX2_IMPLEMENTATION')
