# Copyright 2025 LG Electronics, Inc. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# From: https://en.wikipedia.org/wiki/Ternary_search
def ternary_search(f, left, right, absolute_precision, min_max='max') -> float:
    """Left and right are the current bounds;
    the maximum is between them.
    """
    if abs(right - left) < absolute_precision:
        return (left + right) / 2

    left_third = (2*left + right) / 3
    right_third = (left + 2*right) / 3

    # print('here', left, left_third, right_third, right)
    # print([f(left), f(left_third), f(right_third), f(right)])
    if min_max == 'min':
        if f(left_third) > f(right_third):
            return ternary_search(f, left_third, right, absolute_precision, min_max=min_max)
        else:
            return ternary_search(f, left, right_third, absolute_precision, min_max=min_max)
    elif min_max == 'max':
        if f(left_third) < f(right_third):
            return ternary_search(f, left_third, right, absolute_precision, min_max=min_max)
        else:
            return ternary_search(f, left, right_third, absolute_precision, min_max=min_max)
    else:
        raise ValueError('min_max should be either min or max')


def get_discrete_2d_level_curve(dist, level_curve_prob=2**(-128)):
    min_value, max_value = min(dist), max(dist)

    level_curve = []
    for v1 in range(min_value, max_value + 1):

        best_v2 = max_value
        best_v2_diff = abs(dist[v1] * dist[max_value] - level_curve_prob)
        best_v2_prob = abs(dist[v1] * dist[max_value])

        diff = abs(best_v2 - level_curve_prob)
        for v2 in range(max_value - 1, -1, -1):
            prob = dist[v1] * dist[v2]

            tmp_prob = abs(dist[v1] * dist[v2])
            if level_curve_prob/1.1 < tmp_prob <= 1.1*level_curve_prob:
                level_curve.append((v1, v2))
                continue

            if tmp_prob > 2*level_curve_prob:
                break
    lc = level_curve.copy()
    for v in reversed(lc):
        level_curve.append((v[0], -v[1]))
    level_curve.append(level_curve[0])
    return level_curve
