# Copyright 2025 LG Electronics, Inc. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

import numpy as np
import itertools as it
import pandas as pd
import sys

import matplotlib.pyplot as plt
import seaborn as sns
import scipy
import tqdm
import utils
from proba_util import KyberParameterSet

def gen_cyclic_matrix(row):
    mat = []
    n = len(row)
    for i in range(n):
        mat.append(row[n - i:] + row[:n - i])
    return mat


def center_mod(a, mod):
    a = int(a % mod)

    a2 = a - mod
    if (abs(a) < abs(a2)):
        return a
    return a2

def circular_shift(point, i):
    if i == 0:
        return point
    return point[-i:] + point[:len(point) - i]

# From: https://en.wikipedia.org/wiki/Ternary_search
def ternary_search(f, left, right, absolute_precision) -> float:
    """Left and right are the current bounds;
    the maximum is between them.
    """
    if abs(right - left) < absolute_precision:
        return (left + right) / 2

    left_third = (2*left + right) / 3
    right_third = (left + 2*right) / 3

    # print(left, left_third, right_third, right)
    # print([f(left), f(left_third), f(right_third), f(right)])
    if f(left_third) > f(right_third):
        return ternary_search(f, left_third, right, absolute_precision)
    else:
        return ternary_search(f, left, right_third, absolute_precision)


import math
# From: https://en.wikipedia.org/wiki/Ternary_search
def ternary_search_int(f, left, right) -> float:
    """Left and right are the current bounds;
    the maximum is between them.
    """
    if abs(right - left) <= 1:
    # if math.ceil(right) == math.floor(left):
        if f(right) > f(left):
            return right
        return left

    left_third = round((2*left + right) / 3)
    right_third = round((left + 2*right) / 3)

    print('[Ternary search]', left, left_third, right_third, right, file=sys.stderr)
    # print([f(left), f(left_third), f(right_third), f(right)])
    if f(left_third) < f(right_third):
        return ternary_search_int(f, left_third, right)
    else:
        return ternary_search_int(f, left, right_third)


class MinalCode:

    def init_general(self, q, parameters, lattice_code=False, norm_order=2):
        self.dim = len(parameters)
        self.parameters = parameters
        self.q = q
        self.generator = np.array(gen_cyclic_matrix(parameters))
        self.lattice_code = lattice_code
        self.norm_order = norm_order

    def __init__(self, n, q, beta, norm_order=2):
        assert(n >= 2)
        assert(q > 0)
        assert(q/2 > beta >= 0)

        self.dim = n
        self.beta = beta
        self.parameters = [0] * (n - 2) + [q//2, beta]
        self.q = q
        self.generator = np.array(gen_cyclic_matrix(self.parameters))
        self.norm_order = norm_order
        self.lattice_code = False

    def encode(self, v):
        assert len(v) == self.dim
        assert all([x in [0, 1] for x in v])

        return np.matmul(self.generator, np.array(v))

    def encode_ext(self, v):
        assert len(v) == self.dim

        if not self.lattice_code:
            base = np.array([self.q * (x//2) for x in v])
            # return base + np.matmul(self.generator, np.array([x % 2 for x in v]))
            return base + self.encode(np.array([x % 2 for x in v]))
        else:
            return np.matmul(self.generator, np.array(v))

    def get_neighbors(self, v):
        c = self.encode(np.array(v))

        neighbors = []
        for neighbor_diff in it.product([0, 1, -1, 2, -2], repeat=self.dim):
            neighbor_diff_array = np.array(neighbor_diff)

            if neighbor_diff == (0,) * self.dim:
                continue

            neighbor = self.encode_ext(np.array(v) + neighbor_diff_array)

            neighbors.append((np.array(v) + neighbor_diff_array, neighbor))

        return neighbors

    def minimum_distance(self, ord=None):

        if ord == None:
            ord = self.norm_order

        dists = []
        for word in it.product([0, 1], repeat=self.dim):
            c = self.encode(np.array(word))

            for word2 in it.product([0, 1], repeat=self.dim):
                if word == word2:
                    continue

                c2 = self.encode(np.array(word2))

                diff_mod_q = [center_mod(x - y, self.q) for x, y in zip(c, c2)]
                dists.append(np.linalg.norm(diff_mod_q, ord=ord))

        return min(dists)

    def find_best_beta(q, n_dimensions, norm_order=2):

        def f(beta):
            C = MinalCode(q=q, n=n_dimensions, beta=beta, norm_order=norm_order)
            return C.minimum_distance()

        return ternary_search_int(f, 0, q//2)




    def decode_bruteforce(self, point):
        point = np.array(point)

        dists = {}
        for word in it.product([0, 1], repeat=self.dim):
            codeword = self.encode(word)
            relative = codeword - point
            dists[word] = sum(center_mod(r, self.q)**2 for r in relative)

        return list(min(dists, key=lambda x: dists[x]))

    def decode_by_babai_rounding(self, point):
        w = np.linalg.solve(self.generator, point)
        # print(w)
        return [round(x) % 2 for x in w]

    def test_decode_by_babai_rounding(self):

        errors = {}
        for error in tqdm.tqdm(it.product(range(self.q), repeat=self.dim)):
            if self.decode_by_babai_rounding(error) != self.decode([error]):
                dist = int(np.linalg.norm(error))
                errors[dist] = errors.get(dist, 0) + 1

        return errors

    def decode(self, point):
        point = np.array(point)

        dists = []
        for word in it.product([-2, -1, 0, 1, 2, 3], repeat=self.dim):

            c = self.encode_ext(np.array(word))
            dist = np.linalg.norm(point - c)
            dists.append([word, dist])

        message, _ = min(dists, key=lambda x: x[1])
        return [i % 2 for i in message]

    def __str__(self):
        out = f'MinalCode of dimension {self.dim} and generator matrix:\n'
        for row in self.generator:
            out += f'    {row}\n'

        return out

    def __repr__(self):
        out = f'MinalCode of dimension {self.dim} and generator matrix:\n'
        for row in self.generator:
            out += f'    {row}\n'

        return out

    def get_neighbor_voronoi_cells(self, msg_coordinate):
        msg_coordinate = np.array(msg_coordinate)

        points = []
        center_codeword_idx = None
        # for i, delta_word in enumerate(it.product([-3, -2, -1, 0, 1, 2, 3], repeat=self.dim)):
        for i, delta_word in enumerate(it.product([-2, -1, 0, 1, 2], repeat=self.dim)):
            if all([w == 0 for w in delta_word]):
                center_codeword_idx = i
            points.append(self.encode_ext(msg_coordinate + np.array(delta_word)))

        vor = scipy.spatial.Voronoi(points)

        neighbors_idx = []
        for r in vor.ridge_points:
            r = [int(i) for i in r]

            if center_codeword_idx not in r:
                continue
            r.remove(center_codeword_idx)
            neighbors_idx.append(r[0])


        return [vor.points[i] for i in neighbors_idx]

    def get_dfr_relevant_points_with_weights(self):

        # Auxiliary function to find if there is a circular shift of a vector in a dictionary
        def find_shift_in_dict(c, dfr_estimation_points):
            for i in range(len(c)):
                candidate = tuple(circular_shift(c, i))
                if candidate in dfr_estimation_points:
                    return candidate
            return None

        dfr_estimation_points = {}
        for i, word in enumerate(it.product([0, 1], repeat=self.dim)):
            neighbors = self.get_neighbor_voronoi_cells(word)
            representative_codeword = self.encode(word)

            for neighbor in neighbors:
                c = tuple(map(int, neighbor - representative_codeword))

                # We consider all circular shifts of a vector to be the same point, since the error distribution
                # is the same for all of them (due to the symmetries of the base distributions).
                # Therefore we only need to increment the number of occurrences.
                c_key = find_shift_in_dict(c, dfr_estimation_points)

                if c_key is None:
                    dfr_estimation_points[c] = 1
                else:
                    dfr_estimation_points[c_key] += 1

        # We sort the DFR relevant points with respect to their norm, so that we get the worse
        # results (those that impact the DFR the most) first.
        sorted_points_weights = sorted(dfr_estimation_points.items(), key=lambda x: np.linalg.norm(x[0]))
        sorted_points, weights = zip(*sorted_points_weights)

        return sorted_points, weights

    def print_points_for_dfr_estimation(filepath):


        f = open(filepath, 'w');

        print('#define N_NEIGHBORS 46', file=f)
        print('static int NEIGHBORS[60][N_NEIGHBORS][4] = {', file=f)

        for beta in range(700, 761):
            GeneratorMatrix_4D = [0, 0, 1664, beta]
            C = MinalCode(q=q, parameters=GeneratorMatrix_4D)
            ps, ws = C.get_dfr_relevant_points_with_weights()
            print(len(ps))
            assert(len(ps) == 46)
            print(' {', file=f)
            for c in ps:
                print('  {', end='', file=f)
                for v in c:
                    print(f'{int(v)}, ', file=f, end='')
                print('},', file=f)
            print(' }', file=f),

        f.close()

    # This function should be used to generate file `experiments/setup/exp_kyber_minal4d_dfr_multiple_parameters.csv`
    def find_best_beta_for_kyber4d_params_csv(kyber_params_csv='experiments/setup/exp_kyber4d_dfr_without_betas.csv'):
        # ETA: 2:36.96
        df = pd.read_csv(kyber_params_csv)
        ps = []
        betas = []
        for i, row in df.iterrows():
            assert(row.KYBER_N_BLOCKS_FOR_DU0 == row.KYBER_K)
            assert(row.KYBER_DU1 == 0)
            kps = KyberParameterSet(
                n=256,
                k=row.KYBER_K,
                eta1=row.KYBER_ETA1,
                eta2=row.KYBER_ETA2,
                q=row.KYBER_Q,
                du=row.KYBER_DU0,
                dv=row.KYBER_DV
            )
            print(row, file=sys.stderr)
            error_dist = kps.get_final_error_distribution_1d()
            level_curve = utils.get_discrete_2d_level_curve(error_dist)
            f = lambda p: np.var([np.linalg.norm(x, ord=p) for x in level_curve])
            p = utils.ternary_search(f, 2, 100, 1e-6, min_max='min')
            ps.append(p)
            beta = MinalCode.find_best_beta(row.KYBER_Q, 4, norm_order=p)
            betas.append(beta)

            print(f'Best p-norm = {p} - Best beta = {beta}', file=sys.stderr)

        df['CODE_BETA'] = betas
        df['CODE_PNORM'] = ps

        df.to_csv(sys.stdout, index=False)

        return df


    def generate_4d_dfr_relevant_vectors_for_betas(filepath, kyber_params_csv='experiments/setup/exp_kyber4d_dfr.csv'):

        f = open(filepath, 'w');
        df = pd.read_csv(kyber_params_csv)

        KYBER_Q = 3329
        assert(all([x == KYBER_Q for x in df.KYBER_Q]))

        betas = sorted(set(df.CODE_BETA))
        print(f'{betas=}')

        print(f'#define N_BETAS {len(betas)}', file=f)
        print('#define N_NEIGHBORS 46', file=f)

        print(f'static int INDEX_TO_BETA[N_BETAS] = {{{', '.join(map(str, betas))}}};', file=f)

        print(f'static int NEIGHBORS[N_BETAS][N_NEIGHBORS][4] = {{', file=f)

        weights = []
        for beta in betas:
            C = MinalCode(q=KYBER_Q, n=4, beta=beta)
            ps, ws = C.get_dfr_relevant_points_with_weights()
            print(len(ps))
            assert(len(ps) <= 46)

            if len(ps) < 46:
                padding = (46 - len(ps))
                ps += ((KYBER_Q, KYBER_Q, KYBER_Q, KYBER_Q),) * padding
                ws += (0,) * padding

            print(' {', file=f)
            for c in ps:
                print('  {', end='', file=f)
                for v in c:
                    print(f'{int(v)}, ', file=f, end='')
                print('},', file=f)
            print(' },', file=f),

            weights.append(ws)

        print('};\n', file=f)
        print(f'static int NEIGHBORS_WEIGHTS[N_BETAS][N_NEIGHBORS] = {{', file=f)
        for w in weights:
            print(f'  {{{', '.join(map(str, w))}}},', file=f)

        print('};', file=f)

        f.close()

        return df, weights


# This function is used for generating results/minimum_distance/minimum_distance_p2.csv
def run_minimum_distance_experiment_multiple_dimensions(output='results/minimum_distance/minimum_distance_p2.csv'):

    # The dict below was computed using:
    # best_betas = {n: MinalCode.find_best_beta(3329, n, norm_order=2) for n in range(2, 9)}
    best_betas = {
        2: 446,
        3: 636,
        4: 751,
        5: 832,
        6: 893,
        7: 941,
        8: 981,
    }

    file = open(output, 'w')
    print('n,beta,minimum_distance', file=file)
    for n, best_beta in best_betas.items():
        for beta in range(best_beta - 50, best_beta + 51, 10):
            M = MinalCode(n=n, q=3329, beta=beta)
            print(f'{n},{beta},{M.minimum_distance()}', file=file)
    file.close()


if __name__ == '__main__':

    MinalCode.find_best_beta_for_kyber4d_params_csv()
