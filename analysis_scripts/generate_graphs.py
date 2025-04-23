#!/usr/bin/python3
# Copyright 2025 LG Electronics, Inc. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0


import math
import os
import pickle
import utils

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
import seaborn as sns
import itertools as it
import scipy

from matplotlib import rc
from matplotlib import rcParams

from proba_util import KyberParameterSet

import sys

#TCHES text width: 5.37502in
TCHES_FIG_WIDTH = 5.39
# TCHES_BIG_FIG_HEIGHT = 5.37502
TCHES_SMALL_FIG_HEIGHT = 2.5

TCHES_HALF_PAGE_FIG = 2.2
TCHES_HALF_PAGE_FIG_BIG = 2.4
TCHES_HALF_PAGE_FIG_SMALL = 2

KYBER_Q = 3329


def latexify_tches(fig_width=TCHES_FIG_WIDTH, fig_height=TCHES_SMALL_FIG_HEIGHT, small=False, fs=9):

    sns.reset_orig()

    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)
    rcParams['font.size'] = fs
    from math import sqrt

    if fig_width is None:
        fig_width = 4.8  # approx 12.2 cm

    if fig_height is None:
        golden_mean = (sqrt(5) - 1.0) / 2.0  # Aesthetic ratio
        fig_height = fig_width * golden_mean  # height in inches

    MAX_HEIGHT_INCHES = 8.0
    if fig_height > MAX_HEIGHT_INCHES:
        print("WARNING: fig_height too large:" + fig_height +
              "so will reduce to" + MAX_HEIGHT_INCHES + "inches.")
        fig_height = MAX_HEIGHT_INCHES

    # ticksize = fs
    ticksize = 7
    if small:
        ticksize = 7

    params = {
        'backend': 'pdf',
        'axes.labelsize': fs,  # fontsize for x and y labels (was 10)
        'axes.titlesize': fs,
        'font.size': fs,  # was 10
        'legend.fontsize': fs,  # was 10
        'font.family': 'serif',
        'xtick.labelsize': ticksize,
        'ytick.labelsize': ticksize,
        'ytick.major.width': 0.3,
        'ytick.minor.width': 0.3,
        'xtick.major.width': 0.3,
        'xtick.minor.width': 0.3,
        'text.usetex': True,
        'legend.edgecolor': 'none',
        'legend.frameon': True,
        'legend.framealpha': 0.5,
        'legend.fancybox': False,
        'axes.linewidth': 0.3,
        'axes.linewidth': 0.3,
        'figure.figsize': [fig_width, fig_height],
    }

    rcParams.update(params)



class CodePlots:


    @staticmethod
    def get_distance(a, b, pnorm=2):
        # return math.sqrt((a[0] - b[0])**2 + (a[1] - b[1])**2)
        return np.linalg.norm(np.array(a) - np.array(b), ord=pnorm)

    @staticmethod
    def extend_points(encoded_points):
        points = []
        for p in encoded_points:
            points.append(p)
            points.append((p[0] + KYBER_Q, p[1]))
            points.append((p[0], p[1] + KYBER_Q))
            points.append((p[0] + KYBER_Q, p[1] + KYBER_Q))
            points.append((p[0] - KYBER_Q, p[1]))
            points.append((p[0], p[1] - KYBER_Q))
            points.append((p[0] - KYBER_Q, p[1] - KYBER_Q))
        return points

    @staticmethod
    def get_minimum_distance(points, pnorm=2):
        dists = []
        # points = encoded_points + [(KYBER_Q, KYBER_Q), (0, KYBER_Q), (KYBER_Q, 0)]
        # points = extend_points(encoded_points)

        for a, b in it.combinations(points, 2):
            d = CodePlots.get_distance(a, b, pnorm=pnorm)
            # print(a, b, d)
            dists.append(d)

        return min(dists)

    @staticmethod
    def draw_points():

        # X = 390
        # X = 446
        X = 422
        # X = 2
        Y = X
        encoded_points = [
            (0, 0),
            (X, KYBER_Q//2),
            (KYBER_Q//2, X),
            (KYBER_Q//2 + Y, KYBER_Q//2 + Y),
        ]
        points = CodePlots.extend_points(encoded_points)

        xs, ys = zip(*points)
        plt.scatter(xs, ys)
        ax = plt.gca()
        ax.set_aspect('equal', adjustable='box')
        min_dist = CodePlots.get_minimum_distance(points)
        print(min_dist)
        for p in points:
            circle = plt.Circle(p, min_dist/2, facecolor=None, fill=False, edgecolor='r')
            ax.add_patch(circle)

    def draw_code(X=422, Y=KYBER_Q//2, lattice_code=False, figsize=(TCHES_HALF_PAGE_FIG_SMALL, TCHES_HALF_PAGE_FIG_SMALL)):

        fig, ax = plt.subplots(1, figsize=figsize)


        base = np.array([
            np.array((Y, X)),
            np.array((X, Y)),
        ])

        base_kyber = np.array([
            np.array((0, KYBER_Q)),
            np.array((KYBER_Q, 0)),
        ])
        points = []
        points_idx = []

        if not lattice_code:
            for ik, jk in it.product(range(-3, 3), repeat=2):
                kyber_point = np.matmul(base_kyber, np.array([ik, jk]))
                for i, j in it.product(range(0, 2), repeat=2):
                    points_idx.append((i % 2, j % 2))
                    points.append(kyber_point + np.matmul(base, np.array([i, j])))
        else:
            for i, j in it.product(range(-6, 6), repeat=2):
                points_idx.append((i % 2, j % 2))
                points.append(np.matmul(base, np.array([i, j])))

        xs, ys = zip(*points)
        # plt.scatter(xs, ys, marker='.', s=2, c='k')
        ax = plt.gca()
        ax.set_aspect('equal')
        min_dist = CodePlots.get_minimum_distance(points)
        print(min_dist)

        corepoints = set()
        for ik, jk in it.product(range(2), repeat=2):
            corepoints.add(tuple(np.matmul(base, np.array([ik, jk]))))

        ghost_point_color = '0.5'
        for p in points:
            x, y = p
            if tuple(p) in corepoints:
                plt.scatter([x], [y], marker='.', s=4, c='k')
                circle = plt.Circle(p, min_dist/2, facecolor=None, fill=False,
                                    edgecolor='k', lw=0.7, ls='--')
            else:
                plt.scatter([x], [y], marker='.', s=4, c=ghost_point_color)
                circle = plt.Circle(p, min_dist/2, facecolor=None, fill=False,
                                    edgecolor='k', lw=0.7, ls='--')
                                    # edgecolor=ghost_point_color, lw=0.3, ls='--')
            ax.add_patch(circle)

        # ax.spines['left'].set_position(('data', 0))
        # ax.spines['bottom'].set_position(('data', 0))

        sns.despine()
        # plt.xticks([i*KYBER_Q for i in range(-5, 5)])
        # plt.yticks([i*KYBER_Q for i in range(-5, 5)])
        # plt.xticks([])
        # plt.yticks([])

        plt.xticks([-KYBER_Q/2, 0, KYBER_Q/2, KYBER_Q, 3*KYBER_Q/2], ['$-q/2$', '$0$', '$q/2$', '$q$', '$3q/2$'])
        plt.yticks([-KYBER_Q/2, 0, KYBER_Q/2, KYBER_Q, 3*KYBER_Q/2], ['$-q/2$', '$0$', '$q/2$', '$q$', '$3q/2$'])

        plt.xlim(-KYBER_Q/2, KYBER_Q + KYBER_Q/2 + 1)
        plt.ylim(-KYBER_Q/2, KYBER_Q + KYBER_Q/2 + 1)

        def in_limit(point):
            xlim, ylim = plt.xlim(), plt.ylim()

            if not (xlim[0] - 1 <= point[0] < xlim[1]):
                return False
            if not (ylim[0] - 1<= point[1] < ylim[1]):
                return False
            return True


        for (ik, jk), point in zip(points_idx, points):
            if not in_limit(point):
                continue

            if tuple(point) in corepoints:
                ax.text(x=point[0] + 30, y=point[1] + 30, s=f'{ik}{jk}',  fontsize=7)
            else:
                ax.text(x=point[0] + 30, y=point[1] + 30, s=f'{ik}{jk}',  fontsize=7, c=ghost_point_color)


        # xlabelpad = ax.xaxis.labelpad
        # ylabelpad = ax.yaxis.labelpad
        # for val, letter in [(-KYBER_Q, '-q'), (KYBER_Q, 'q')]:
        #     ax.annotate('${}$'.format(letter), xy=(val, 0), xytext=(0, -xlabelpad),
        #                 # xycoords=('data', 'axes fraction'),
        #                 textcoords='offset points',
        #                 ha='right' if val > 0 else 'left', va='top')
        #     ax.annotate('${}$'.format(letter), xy=(0, val), xytext=(-ylabelpad, 0),
        #                 # xycoords=('data', 'axes fraction'),
        #                 textcoords='offset points',
        #                 ha='right', va='top' if val > 0 else 'bottom')

        rect = plt.Rectangle((0, 0), KYBER_Q, KYBER_Q, facecolor=(0,0,0,0.1), edgecolor='k', lw=0.3)
        ax.add_patch(rect)

        plt.tight_layout()
        plt.subplots_adjust(top=0.90, bottom=0.1, left=0.2, right=0.9)



    def draw_code_with_voronoi(X=422, Y=KYBER_Q//2, figsize=(TCHES_HALF_PAGE_FIG_SMALL, TCHES_HALF_PAGE_FIG_SMALL)):

        fig, ax = plt.subplots(1, figsize=figsize)

        base = np.array([
            np.array((Y, X)),
            np.array((X, Y)),
        ])

        base_kyber = np.array([
            np.array((0, KYBER_Q)),
            np.array((KYBER_Q, 0)),
        ])
        points = []
        points_idx = []

        for ik, jk in it.product(range(-3, 3), repeat=2):
            kyber_point = np.matmul(base_kyber, np.array([ik, jk]))
            for i, j in it.product(range(0, 2), repeat=2):
                points_idx.append((i % 2, j % 2))
                points.append(kyber_point + np.matmul(base, np.array([i, j])))

        vor = scipy.spatial.Voronoi(points)
        scipy.spatial.voronoi_plot_2d(vor, ax=ax, show_vertices=False, show_points=False, line_colors='k', line_width=0.5)

        # return
        xs, ys = zip(*points)
        # plt.scatter(xs, ys, marker='.', s=2, c='k')
        ax = plt.gca()
        ax.set_aspect('equal')
        min_dist = CodePlots.get_minimum_distance(points)
        print(min_dist)

        corepoints = set()
        for ik, jk in it.product(range(2), repeat=2):
            corepoints.add(tuple(np.matmul(base, np.array([ik, jk]))))

        ghost_point_color = '0.5'
        for p in points:
            x, y = p
            if tuple(p) in corepoints:
                plt.scatter([x], [y], marker='.', s=4, c='k')
                circle = plt.Circle(p, min_dist/2, facecolor=None, fill=False,
                                    edgecolor='k', lw=0.7, ls='--')
            else:
                plt.scatter([x], [y], marker='.', s=4, c=ghost_point_color)
                circle = plt.Circle(p, min_dist/2, facecolor=None, fill=False,
                                    edgecolor='k', lw=0.7, ls='--')
                                    # edgecolor=ghost_point_color, lw=0.3, ls='--')
            ax.add_patch(circle)

        ax.spines['left'].set_position(('data', 0))
        ax.spines['bottom'].set_position(('data', 0))

        sns.despine()
        # plt.xticks([i*KYBER_Q for i in range(-5, 5)])
        # plt.yticks([i*KYBER_Q for i in range(-5, 5)])
        plt.xticks([])
        plt.yticks([])
        plt.xlim(-KYBER_Q/2, KYBER_Q + KYBER_Q/2 + 1)
        plt.ylim(-KYBER_Q/2, KYBER_Q + KYBER_Q/2 + 1)

        def in_limit(point):
            xlim, ylim = plt.xlim(), plt.ylim()

            if not (xlim[0] - 1 <= point[0] < xlim[1]):
                return False
            if not (ylim[0] - 1<= point[1] < ylim[1]):
                return False
            return True


        for (ik, jk), point in zip(points_idx, points):
            if not in_limit(point):
                continue

            if tuple(point) in corepoints:
                ax.text(x=point[0] + 30, y=point[1] + 30, s=f'{ik}{jk}',  fontsize=7)
            else:
                ax.text(x=point[0] + 30, y=point[1] + 30, s=f'{ik}{jk}',  fontsize=7, c=ghost_point_color)


        xlabelpad = ax.xaxis.labelpad
        ylabelpad = ax.yaxis.labelpad
        for val, letter in [(-KYBER_Q, '-q'), (KYBER_Q, 'q')]:
            ax.annotate('${}$'.format(letter), xy=(val, 0), xytext=(0, -xlabelpad),
                        # xycoords=('data', 'axes fraction'),
                        textcoords='offset points',
                        ha='right' if val > 0 else 'left', va='top')
            ax.annotate('${}$'.format(letter), xy=(0, val), xytext=(-ylabelpad, 0),
                        # xycoords=('data', 'axes fraction'),
                        textcoords='offset points',
                        ha='right', va='top' if val > 0 else 'bottom')

        rect = plt.Rectangle((0, 0), KYBER_Q, KYBER_Q, facecolor=(0,0,0,0.1), edgecolor='k', lw=0.3)
        ax.add_patch(rect)
        plt.tight_layout()


    def get_line_equations_for_good_code(X=422, Q=KYBER_Q):

        def eq_dist(p1, p2):
            a, b = p1
            c, d = p2
            out = (2*(a-c), 2*(b-d), c**2 + d**2 - a**2 - b**2)
            if out[1] < 0:
                out = [-o for o in out]

            return out

        # Q = KYBER_Q
        base = np.array([
            np.array([Q//2, X]),
            np.array([X, Q//2]),
        ])
        zero = np.array((0, 0))

        BS = [
            [base[0] + base[1] - np.array((Q, 0)), base[1], base[0] + base[1]],
            [base[0] - np.array((Q, 0)), [0, 0], base[0]],
        ]
        BS.append([b - np.array((0, Q)) for b in BS[0]])

        # return BS
        lines = []
        lines.append([-1, 1, 0])
        lines.append(eq_dist(BS[1][2], BS[0][2]))
        lines.append(eq_dist(BS[1][1], BS[2][1]))
        lines.append(eq_dist(BS[1][1], BS[1][2]))
        lines.append(eq_dist(BS[2][0], BS[1][1]))
        lines.append(eq_dist(BS[2][0], BS[2][1]))



        for l in lines:
            print(f'{l[1]}*y = {-l[0]}*x + ({-l[2]})')
        return lines


    def draw_good_2d_code_decoding_lines(X=422, figsize=(TCHES_HALF_PAGE_FIG, TCHES_HALF_PAGE_FIG)):

        fig, ax = plt.subplots(1, figsize=figsize)

        base = np.array([
            np.array((KYBER_Q//2, X)),
            np.array((X, KYBER_Q//2)),
        ])

        base_kyber = np.array([
            np.array((0, KYBER_Q)),
            np.array((KYBER_Q, 0)),
        ])
        points = []
        points_idx = []

        for ik, jk in it.product(range(-3, 3), repeat=2):
            kyber_point = np.matmul(base_kyber, np.array([ik, jk]))
            for i, j in it.product(range(0, 2), repeat=2):
                points_idx.append((i % 2, j % 2))
                points.append(kyber_point + np.matmul(base, np.array([i, j])))

        xs, ys = zip(*points)
        # plt.scatter(xs, ys, marker='.', s=2, c='k')
        ax = plt.gca()
        ax.set_aspect('equal')
        min_dist = CodePlots.get_minimum_distance(points)
        print(min_dist)

        corepoints = set()
        for ik, jk in it.product(range(2), repeat=2):
            corepoints.add(tuple(np.matmul(base, np.array([ik, jk]))))

        ghost_point_color = '0.5'
        for p in points:
            x, y = p
            # if tuple(p) in corepoints:
            # plt.scatter([x], [y], marker='.', s=4, c='k')
            # circle = plt.Circle(p, min_dist/2, facecolor=None, fill=False,
            #                     edgecolor='k', lw=0.7, ls='--')
            # else:
            plt.scatter([x], [y], marker='.', s=4, c=ghost_point_color)
            circle = plt.Circle(p, min_dist/2, facecolor=None, fill=False,
                                # edgecolor='k', lw=0.7, ls='--')
                                edgecolor=ghost_point_color, lw=0.3, ls='--')
            ax.add_patch(circle)

        ax.spines['left'].set_position(('data', 0))
        ax.spines['bottom'].set_position(('data', 0))

        sns.despine()
        # plt.xticks([i*KYBER_Q for i in range(-5, 5)])
        # plt.yticks([i*KYBER_Q for i in range(-5, 5)])
        plt.xticks([])
        plt.yticks([])
        plt.xlim(-KYBER_Q/2 - 500, KYBER_Q/2 + 500)
        plt.ylim(-KYBER_Q/2 -500, KYBER_Q/2 + 500)

        def in_limit(point):
            xlim, ylim = plt.xlim(), plt.ylim()

            if not (xlim[0] - 1 <= point[0] < xlim[1]):
                return False
            if not (ylim[0] - 1<= point[1] < ylim[1]):
                return False
            return True


        for (ik, jk), point in zip(points_idx, points):
            if not in_limit(point):
                continue

            # if tuple(point) in corepoints:
            #     ax.text(x=point[0] + 30, y=point[1] + 30, s=f'{ik}{jk}')
            # else:
            ax.text(x=point[0] + 30, y=point[1] + 30, s=f'{ik}{jk}', fontsize=7, c=ghost_point_color)


        xlabelpad = ax.xaxis.labelpad
        ylabelpad = ax.yaxis.labelpad
        for val, letter in [(-KYBER_Q/2, '-q/2'), (KYBER_Q/2, 'q/2')]:
            ax.annotate('${}$'.format(letter), xy=(val, 0), xytext=(0, +xlabelpad),
                        # xycoords=('data', 'axes fraction'),
                        textcoords='offset points',
                        ha='left' if val > 0 else 'right', va='center')
            # ax.annotate('${}$'.format(letter), xy=(0, val), xytext=(-ylabelpad, 0),
            #             # xycoords=('data', 'axes fraction'),
            #             textcoords='offset points',
            #             ha='right', va='top' if val > 0 else 'bottom')

        rect = plt.Rectangle((-KYBER_Q/2, -KYBER_Q/2), KYBER_Q, KYBER_Q, fill=None, edgecolor='k', lw=0.3)
        triangle = plt.Polygon([(-KYBER_Q/2, -KYBER_Q/2), (KYBER_Q/2, -KYBER_Q/2), (KYBER_Q/2, KYBER_Q/2)], facecolor=(0,0,0,0.1), edgecolor='k', lw=0.3)
        ax.add_patch(rect)
        ax.add_patch(triangle)

        decoding_lines = CodePlots.get_line_equations_for_good_code(X)


           # ax.annotate('${}$'.format(letter), xy=(val, 0), xytext=(0, -xlabelpad),
           #              # xycoords=('data', 'axes fraction'),
           #              textcoords='offset points',
           #              ha='right' if val > 0 else 'left', va='top')
           #  ax.annotate('${}$'.format(letter), xy=(0, val), xytext=(-ylabelpad, 0),
           #              # xycoords=('data', 'axes fraction'),
           #              textcoords='offset points',
           #              ha='right', va='top' if val > 0 else 'bottom')
        def put_label_if_in_limit(i, x_text, y_text, get_x_from_y, get_y_from_x):
            if i == 0:
                return

            if (y_text <= -KYBER_Q/2 and x_text >= -KYBER_Q/2):
                pady = 200
                padx = 100
                # return ax.text(x=x_text, y=y_text, s=fr'$\ell_{i}$')
                ax.annotate(fr'$\ell_{i}$', xy=(x_text, y_text),
                            xytext=(padx + get_x_from_y(y_text - pady), y_text - pady)),

            if (x_text >= KYBER_Q/2):
                padx = 200
                pady = 100

                ax.annotate(fr'$\ell_{i}$', xy=(x_text, y_text),
                            xytext=(x_text + padx, get_y_from_x(x_text + padx) + pady)),
                # return ax.text(x=x_text, y=y_text, s=fr'$\ell_{i}$')
            return


        for i, l in enumerate(decoding_lines):
            a_x, a_y, a_c = l
            x0 = -KYBER_Q
            y0 = (-a_x * x0 - a_c)/a_y
            x1 = KYBER_Q
            y1 = (-a_x * x1 - a_c)/a_y
            plt.axline((x0, y0), (x1, y1), lw=0.3, c='k')

            def get_x_from_y(y_text):
                return (-a_y * y_text - a_c) / a_x

            def get_y_from_x(x_text):
                return (-a_x * x_text - a_c) / a_y

            for y_text in [-KYBER_Q/2, KYBER_Q/2]:
                x_text = get_x_from_y(y_text)
                put_label_if_in_limit(i, x_text, y_text, get_x_from_y, get_y_from_x)

            for x_text in [-KYBER_Q/2, KYBER_Q/2]:
                y_text = get_y_from_x(x_text)
                put_label_if_in_limit(i, x_text, y_text, get_x_from_y, get_y_from_x)

        plt.tight_layout()

    def draw_good_2d_code_decoding_lines_voronoi(X=422, figsize=(TCHES_HALF_PAGE_FIG, TCHES_HALF_PAGE_FIG)):

        fig, ax = plt.subplots(1, figsize=figsize)

        base = np.array([
            np.array((KYBER_Q//2, X)),
            np.array((X, KYBER_Q//2)),
        ])

        base_kyber = np.array([
            np.array((0, KYBER_Q)),
            np.array((KYBER_Q, 0)),
        ])
        points = []
        points_idx = []

        for ik, jk in it.product(range(-3, 3), repeat=2):
            kyber_point = np.matmul(base_kyber, np.array([ik, jk]))
            for i, j in it.product(range(0, 2), repeat=2):
                points_idx.append((i % 2, j % 2))
                points.append(kyber_point + np.matmul(base, np.array([i, j])))

        vor = scipy.spatial.Voronoi(points)
        scipy.spatial.voronoi_plot_2d(vor, ax=ax, show_vertices=False, show_points=False, line_colors='0.4', line_width=0.3)

        xs, ys = zip(*points)
        # plt.scatter(xs, ys, marker='.', s=2, c='k')
        ax = plt.gca()
        ax.set_aspect('equal')
        min_dist = CodePlots.get_minimum_distance(points)
        print(min_dist)

        corepoints = set()
        for ik, jk in it.product(range(2), repeat=2):
            corepoints.add(tuple(np.matmul(base, np.array([ik, jk]))))

        ghost_point_color = '0.6'
        for p in points:
            x, y = p
            # if tuple(p) in corepoints:
            # plt.scatter([x], [y], marker='.', s=4, c='k')
            # circle = plt.Circle(p, min_dist/2, facecolor=None, fill=False,
            #                     edgecolor='k', lw=0.7, ls='--')
            # else:
            plt.scatter([x], [y], marker='.', s=4, c=ghost_point_color)
            circle = plt.Circle(p, min_dist/2, facecolor=None, fill=False,
                                # edgecolor='k', lw=0.7, ls='--')
                                edgecolor=ghost_point_color, lw=0.3, ls='--')
            ax.add_patch(circle)

        # ax.spines['left'].set_position(('data', 0))
        # ax.spines['bottom'].set_position(('data', 0))

        sns.despine()
        # plt.xticks([i*KYBER_Q for i in range(-5, 5)])
        # plt.yticks([i*KYBER_Q for i in range(-5, 5)])
        # plt.xticks([])
        # plt.yticks([])

        plt.xticks([-KYBER_Q/2, 0, KYBER_Q/2], ['$-q/2$', '$0$', '$q/2$'])
        plt.yticks([-KYBER_Q/2, 0, KYBER_Q/2], ['$-q/2$', '$0$', '$q/2$'])

        plt.xlim(-KYBER_Q/2 - 500, KYBER_Q/2 + 500)
        plt.ylim(-KYBER_Q/2 -500, KYBER_Q/2 + 500)

        def in_limit(point):
            xlim, ylim = plt.xlim(), plt.ylim()

            if not (xlim[0] - 1 <= point[0] < xlim[1]):
                return False
            if not (ylim[0] - 1<= point[1] < ylim[1]):
                return False
            return True


        for (ik, jk), point in zip(points_idx, points):
            if not in_limit(point):
                continue

            # if tuple(point) in corepoints:
            #     ax.text(x=point[0] + 30, y=point[1] + 30, s=f'{ik}{jk}')
            # else:
            ax.text(x=point[0] + 30, y=point[1] + 30, s=f'{ik}{jk}', fontsize=6, c=ghost_point_color)

        ax.text(x=KYBER_Q/2+X, y=KYBER_Q/2+X, s=f'$A$', fontsize=7, c='k', horizontalalignment='right', verticalalignment='bottom')
        ax.text(x=0, y=0, s=f'$B$', fontsize=7, c='k', horizontalalignment='right', verticalalignment='bottom')
        ax.text(x=KYBER_Q/2+X- KYBER_Q, y=KYBER_Q/2+X - KYBER_Q, s=f'$C$', fontsize=7, c='k', horizontalalignment='right', verticalalignment='bottom')
        ax.text(x=KYBER_Q/2, y=X, s=f'$D$', fontsize=7, c='k', horizontalalignment='right', verticalalignment='bottom')
        ax.text(x=X, y=-KYBER_Q/2, s=f'$E$', fontsize=7, c='k', horizontalalignment='right', verticalalignment='bottom')
        ax.text(x=KYBER_Q/2+X, y=-KYBER_Q/2+X, s=f'$F$', fontsize=7, c='k', horizontalalignment='right', verticalalignment='bottom')

        plt.scatter([0], [0], marker='.', s=4, c='k')
        plt.scatter([KYBER_Q/2], [X], marker='.', s=4, c='k')
        plt.scatter([X], [-KYBER_Q/2], marker='.', s=4, c='k')
        plt.scatter([KYBER_Q/2+X], [-KYBER_Q/2+X], marker='.', s=4, c='k')
        plt.scatter([KYBER_Q/2+X], [KYBER_Q/2+X], marker='.', s=4, c='k')
        plt.scatter([KYBER_Q/2+X - KYBER_Q], [KYBER_Q/2+X - KYBER_Q], marker='.', s=4, c='k')

        # xlabelpad = ax.xaxis.labelpad
        # ylabelpad = ax.yaxis.labelpad
        # for val, letter in [(-KYBER_Q/2, '-q/2'), (KYBER_Q/2, 'q/2')]:
        #     ax.annotate('${}$'.format(letter), xy=(val, 0), xytext=(0, +xlabelpad),
        #                 # xycoords=('data', 'axes fraction'),
        #                 textcoords='offset points',
        #                 ha='left' if val > 0 else 'right', va='center')
        #     # ax.annotate('${}$'.format(letter), xy=(0, val), xytext=(-ylabelpad, 0),
        #     #             # xycoords=('data', 'axes fraction'),
        #     #             textcoords='offset points',
        #     #             ha='right', va='top' if val > 0 else 'bottom')

        # rect = plt.Rectangle((-KYBER_Q/2, -KYBER_Q/2), KYBER_Q, KYBER_Q, fill=None, edgecolor='k', lw=0.3)
        triangle = plt.Polygon([(-KYBER_Q/2, -KYBER_Q/2), (KYBER_Q/2, -KYBER_Q/2), (KYBER_Q/2, KYBER_Q/2)], facecolor=(0,0,0,0.1), edgecolor='k', lw=0)
        # ax.add_patch(rect)
        ax.add_patch(triangle)

        plt.axhline(y=-KYBER_Q/2, lw=0.3, c='k', ls=':')
        plt.axhline(y=KYBER_Q/2, lw=0.3, c='k', ls=':')
        plt.axvline(x=-KYBER_Q/2, lw=0.3, c='k', ls=':')
        plt.axvline(x=KYBER_Q/2, lw=0.3, c='k', ls=':')


        decoding_lines = CodePlots.get_line_equations_for_good_code(X)


           # ax.annotate('${}$'.format(letter), xy=(val, 0), xytext=(0, -xlabelpad),
           #              # xycoords=('data', 'axes fraction'),
           #              textcoords='offset points',
           #              ha='right' if val > 0 else 'left', va='top')
           #  ax.annotate('${}$'.format(letter), xy=(0, val), xytext=(-ylabelpad, 0),
           #              # xycoords=('data', 'axes fraction'),
           #              textcoords='offset points',
           #              ha='right', va='top' if val > 0 else 'bottom')
        def put_label_if_in_limit(i, x_text, y_text, get_x_from_y, get_y_from_x):
            if i == 0:
                return

            if (y_text <= -KYBER_Q/2 and x_text >= -KYBER_Q/2):
                pady = 300
                padx = 50
                # return ax.text(x=x_text, y=y_text, s=fr'$\ell_{i}$')
                ax.annotate(fr'$\ell_{i}$', xy=(x_text, y_text),
                            xytext=(padx + get_x_from_y(y_text - pady), y_text - pady)),

            if (x_text >= KYBER_Q/2):
                padx = 200
                pady = 100

                ax.annotate(fr'$\ell_{i}$', xy=(x_text, y_text),
                            xytext=(x_text + padx, get_y_from_x(x_text + padx) + pady)),
                # return ax.text(x=x_text, y=y_text, s=fr'$\ell_{i}$')
            return


        for i, l in enumerate(decoding_lines):
            if i == 0:
                continue
            a_x, a_y, a_c = l
            x0 = -KYBER_Q
            y0 = (-a_x * x0 - a_c)/a_y
            x1 = KYBER_Q
            y1 = (-a_x * x1 - a_c)/a_y
            plt.axline((x0, y0), (x1, y1), lw=0.4, c='k', ls='--')

            def get_x_from_y(y_text):
                return (-a_y * y_text - a_c) / a_x

            def get_y_from_x(x_text):
                return (-a_x * x_text - a_c) / a_y

            for y_text in [-KYBER_Q/2, KYBER_Q/2]:
                x_text = get_x_from_y(y_text)
                put_label_if_in_limit(i, x_text, y_text, get_x_from_y, get_y_from_x)

            for x_text in [-KYBER_Q/2, KYBER_Q/2]:
                y_text = get_y_from_x(x_text)
                put_label_if_in_limit(i, x_text, y_text, get_x_from_y, get_y_from_x)

        plt.tight_layout()
        plt.subplots_adjust(left=0.2, right=0.93, top=0.93, bottom=0.1)

    def draw_good_2d_code_symmetry(X=422, figsize=(TCHES_HALF_PAGE_FIG, TCHES_HALF_PAGE_FIG)):

        fig, ax = plt.subplots(1, figsize=figsize)

        base = np.array([
            np.array((KYBER_Q//2, X)),
            np.array((X, KYBER_Q//2)),
        ])

        base_kyber = np.array([
            np.array((0, KYBER_Q)),
            np.array((KYBER_Q, 0)),
        ])
        points = []
        points_idx = []

        for ik, jk in it.product(range(-3, 3), repeat=2):
            kyber_point = np.matmul(base_kyber, np.array([ik, jk]))
            for i, j in it.product(range(0, 2), repeat=2):
                points_idx.append((i % 2, j % 2))
                points.append(kyber_point + np.matmul(base, np.array([i, j])))

        xs, ys = zip(*points)
        # plt.scatter(xs, ys, marker='.', s=2, c='k')
        ax = plt.gca()
        ax.set_aspect('equal')
        min_dist = CodePlots.get_minimum_distance(points)
        print(min_dist)

        corepoints = set()
        for ik, jk in it.product(range(2), repeat=2):
            corepoints.add(tuple(np.matmul(base, np.array([ik, jk]))))

        ghost_point_color = '0.5'
        for p in points:
            x, y = p
            # if tuple(p) in corepoints:
            # plt.scatter([x], [y], marker='.', s=4, c='k')
            # circle = plt.Circle(p, min_dist/2, facecolor=None, fill=False,
            #                     edgecolor='k', lw=0.7, ls='--')
            # else:
            plt.scatter([x], [y], marker='.', s=4, c=ghost_point_color)
            circle = plt.Circle(p, min_dist/2, facecolor=None, fill=False,
                                # edgecolor='k', lw=0.7, ls='--')
                                edgecolor=ghost_point_color, lw=0.3, ls='--')
            ax.add_patch(circle)

        ax.spines['left'].set_position(('data', 0))
        ax.spines['bottom'].set_position(('data', 0))

        sns.despine()
        # plt.xticks([i*KYBER_Q for i in range(-5, 5)])
        # plt.yticks([i*KYBER_Q for i in range(-5, 5)])
        plt.xticks([])
        plt.yticks([])
        plt.xlim(-KYBER_Q/2 - 500, KYBER_Q/2 + 500)
        plt.ylim(-KYBER_Q/2 -500, KYBER_Q/2 + 500)

        def in_limit(point):
            xlim, ylim = plt.xlim(), plt.ylim()

            if not (xlim[0] - 1 <= point[0] < xlim[1]):
                return False
            if not (ylim[0] - 1<= point[1] < ylim[1]):
                return False
            return True


        for (ik, jk), point in zip(points_idx, points):
            if not in_limit(point):
                continue

            # if tuple(point) in corepoints:
            #     ax.text(x=point[0] + 30, y=point[1] + 30, s=f'{ik}{jk}')
            # else:
            ax.text(x=point[0] + 30, y=point[1] + 30, s=f'{ik}{jk}', fontsize=7, c=ghost_point_color)


        xlabelpad = ax.xaxis.labelpad
        ylabelpad = ax.yaxis.labelpad
        for val, letter in [(-KYBER_Q/2, '-q/2'), (KYBER_Q/2, 'q/2')]:
            ax.annotate('${}$'.format(letter), xy=(val, 0), xytext=(0, +xlabelpad),
                        # xycoords=('data', 'axes fraction'),
                        textcoords='offset points',
                        ha='left' if val > 0 else 'right', va='center')
            ax.annotate('${}$'.format(letter), xy=(0, val), xytext=(-ylabelpad, 0),
                        # xycoords=('data', 'axes fraction'),
                        textcoords='offset points',
                        ha='right', va='bottom' if val > 0 else 'top')

        rect = plt.Rectangle((-KYBER_Q/2, -KYBER_Q/2), KYBER_Q, KYBER_Q, fill=None, edgecolor='k', lw=0.3)
        # triangle = plt.Polygon([(-KYBER_Q/2, -KYBER_Q/2), (KYBER_Q/2, -KYBER_Q/2), (KYBER_Q/2, KYBER_Q/2)], facecolor=(0,0,0,0.1), edgecolor='k', lw=0.3)
        triangle = plt.Polygon([(-KYBER_Q/2, -KYBER_Q/2), (KYBER_Q/2, -KYBER_Q/2), (KYBER_Q/2, KYBER_Q/2)], fill=None, edgecolor='k', lw=0.3)
        ax.add_patch(rect)
        ax.add_patch(triangle)

        target = np.array([(-KYBER_Q/2 + 800), 500])
        target_reflected = np.array([500, (-KYBER_Q/2 + 800)])
        c_reflected = np.array([X, -KYBER_Q/2])
        c_decoded = np.array([-KYBER_Q/2, X])

        arrow_1 = patches.FancyArrowPatch(target, target_reflected, lw=0.3, color='k', arrowstyle="simple,head_width=1,head_length=2", connectionstyle="arc3,rad=.2")
        arrow_2 = patches.FancyArrowPatch(target_reflected, c_reflected, lw=0.3, color='k', arrowstyle="simple,head_width=1,head_length=2", connectionstyle="arc3,rad=.2")
        arrow_3 = patches.FancyArrowPatch(c_reflected, c_decoded, lw=0.3, color='k', arrowstyle="simple,head_width=1,head_length=2", connectionstyle="arc3,rad=-.2")

        ax.add_patch(arrow_1)
        ax.add_patch(arrow_2)
        ax.add_patch(arrow_3)
        ax.text(x=target[0] + 30, y=target[1] + 30, s=r'$\mathbf{x}$')
        plt.scatter(*target, marker='.', s=4, c='k')
        plt.scatter(*target_reflected, marker='.', s=4, c='k')
        plt.scatter(*c_reflected, marker='.', s=4, c='k')
        plt.scatter(*c_decoded, marker='.', s=4, c='k')

        plt.tight_layout()
        plt.subplots_adjust(left=0.2, right=0.93, top=0.93, bottom=0.1)


    def draw_good_2d_code_symmetry_voronoi(X=422, figsize=(TCHES_HALF_PAGE_FIG, TCHES_HALF_PAGE_FIG)):

        fig, ax = plt.subplots(1, figsize=figsize)

        base = np.array([
            np.array((KYBER_Q//2, X)),
            np.array((X, KYBER_Q//2)),
        ])

        base_kyber = np.array([
            np.array((0, KYBER_Q)),
            np.array((KYBER_Q, 0)),
        ])
        points = []
        points_idx = []

        for ik, jk in it.product(range(-3, 3), repeat=2):
            kyber_point = np.matmul(base_kyber, np.array([ik, jk]))
            for i, j in it.product(range(0, 2), repeat=2):
                points_idx.append((i % 2, j % 2))
                points.append(kyber_point + np.matmul(base, np.array([i, j])))

        vor = scipy.spatial.Voronoi(points)
        scipy.spatial.voronoi_plot_2d(vor, ax=ax, show_vertices=False, show_points=False, line_colors='0.4', line_width=0.3)

        xs, ys = zip(*points)
        # plt.scatter(xs, ys, marker='.', s=2, c='k')
        ax = plt.gca()
        ax.set_aspect('equal')
        min_dist = CodePlots.get_minimum_distance(points)
        print(min_dist)

        corepoints = set()
        for ik, jk in it.product(range(2), repeat=2):
            corepoints.add(tuple(np.matmul(base, np.array([ik, jk]))))

        ghost_point_color = '0.5'
        for p in points:
            x, y = p
            # if tuple(p) in corepoints:
            # plt.scatter([x], [y], marker='.', s=4, c='k')
            # circle = plt.Circle(p, min_dist/2, facecolor=None, fill=False,
            #                     edgecolor='k', lw=0.7, ls='--')
            # else:
            plt.scatter([x], [y], marker='.', s=4, c=ghost_point_color)
            circle = plt.Circle(p, min_dist/2, facecolor=None, fill=False,
                                # edgecolor='k', lw=0.7, ls='--')
                                edgecolor=ghost_point_color, lw=0.3, ls='--')
            ax.add_patch(circle)

        # ax.spines['left'].set_position(('data', 0))
        # ax.spines['bottom'].set_position(('data', 0))

        sns.despine()
        # plt.xticks([i*KYBER_Q for i in range(-5, 5)])
        # plt.yticks([i*KYBER_Q for i in range(-5, 5)])
        # plt.xticks([])
        # plt.yticks([])
        # plt.axhline(y=-KYBER_Q/2, lw=0.3, c='k', ls=':')
        # plt.axhline(y=KYBER_Q/2, lw=0.3, c='k', ls=':')
        # plt.axvline(x=-KYBER_Q/2, lw=0.3, c='k', ls=':')
        # plt.axvline(x=KYBER_Q/2, lw=0.3, c='k', ls=':')

        plt.xticks([-KYBER_Q/2, 0, KYBER_Q/2], ['$-q/2$', '$0$', '$q/2$'])
        plt.yticks([-KYBER_Q/2, 0, KYBER_Q/2], ['$-q/2$', '$0$', '$q/2$'])
        plt.xlim(-KYBER_Q/2 - 500, KYBER_Q/2 + 500)
        plt.ylim(-KYBER_Q/2 -500, KYBER_Q/2 + 500)

        def in_limit(point):
            xlim, ylim = plt.xlim(), plt.ylim()

            if not (xlim[0] - 1 <= point[0] < xlim[1]):
                return False
            if not (ylim[0] - 1<= point[1] < ylim[1]):
                return False
            return True


        for (ik, jk), point in zip(points_idx, points):
            if not in_limit(point):
                continue

            # if tuple(point) in corepoints:
            #     ax.text(x=point[0] + 30, y=point[1] + 30, s=f'{ik}{jk}')
            # else:
            ax.text(x=point[0] + 30, y=point[1] + 30, s=f'{ik}{jk}', fontsize=7, c=ghost_point_color)


        # xlabelpad = ax.xaxis.labelpad
        # ylabelpad = ax.yaxis.labelpad
        # for val, letter in [(-KYBER_Q/2, '-q/2'), (KYBER_Q/2, 'q/2')]:
        #     ax.annotate('${}$'.format(letter), xy=(val, 0), xytext=(0, +xlabelpad),
        #                 # xycoords=('data', 'axes fraction'),
        #                 textcoords='offset points',
        #                 ha='left' if val > 0 else 'right', va='center')
        #     ax.annotate('${}$'.format(letter), xy=(0, val), xytext=(-ylabelpad, 0),
        #                 # xycoords=('data', 'axes fraction'),
        #                 textcoords='offset points',
        #                 ha='right', va='bottom' if val > 0 else 'top')

        rect = plt.Rectangle((-KYBER_Q/2, -KYBER_Q/2), KYBER_Q, KYBER_Q, fill=None, edgecolor='k', lw=0.5, ls='--')
        # triangle = plt.Polygon([(-KYBER_Q/2, -KYBER_Q/2), (KYBER_Q/2, -KYBER_Q/2), (KYBER_Q/2, KYBER_Q/2)], facecolor=(0,0,0,0.1), edgecolor='k', lw=0.3)
        # triangle = plt.Polygon([(-KYBER_Q/2, -KYBER_Q/2), (KYBER_Q/2, -KYBER_Q/2), (KYBER_Q/2, KYBER_Q/2)], fill=None, edgecolor='k', lw=0.5, ls='--')
        ax.add_patch(rect)
        ax.plot((-KYBER_Q/2, KYBER_Q/2), (-KYBER_Q/2, KYBER_Q/2), lw=0.5, ls='--', color='k')
        # ax.add_patch(triangle)

        target = np.array([(-KYBER_Q/2 + 800), 500])
        target_reflected = np.array([500, (-KYBER_Q/2 + 800)])
        c_reflected = np.array([X, -KYBER_Q/2])
        c_decoded = np.array([-KYBER_Q/2, X])

        arrow_1 = patches.FancyArrowPatch(target, target_reflected, lw=0.3, color='k', arrowstyle="simple,head_width=1,head_length=2", connectionstyle="arc3,rad=.2")
        arrow_2 = patches.FancyArrowPatch(target_reflected, c_reflected, lw=0.3, color='k', arrowstyle="simple,head_width=1,head_length=2", connectionstyle="arc3,rad=.2")
        arrow_3 = patches.FancyArrowPatch(c_reflected, c_decoded, lw=0.3, color='k', arrowstyle="simple,head_width=1,head_length=2", connectionstyle="arc3,rad=-.2")

        ax.add_patch(arrow_1)
        ax.add_patch(arrow_2)
        ax.add_patch(arrow_3)
        ax.text(x=target[0] + 30, y=target[1] + 30, s=r'$\mathbf{x}$')
        plt.scatter(*target, marker='.', s=4, c='k')
        plt.scatter(*target_reflected, marker='.', s=4, c='k')
        plt.scatter(*c_reflected, marker='.', s=4, c='k')
        plt.scatter(*c_decoded, marker='.', s=4, c='k')

        plt.tight_layout()
        plt.subplots_adjust(left=0.2, right=0.93, top=0.93, bottom=0.1)


    @staticmethod
    def draw_kyber_code(figsize):
        CodePlots.draw_code(X=0, figsize=figsize)

    @staticmethod
    def draw_good_2d_code(figsize=(TCHES_HALF_PAGE_FIG, TCHES_HALF_PAGE_FIG)):
        CodePlots.draw_code(X=446, figsize=figsize)

    @staticmethod
    def draw_lattice_code(figsize):
        CodePlots.draw_code(X=446, figsize=figsize, lattice_code=True)


    @staticmethod
    def draw_codes(figures_dir):
        CodePlots.draw_kyber_code(figsize=(TCHES_FIG_WIDTH/3, TCHES_FIG_WIDTH/3 * 0.92))
        plt.savefig(os.path.join(figures_dir, 'codes/kyber_code_2cols.pdf'))

        CodePlots.draw_good_2d_code(figsize=(TCHES_FIG_WIDTH/3, TCHES_FIG_WIDTH/3 * 0.92))
        plt.savefig(os.path.join(figures_dir, 'codes/kyber_better_2d_code_2cols.pdf'))

        CodePlots.draw_lattice_code(figsize=(TCHES_FIG_WIDTH/3, TCHES_FIG_WIDTH/3 * 0.92))
        plt.savefig(os.path.join(figures_dir, 'codes/lattice_code_2cols.pdf'))

    def draw_code_decoding(figures_dir):
        CodePlots.draw_good_2d_code_decoding_lines(X=446, figsize=(TCHES_HALF_PAGE_FIG_SMALL, TCHES_HALF_PAGE_FIG_SMALL))
        plt.savefig(os.path.join(figures_dir, 'codes/kyber_code_decoding_2cols.pdf'))

        CodePlots.draw_good_2d_code_symmetry(X=446, figsize=(TCHES_HALF_PAGE_FIG_SMALL, TCHES_HALF_PAGE_FIG_SMALL))
        plt.savefig(os.path.join(figures_dir, 'codes/kyber_code_decoding_symmetry_2cols.pdf'))

    def draw_code_decoding_voronoi(figures_dir):
        CodePlots.draw_good_2d_code_decoding_lines_voronoi(X=446, figsize=(TCHES_FIG_WIDTH/2.5, TCHES_FIG_WIDTH/2.5 * 0.92))
        plt.savefig(os.path.join(figures_dir, 'codes/kyber_code_decoding_2cols.pdf'))

        CodePlots.draw_good_2d_code_symmetry_voronoi(X=446, figsize=(TCHES_FIG_WIDTH/2.5, TCHES_FIG_WIDTH/2.5 * 0.92))
        plt.savefig(os.path.join(figures_dir, 'codes/kyber_code_decoding_symmetry_2cols.pdf'))

    def draw_code_with_q_grid(X=422, Y=KYBER_Q//2, lattice_code=False, figsize=(TCHES_HALF_PAGE_FIG, TCHES_HALF_PAGE_FIG)):

        fig, ax = plt.subplots(1, figsize=figsize)


        base = np.array([
            np.array((Y, X)),
            np.array((X, Y)),
        ])

        base_kyber = np.array([
            np.array((0, KYBER_Q)),
            np.array((KYBER_Q, 0)),
        ])
        points = []
        points_idx = []

        if not lattice_code:
            for ik, jk in it.product(range(-3, 3), repeat=2):
                kyber_point = np.matmul(base_kyber, np.array([ik, jk]))
                for i, j in it.product(range(0, 2), repeat=2):
                    points_idx.append((i % 2, j % 2))
                    points.append(kyber_point + np.matmul(base, np.array([i, j])))
        else:
            for i, j in it.product(range(-6, 6), repeat=2):
                points_idx.append((i % 2, j % 2))
                points.append(np.matmul(base, np.array([i, j])))

        xs, ys = zip(*points)
        # plt.scatter(xs, ys, marker='.', s=2, c='k')
        ax = plt.gca()
        ax.set_aspect('equal')
        min_dist = CodePlots.get_minimum_distance(points)
        print(min_dist)

        corepoints = set()
        for ik, jk in it.product(range(2), repeat=2):
            corepoints.add(tuple(np.matmul(base, np.array([ik, jk]))))

        ghost_point_color = '0.5'
        for p in points:
            x, y = p
            if tuple(p) in corepoints:
                plt.scatter([x], [y], marker='.', s=4, c='k')
                circle = plt.Circle(p, min_dist/2, facecolor=None, fill=False,
                                    edgecolor='k', lw=0.3, ls='--')
            else:
                plt.scatter([x], [y], marker='.', s=4, c=ghost_point_color)
                circle = plt.Circle(p, min_dist/2, facecolor=None, fill=False,
                                    edgecolor='k', lw=0.3, ls='--')
                                    # edgecolor=ghost_point_color, lw=0.3, ls='--')
            ax.add_patch(circle)

        ax.spines['left'].set_position(('data', 0))
        ax.spines['bottom'].set_position(('data', 0))

        sns.despine()
        # plt.xticks([i*KYBER_Q for i in range(-5, 5)])
        # plt.yticks([i*KYBER_Q for i in range(-5, 5)])
        plt.xticks([])
        plt.yticks([])
        plt.xlim(-KYBER_Q/2, 5*KYBER_Q/2)
        plt.ylim(-KYBER_Q/2, 5*KYBER_Q/2)

        def in_limit(point):
            xlim, ylim = plt.xlim(), plt.ylim()

            if not (xlim[0] - 1 <= point[0] < xlim[1]):
                return False
            if not (ylim[0] - 1<= point[1] < ylim[1]):
                return False
            return True


        for (ik, jk), point in zip(points_idx, points):
            if not in_limit(point):
                continue

            if tuple(point) in corepoints:
                ax.text(x=point[0] + 30, y=point[1] + 30, s=f'{ik}{jk}',  fontsize=7)
            else:
                ax.text(x=point[0] + 30, y=point[1] + 30, s=f'{ik}{jk}', fontsize=7, c=ghost_point_color)


        xlabelpad = ax.xaxis.labelpad
        ylabelpad = ax.yaxis.labelpad
        for val, letter in [(-KYBER_Q, '-q'), (KYBER_Q, 'q'), (2*KYBER_Q, '2q')]:
            ax.annotate('${}$'.format(letter), xy=(val, 0), xytext=(0, -xlabelpad),
                        # xycoords=('data', 'axes fraction'),
                        textcoords='offset points',
                        ha='right' if val > 0 else 'left', va='top')
            ax.annotate('${}$'.format(letter), xy=(0, val), xytext=(-ylabelpad, 0),
                        # xycoords=('data', 'axes fraction'),
                        textcoords='offset points',
                        ha='right', va='top' if val > 0 else 'bottom')

        plt.axhline(y=KYBER_Q, lw=0.3, c='k', ls=':')
        plt.axhline(y=2*KYBER_Q, lw=0.3, c='k', ls=':')
        plt.axvline(x=KYBER_Q, lw=0.3, c='k', ls=':')
        plt.axvline(x=2*KYBER_Q, lw=0.3, c='k', ls=':')
        rect = plt.Rectangle((0, 0), KYBER_Q, KYBER_Q, facecolor=(0,0,0,0.1), edgecolor='k', lw=0.3)
        ax.add_patch(rect)
        # rect = plt.Rectangle((KYBER_Q, KYBER_Q), KYBER_Q, KYBER_Q, fill=None, edgecolor='k', lw=0.3)
        # ax.add_patch(rect)
        # rect = plt.Rectangle((0, KYBER_Q), KYBER_Q, KYBER_Q, fill=None, edgecolor='k', lw=0.3)
        # ax.add_patch(rect)
        # rect = plt.Rectangle((KYBER_Q, 0), KYBER_Q, KYBER_Q, fill=None, edgecolor='k', lw=0.3)
        # ax.add_patch(rect)

        plt.tight_layout()


    def draw_code_with_q_grid_voronoi(X=422, Y=KYBER_Q//2, lattice_code=False, figsize=(TCHES_HALF_PAGE_FIG, TCHES_HALF_PAGE_FIG)):

        fig, ax = plt.subplots(1, figsize=figsize)


        base = np.array([
            np.array((Y, X)),
            np.array((X, Y)),
        ])

        base_kyber = np.array([
            np.array((0, KYBER_Q)),
            np.array((KYBER_Q, 0)),
        ])
        points = []
        points_idx = []

        if not lattice_code:
            for ik, jk in it.product(range(-3, 3), repeat=2):
                kyber_point = np.matmul(base_kyber, np.array([ik, jk]))
                for i, j in it.product(range(0, 2), repeat=2):
                    points_idx.append((i % 2, j % 2))
                    points.append(kyber_point + np.matmul(base, np.array([i, j])))
        else:
            for i, j in it.product(range(-6, 6), repeat=2):
                points_idx.append((i % 2, j % 2))
                points.append(np.matmul(base, np.array([i, j])))

        vor = scipy.spatial.Voronoi(points)
        scipy.spatial.voronoi_plot_2d(vor, ax=ax, show_vertices=False, show_points=False, line_colors='k', line_width=0.3)

        xs, ys = zip(*points)
        # plt.scatter(xs, ys, marker='.', s=2, c='k')
        ax = plt.gca()
        ax.set_aspect('equal')
        min_dist = CodePlots.get_minimum_distance(points)
        print(min_dist)

        corepoints = set()
        for ik, jk in it.product(range(2), repeat=2):
            corepoints.add(tuple(np.matmul(base, np.array([ik, jk]))))

        ghost_point_color = '0.5'
        for p in points:
            x, y = p
            if tuple(p) in corepoints:
                plt.scatter([x], [y], marker='.', s=4, c='k')
                circle = plt.Circle(p, min_dist/2, facecolor=None, fill=False,
                                    edgecolor='k', lw=0.3, ls='--')
            else:
                plt.scatter([x], [y], marker='.', s=4, c=ghost_point_color)
                circle = plt.Circle(p, min_dist/2, facecolor=None, fill=False,
                                    edgecolor='k', lw=0.3, ls='--')
                                    # edgecolor=ghost_point_color, lw=0.3, ls='--')
            ax.add_patch(circle)

        # return
        # ax.spines['left'].set_position(('data', 0))
        # ax.spines['bottom'].set_position(('data', 0))

        sns.despine()
        # plt.xticks([i*KYBER_Q for i in range(-5, 5)])
        # plt.yticks([i*KYBER_Q for i in range(-5, 5)])
        plt.xticks([0, KYBER_Q, 2*KYBER_Q], ['$0$', '$q$', '$2q$'])
        plt.yticks([0, KYBER_Q, 2*KYBER_Q], ['$0$', '$q$', '$2q$'])
        plt.xlim(-KYBER_Q/2, 5*KYBER_Q/2)
        plt.ylim(-KYBER_Q/2, 5*KYBER_Q/2)

        def in_limit(point):
            xlim, ylim = plt.xlim(), plt.ylim()

            if not (xlim[0] - 1 <= point[0] < xlim[1]):
                return False
            if not (ylim[0] - 1<= point[1] < ylim[1]):
                return False
            return True


        for (ik, jk), point in zip(points_idx, points):
            if not in_limit(point):
                continue

            if tuple(point) in corepoints:
                ax.text(x=point[0] + 30, y=point[1] + 30, s=f'{ik}{jk}',  fontsize=7)
            else:
                ax.text(x=point[0] + 30, y=point[1] + 30, s=f'{ik}{jk}', fontsize=7, c=ghost_point_color)


        # xlabelpad = ax.xaxis.labelpad
        # ylabelpad = ax.yaxis.labelpad
        # for val, letter in [(-KYBER_Q, '-q'), (KYBER_Q, 'q'), (2*KYBER_Q, '2q')]:
        #     ax.annotate('${}$'.format(letter), xy=(val, 0), xytext=(0, -xlabelpad),
        #                 # xycoords=('data', 'axes fraction'),
        #                 textcoords='offset points',
        #                 ha='right' if val > 0 else 'left', va='top')
        #     ax.annotate('${}$'.format(letter), xy=(0, val),
        #                 # xytext=(-ylabelpad, 0),
        #                 # xytext=(-ylabelpad, 0),
        #                 # xycoords=('data', 'axes fraction'),
        #                 xycoords='axes points',
        #                 textcoords='offset points',
        #                 ha='right', va='top' if val > 0 else 'bottom')

        plt.axhline(y=0, lw=0.3, c='k', ls=':')
        plt.axhline(y=KYBER_Q, lw=0.3, c='k', ls=':')
        plt.axhline(y=2*KYBER_Q, lw=0.3, c='k', ls=':')
        plt.axvline(x=0, lw=0.3, c='k', ls=':')
        plt.axvline(x=KYBER_Q, lw=0.3, c='k', ls=':')
        plt.axvline(x=2*KYBER_Q, lw=0.3, c='k', ls=':')
        rect = plt.Rectangle((0, 0), KYBER_Q, KYBER_Q, facecolor=(0,0,0,0.1), edgecolor='k', lw=0.3)
        ax.add_patch(rect)
        # rect = plt.Rectangle((KYBER_Q, KYBER_Q), KYBER_Q, KYBER_Q, fill=None, edgecolor='k', lw=0.3)
        # ax.add_patch(rect)
        # rect = plt.Rectangle((0, KYBER_Q), KYBER_Q, KYBER_Q, fill=None, edgecolor='k', lw=0.3)
        # ax.add_patch(rect)
        # rect = plt.Rectangle((KYBER_Q, 0), KYBER_Q, KYBER_Q, fill=None, edgecolor='k', lw=0.3)
        # ax.add_patch(rect)

        plt.tight_layout()

    def draw_code_with_q_grid_voronoi_pnorm(X=422, Y=KYBER_Q//2, lattice_code=False, pnorm=2,
                                            figsize=(TCHES_HALF_PAGE_FIG, TCHES_HALF_PAGE_FIG)):

        fig, ax = plt.subplots(1, figsize=figsize)


        base = np.array([
            np.array((Y, X)),
            np.array((X, Y)),
        ])

        base_kyber = np.array([
            np.array((0, KYBER_Q)),
            np.array((KYBER_Q, 0)),
        ])
        points = []
        points_idx = []

        if not lattice_code:
            for ik, jk in it.product(range(-3, 3), repeat=2):
                kyber_point = np.matmul(base_kyber, np.array([ik, jk]))
                for i, j in it.product(range(0, 2), repeat=2):
                    points_idx.append((i % 2, j % 2))
                    points.append(kyber_point + np.matmul(base, np.array([i, j])))
        else:
            for i, j in it.product(range(-6, 6), repeat=2):
                points_idx.append((i % 2, j % 2))
                points.append(np.matmul(base, np.array([i, j])))

        # vor = scipy.spatial.Voronoi(points)
        # scipy.spatial.voronoi_plot_2d(vor, ax=ax, show_vertices=False, show_points=False, line_colors='k', line_width=0.3)

        xs, ys = zip(*points)
        # plt.scatter(xs, ys, marker='.', s=2, c='k')
        ax = plt.gca()
        ax.set_aspect('equal')
        min_dist = CodePlots.get_minimum_distance(points, pnorm=pnorm)
        print(min_dist)

        corepoints = set()
        for ik, jk in it.product(range(2), repeat=2):
            corepoints.add(tuple(np.matmul(base, np.array([ik, jk]))))

        ghost_point_color = '0.5'
        for p in points:
            x, y = p
            if tuple(p) in corepoints:
                plt.scatter([x], [y], marker='.', s=4, c='k')
            else:
                plt.scatter([x], [y], marker='.', s=4, c=ghost_point_color)

            CodePlots.plot_circle_pnorm(p, min_dist/2, p=pnorm, ax=ax, color='k', lw=0.3, ls='--')

            # circle = plt.Circle(p, min_dist/2, facecolor=None, fill=False,
            #                     edgecolor='k', lw=0.3, ls='--')
            # ax.add_patch(circle)



        # return
        # ax.spines['left'].set_position(('data', 0))
        # ax.spines['bottom'].set_position(('data', 0))

        sns.despine()
        # plt.xticks([i*KYBER_Q for i in range(-5, 5)])
        # plt.yticks([i*KYBER_Q for i in range(-5, 5)])
        plt.xticks([0, KYBER_Q, 2*KYBER_Q], ['$0$', '$q$', '$2q$'])
        plt.yticks([0, KYBER_Q, 2*KYBER_Q], ['$0$', '$q$', '$2q$'])
        plt.xlim(-0.6* KYBER_Q/2, 4.6*KYBER_Q/2)
        plt.ylim(-0.6* KYBER_Q/2, 4.6*KYBER_Q/2)

        def in_limit(point):
            xlim, ylim = plt.xlim(), plt.ylim()

            if not (xlim[0] - 1 <= point[0] < xlim[1]):
                return False
            if not (ylim[0] - 1<= point[1] < ylim[1]):
                return False
            return True


        for (ik, jk), point in zip(points_idx, points):
            if not in_limit(point):
                continue

            if tuple(point) in corepoints:
                ax.text(x=point[0] + 30, y=point[1] + 30, s=f'{ik}{jk}',  fontsize=7)
            else:
                ax.text(x=point[0] + 30, y=point[1] + 30, s=f'{ik}{jk}', fontsize=7, c=ghost_point_color)

        plt.axhline(y=0, lw=0.3, c='k', ls=':')
        plt.axhline(y=KYBER_Q, lw=0.3, c='k', ls=':')
        plt.axhline(y=2*KYBER_Q, lw=0.3, c='k', ls=':')
        plt.axvline(x=0, lw=0.3, c='k', ls=':')
        plt.axvline(x=KYBER_Q, lw=0.3, c='k', ls=':')
        plt.axvline(x=2*KYBER_Q, lw=0.3, c='k', ls=':')
        rect = plt.Rectangle((0, 0), KYBER_Q, KYBER_Q, facecolor=(0,0,0,0.1), edgecolor='k', lw=0.3)
        ax.add_patch(rect)

        plt.tight_layout()


    def plot_circle_pnorm(center, radius, p=2, ax=None, color='k', lw=0.3, ls='--'):
        # fig, ax = plt.subplots(1, figsize=(TCHES_HALF_PAGE_FIG, TCHES_HALF_PAGE_FIG_BIG))

        # if ax is None:
        #     fig, ax = plt.subplots(1)

        # circle = plt.Circle((0, 0), 10, facecolor=None, fill=True, edgecolor='r')

        # ax.set_aspect(1)

        # ax.add_patch(circle)
        # plt.xlim(-50, 50)
        # plt.ylim(-50, 50)

        center_x, center_y = center
        x = np.linspace(center_x - radius, center_x + radius, 500)
        y = np.linspace(center_y - radius, center_y + radius, 500)
        X, Y = np.meshgrid(x, y)

        # F = X**2 + Y**2 -9
        F = (np.abs(X - center_x)**p + np.abs(Y - center_y)**p - radius**p)
        ax.contour(X, Y,F,[0], colors=color, linewidths=lw, linestyles=ls)

        # F = X**2 + Y**2 -100
        # ax.contour(X,Y,F,[0])
        return ax


    def draw_pair_of_good_codes(figures_dir):
        X, Y = 200, 1800
        CodePlots.draw_code_with_q_grid(X=X, Y=Y, figsize=(TCHES_HALF_PAGE_FIG_SMALL, TCHES_HALF_PAGE_FIG_SMALL))
        plt.savefig(os.path.join(figures_dir, f'codes/candidate_code_{X}_{Y}.pdf'))
        # ../paper/figs/codes/candidate_code_200_1800.pdf

        X, Y = 446, 1664
        CodePlots.draw_code_with_q_grid(X=X, Y=Y, figsize=(TCHES_HALF_PAGE_FIG_SMALL, TCHES_HALF_PAGE_FIG_SMALL))
        plt.savefig(os.path.join(figures_dir, f'codes/candidate_code_{X}_{Y}.pdf'))
        # ../paper/figs/codes/candidate_code_446_1664.pdf


    def draw_pair_of_good_codes_voronoi(figures_dir):
        X, Y = 200, 1800
        CodePlots.draw_code_with_q_grid_voronoi(X=X, Y=Y, figsize=(TCHES_HALF_PAGE_FIG_SMALL, TCHES_HALF_PAGE_FIG_SMALL))
        plt.savefig(os.path.join(figures_dir, f'codes/candidate_code_{X}_{Y}.pdf'))
        # ../paper/figs/codes/candidate_code_200_1800.pdf

        X, Y = 446, 1664
        CodePlots.draw_code_with_q_grid_voronoi(X=X, Y=Y, figsize=(TCHES_HALF_PAGE_FIG_SMALL, TCHES_HALF_PAGE_FIG_SMALL))
        plt.savefig(os.path.join(figures_dir, f'codes/candidate_code_{X}_{Y}.pdf'))
        # ../paper/figs/codes/candidate_code_446_1664.pdf


    def draw_pair_of_good_codes_voronoi_pnorm(figures_dir):
        X, Y, pnorm = 446, 1664, 2
        CodePlots.draw_code_with_q_grid_voronoi_pnorm(X=X, Y=Y, pnorm=2, figsize=(TCHES_HALF_PAGE_FIG_SMALL, TCHES_HALF_PAGE_FIG_SMALL))
        plt.savefig(os.path.join(figures_dir, f'codes/candidate_code_pnorm{pnorm}_{X}_{Y}.pdf'))
        # ../paper/figs/codes/candidate_code_200_1800.pdf

        X, Y, pnorm = 339, 1664, 3
        CodePlots.draw_code_with_q_grid_voronoi_pnorm(X=X, Y=Y, pnorm=pnorm, figsize=(TCHES_HALF_PAGE_FIG_SMALL, TCHES_HALF_PAGE_FIG_SMALL))
        plt.savefig(os.path.join(figures_dir, f'codes/candidate_code_pnorm{pnorm}_{X}_{Y}.pdf'))

        # X, Y, pnorm = 264, 1664, 5
        X, Y, pnorm = 111, 1664, 10
        CodePlots.draw_code_with_q_grid_voronoi_pnorm(X=X, Y=Y, pnorm=pnorm, figsize=(TCHES_HALF_PAGE_FIG_SMALL, TCHES_HALF_PAGE_FIG_SMALL))
        plt.savefig(os.path.join(figures_dir, f'codes/candidate_code_pnorm{pnorm}_{X}_{Y}.pdf'))
        # ../paper/figs/codes/candidate_code_446_1664.pdf



class DistributionPlots:


    def select_kyber_parameters_from_df(df, sec_level):
        parameters = {
            1: {'q': 3329, 'level': 1, 'k': 2, 'eta1': 3, 'eta2': 2, 'n_blocks_for_du0': 2,
                'du0': 10, 'du1': 0,  'dv': 4},
            3: {'q': 3329, 'level': 3, 'k': 3, 'eta1': 2, 'eta2': 2, 'n_blocks_for_du0': 3,
                'du0': 10, 'du1': 0,  'dv': 4},
            5: {'q': 3329, 'level': 5, 'k': 4, 'eta1': 2, 'eta2': 2, 'n_blocks_for_du0': 4,
                'du0': 11, 'du1': 0,  'dv': 5},
        }
        return  df[
                    (df.KYBER_Q == parameters[sec_level]['q']) &
                    (df.KYBER_LEVEL == parameters[sec_level]['level']) &
                    (df.KYBER_K == parameters[sec_level]['k']) &
                    (df.KYBER_ETA1 == parameters[sec_level]['eta1']) &
                    (df.KYBER_ETA2 == parameters[sec_level]['eta2']) &
                    (df.KYBER_N_BLOCKS_FOR_DU0 == parameters[sec_level]['n_blocks_for_du0']) &
                    (df.KYBER_DU0 == parameters[sec_level]['du0']) &
                    (df.KYBER_DU1 == parameters[sec_level]['du1']) &
                    (df.KYBER_DV == parameters[sec_level]['dv'])
                ]

    def get_discrete_2d_level_curve2(dist, level_curve_prob=2**(-128)):
        min_value, max_value = min(dist), max(dist)

        level_curve = []
        for v1 in range(min_value, max_value + 1):

            for v2 in range(min_value, max_value + 1):
                prob = dist[v1] * dist[v2]

                if level_curve_prob/2 <= prob <= 2*level_curve_prob:
                    level_curve.append((v1, v2))
        return level_curve

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


    def plot_shape_of_distributions_different_dv(figures_dir):
        fig, axs = plt.subplots(1, 7, figsize=(TCHES_FIG_WIDTH, 0.7))

        for i, dv in enumerate(range(0, 7)):
        # for i, dv in enumerate(range(0, 1)):
            ax = axs[i]
            # p = min({p: np.var([np.linalg.norm(x, ord=p) for x in points]) for p in np.linspace(2, 3, 100)}.items(), key=lambda x: x[1])
            ks = KyberParameterSet(256, 4, 2, 2, 3329, 12, dv)
            error_dist = ks.get_final_error_distribution_1d()
            # xs, ys = zip(*error_dist.items())
            # plt.scatter(xs, ys)
            # print(len(DistributionPlots.get_discrete_2d_level_curve2(error_dist)))
            level_curve = DistributionPlots.get_discrete_2d_level_curve(error_dist)
            f = lambda p: np.var([np.linalg.norm(x, ord=p) for x in level_curve])
            p = utils.ternary_search(f, 2, 100, 1e-1, min_max='min')
            print(dv, p, len(level_curve))
            # ({p: np.var([np.linalg.norm(x, ord=p) for x in points]) for p in np.linspace(2, 3, 100)}.items(), key=lambda x: x[1])
            # print()
            xs, ys = zip(*level_curve)
            ax.set_aspect('equal', adjustable='box')
            ax.plot(xs, ys, lw=0.6, c='k')
            # ax.set_title(f'$d_v = {i}$\n$p = {p:.2f}$')
            # ax.set_title(f'dv = {i}')
            ax.axes.set_axis_off()
            sns.despine()

            t = ax.annotate(f'$d_v = {i}$\n$p = {p:.2f}$', (0.5, 0.5), xycoords='axes fraction', fontsize=6, ha='center', va='center')


        plt.tight_layout()
        plt.subplots_adjust(top=0.95, bottom=0.05, left=0.1, right=0.9, hspace=0.2, wspace=0.254)

        # plt.subplots_adjust(top=0.629, bottom=0.0, left=0.028, right=0.972, hspace=0.2, wspace=0.254)

        plt.savefig(os.path.join(figures_dir, 'codes/dist_shape_level5.pdf'))


class DFRPlots():

    def ct_size():
        def get_ciphertext_size(k, du, dv):
            return (256 // 8) * (k * du + dv)

    def get_ciphertext_size_ext(k, n_blocks_for_du0, du0, du1, dv):
        assert(n_blocks_for_du0 <= k)
        return (256 // 8) * (n_blocks_for_du0 * du0 + (k - n_blocks_for_du0) * du1 + dv)


    def plot_effect_of_dv(figures_dir, level=5, datafile='results/kyber2d_dfr/result_kyber2d_dfr_effect_of_dv_in_level5.csv'):
        df = pd.read_csv(datafile)

        fig, ax = plt.subplots(1, figsize=(TCHES_FIG_WIDTH, 2.2))

        styles = [
            {'ls': '--', 'c': '0.6'},
            {'ls': '-', 'c': '0.6'},
            {'ls': '-.',  'c': '0.4'},
            {'ls': '--',  'c': '0.2'},
            {'ls': '-',  'c': '0'},
        ]
        plt.plot([0], [0], marker='*', c='k', mfc='none', ms=8, lw=0, mew=0.5, label=r'Best code')


        df = df[(3 <= df.KYBER_DV) & (df.KYBER_DV < 8)]
        for i, (idx, g) in enumerate(df.groupby('KYBER_DV')):
            best_2d_code_row = g[g.DFR == g.DFR.min()].iloc[0]

            sns.lineplot(data=g, x='CODE_BETA', y='DFR', label=f'$d_v = {idx}$', lw=0.6, **styles[i])
            plt.plot([best_2d_code_row.CODE_BETA], [best_2d_code_row.DFR], marker='*', c='k', mfc='none', ms=8, mew=0.5, lw=0)
            # plt.plot([best_2d_code_row.CODE_BETA], [best_2d_code_row.DFR], marker='.', ms=1, c='k')

        plt.yscale('log', base=2)
        plt.xlabel(r'Code parameter $\beta$')
        plt.ylabel(r'Decryption failure rate')
        sns.despine()
        plt.ylim(2**-210, 2**-50)
        plt.xlim(-10, 510)

        # plt.plot([0], [0], marker='.', ms=10, c='black', lw=0, label='(d_u, d_v)')
        # plt.legend(loc='upper right')
        legend = ax.legend(fontsize=7, loc='upper right')
        legend.get_frame().set_linewidth(0.3)
        legend.get_frame().set_edgecolor("k")
        plt.tight_layout()
        plt.savefig(os.path.join(figures_dir, 'compression/effect_of_dv_in_level5.pdf'))


    def plot_effect_of_dimension_and_beta(figures_dir, datafile='results/minimum_distance/minimum_distance_p2.csv'):
        full_df = pd.read_csv(datafile)

        # sns.lineplot(data=df, x='beta', y='minimum_distance', hue='n')


        fig, ax = plt.subplots(1, figsize=(TCHES_FIG_WIDTH / 2, 2))

        styles = {
            2: {'ls': '-', 'c': '0.6'},
            3: {'ls': '--', 'c': '0.5'},
            4: {'ls': '-',  'c': '0.4'},
            5: {'ls': '--',  'c': '0.3'},
            6: {'ls': '-',  'c': '0.2'},
            7: {'ls': '--',  'c': '0.1'},
            8: {'ls': '-',  'c': '0'},
        }

        for n in range(2, 9):
            df = full_df[full_df.n == n][:-3]

            sns.lineplot(data=df, x='beta', y='minimum_distance', ax=ax, **styles[n], lw=0.3)

            y_peak = max(df.minimum_distance)
            x_peak = df[df.minimum_distance == y_peak].iloc[0].beta

            x_peak = df.iloc[-1].beta
            y_peak = df.iloc[-1].minimum_distance

            t = plt.annotate(rf'$\mu = {n}$', (x_peak, y_peak - 5), fontsize=6, ha='left')

            best_2d_code_row = df[df.minimum_distance == df.minimum_distance.max()].iloc[0]
            plt.plot([best_2d_code_row.beta], [best_2d_code_row.minimum_distance], marker='*', c='k', mfc='none', ms=8, mew=0.5, lw=0)

            print(df)

        plt.plot([0], [0], marker='*', c='k', mfc='none', ms=8, lw=0, mew=0.5, label=r'Max minimum distance')

        legend = ax.legend(fontsize=7, loc='upper left')
        legend.get_frame().set_linewidth(0.3)
        legend.get_frame().set_edgecolor("k")

        plt.xlim(390, 1050)
        plt.ylim(1680, 1980)

        sns.despine()
        plt.xlabel(r'Code parameter $\beta$')
        plt.ylabel(r'Minimum distance')
        plt.tight_layout()
        plt.subplots_adjust(left=0.19, right=0.95, top=0.96, bottom=0.2)
        plt.savefig(os.path.join(figures_dir, 'compression/effect_of_dimension_and_beta.pdf'))


    def plot_effect_of_dv_dense(figures_dir, level=5, datafile='results/kyber2d_dfr/result_kyber2d_dfr_effect_of_dv.csv'):
        df = pd.read_csv(datafile)

        df = df[(df.KYBER_LEVEL == 5) &
                (df.KYBER_Q == 3329) &
                (df.KYBER_ETA1 == 2) &
                (df.KYBER_ETA2 == 2) &
                (df.KYBER_DU0 == 11)]

        fig, ax = plt.subplots(1, figsize=(TCHES_FIG_WIDTH / 2, 2))

        styles = [
            {'ls': '--', 'c': '0.6'},
            {'ls': '-', 'c': '0.6'},
            {'ls': '-.',  'c': '0.4'},
            {'ls': '--',  'c': '0.2'},
            {'ls': '-',  'c': '0'},
            {'ls': '-',  'c': '0'},
            {'ls': '-',  'c': '0'},
            {'ls': '-',  'c': '0'},
        ]
        plt.plot([0], [0], marker='*', c='k', mfc='none', ms=8, lw=0, mew=0.5, label=r'Best 2D code')

        adjust_for_dv = {
            # 7: (0, 1/4),
            8: (0, 1/4),
        }

        df = df[(4 <= df.KYBER_DV) & (df.KYBER_DV < 8)]
        for i, (dv, g) in enumerate(df.groupby('KYBER_DV')):
            best_2d_code_row = g[g.DFR == g.DFR.min()].iloc[0]

            # sns.lineplot(data=g, x='CODE_BETA', y='DFR', label=f'$d_v = {dv}$', lw=0.6, **styles[i])
            g = g[g.CODE_BETA <= best_2d_code_row.CODE_BETA + 30]
            sns.lineplot(data=g, x='CODE_BETA', y='DFR', label=f'', lw=0.6, **styles[i])
            plt.plot([best_2d_code_row.CODE_BETA], [best_2d_code_row.DFR], marker='*', c='k', mfc='none', ms=8, mew=0.5, lw=0)
            # plt.plot([best_2d_code_row.CODE_BETA], [best_2d_code_row.DFR], marker='.', ms=1, c='k')

            u = g[g.DFR < 2**-150]
            # last_row = u[u.CODE_BETA == u.CODE_BETA.max()].iloc[0]
            # (x, y) = last_row.CODE_BETA, last_row.DFR
            label_row = u[u.CODE_BETA == 0].iloc[0]
            (x, y) = label_row.CODE_BETA, label_row.DFR

            if dv in adjust_for_dv:
                x += adjust_for_dv[dv][0]
                y *= adjust_for_dv[dv][1]

            # t = plt.annotate(rf'$d_v = {dv}$', (x + 10, y), fontsize=7)
            t = plt.annotate(rf'$d_v = {dv}$', (0, y*2), fontsize=6)


        plt.yscale('log', base=2)
        plt.xlabel(r'Code parameter $\beta$')
        plt.ylabel(r'DFR')
        sns.despine()
        plt.ylim(2**-206, 2**-150)
        # plt.ylim(-206, -150)
        plt.xlim(-10, 465)

        # plt.plot([0], [0], marker='.', ms=10, c='black', lw=0, label='(d_u, d_v)')
        # plt.legend(loc='upper right')
        legend = ax.legend(fontsize=7, loc='lower left')
        legend.get_frame().set_linewidth(0.3)
        legend.get_frame().set_edgecolor("k")
        plt.tight_layout()
        plt.subplots_adjust(left=0.21, right=0.98, top=0.96, bottom=0.2)

        plt.savefig(os.path.join(figures_dir, 'compression/effect_of_dv_in_level5.pdf'))



    def plot_compression_simple_factors_2d_4d(figures_dir,
                                              datafile='results/kyber2d_dfr/result_kyber2d_dfr.csv',
                                              datafile4d='results/kyber4d_dfr/result_kyber4d_dfr.csv'):

        df2d = pd.read_csv(datafile)
        df4d = pd.read_csv(datafile4d)


        # Filter data frame for only simple compression parameters (du, dv)
        df2d = df2d[df2d.KYBER_K == df2d.KYBER_N_BLOCKS_FOR_DU0]

        df2d['CODE_DIMENSION'] = 2

        df = pd.concat([df2d, df4d])
        # plt.figure(figsize=(TCHES_FIG_WIDTH, 3))
        fig, ax = plt.subplots(1, figsize=(TCHES_FIG_WIDTH, 2.7))

        print(df)

        def get_point(row):
            ct_size = DFRPlots.get_ciphertext_size_ext(k=row.KYBER_K,
                                                       n_blocks_for_du0=row.KYBER_N_BLOCKS_FOR_DU0,
                                                       du0=row.KYBER_DU0,
                                                       du1=row.KYBER_DU1,
                                                       dv=row.KYBER_DV)
            dfr = row.DFR
            return (ct_size, dfr)

        def get_point_label(row):
            vec_du = ([round(row.KYBER_DU0)] *  round(row.KYBER_N_BLOCKS_FOR_DU0) +
                      [round(row.KYBER_DU1)] * round(row.KYBER_K - row.KYBER_N_BLOCKS_FOR_DU0))
            vec_du_as_str = r','.join((str(s) for s in vec_du))
            if (row.KYBER_ETA2 == 0):
                label = rf'${{({vec_du[0]},{round(row.KYBER_DV)})}}^\star$'
            elif (row.KYBER_ETA2 == 1):
                label = rf'${{({vec_du[0]},{round(row.KYBER_DV)})}}^\dagger$'
            else:
                label = f'$({vec_du[0]},{round(row.KYBER_DV)})$'
            return (round(row.KYBER_LEVEL), label)


        labels_to_points = {}
        for idx, g in df.groupby(['KYBER_LEVEL', 'KYBER_N_BLOCKS_FOR_DU0',  'KYBER_DU0',  'KYBER_DU1', 'KYBER_DV', 'KYBER_ETA2']):
            # sns.lineplot(data=g, x='CODE_BETA', y='DFR', label=idx)
            # plt.yscale('log', base=2)

            g2 = g[g.CODE_DIMENSION == 2]
            try:
                best_2d_code_row = g2[(g2.DFR == g2.DFR.min())].iloc[0]
            except IndexError:
                continue

            g4 = g[g.CODE_DIMENSION == 4]
            best_4d_code_row = g4[(g4.DFR == g4.DFR.min())].iloc[0]

            kyber_code_row = g[g.CODE_BETA == 0].iloc[0]

            best_4d_code_point = get_point(best_4d_code_row)
            best_2d_code_point = get_point(best_2d_code_row)
            kyber_code_point = get_point(kyber_code_row)



            # plt.annotate(get_point_label(best_2d_code_row),
            #              best_2d_code_point, fontsize=6)
            plt.plot([best_2d_code_point[0]], [best_2d_code_point[1]], marker='P', mfc='k', mew=0.0, ms=3, c='k')
            plt.plot([best_4d_code_point[0]], [best_4d_code_point[1]], marker='x', mfc='none', mew=0.6, ms=4, c='k')

            # plt.annotate(get_point_label(kyber_code_row),
            #              kyber_code_point, fontsize=6)
            plt.plot([kyber_code_point[0]], [kyber_code_point[1]], marker='.', ms=2, c='k')

            l_4d = get_point_label(best_4d_code_row)
            l_2d = get_point_label(best_2d_code_row)
            l_kyber = get_point_label(kyber_code_row)
            if l_4d not in labels_to_points:
                labels_to_points[l_2d] = []
            if l_2d not in labels_to_points:
                labels_to_points[l_2d] = []
            if l_kyber not in labels_to_points:
                labels_to_points[l_kyber] = []
            labels_to_points[l_4d].append(best_4d_code_point)
            labels_to_points[l_2d].append(best_2d_code_point)
            labels_to_points[l_kyber].append(kyber_code_point)


        label_corrections = {
            (1, '$(9,3)$'): (0, 0),
            (1, '$(9,4)$'): (0, 0),
            (1, '$(9,5)$'): (0, 0),
            (1, '$(10,3)$'): (0, 0),
            (1, '$(10,4)$'): (0, 0),
            (1, r'${(9,5)}^\star$'): (0, 0),
            (3, '$(9,3)$'): (0, 0),
            (3, '$(9,4)$'): (0, 0),
            (3, '$(9,5)$'): (0, 0),
            (3, r'${(9,5)}^\star$'): (0, -12),
            (3, r'${(10,3)}^\dagger$'): (0, -1),
            (3, '$(10,3)$'): (0, -2),
            (3, '$(10,4)$'): (0, 0),
            (5, '$(10,4)$'): (0, 6),
            (5, '$(10,5)$'): (0, 12),
            (5, '$(10,6)$'): (0, 13),
            (5, '$(11,4)$'): (0, -6),
            (5, '$(11,5)$'): (0, 2),
        }
        marked = {
            (5, '$(11,5)$'),
            (5, '$(10,6)$'),
            (3, '$(10,4)$'),
            (1, '$(10,4)$'),
        }
        print(labels_to_points)
        for (level, l), points in labels_to_points.items():
            if (level, l) not in label_corrections:
                label_corrections[level, l] = (0, 0)
            x = points[0][0] + 30 + label_corrections[level, l][0]
            y = 2**np.mean([np.log2(p[1]) for p in points])*2**label_corrections[level, l][1]

            # if (level, l) in marked:
            #     t = plt.annotate(r'\underline{' + l + r'}', (x, y), fontsize=7)
            # else:
            #     plt.annotate(l, (x, y), fontsize=7, c='0.3')
            plt.annotate(l, (x, y), fontsize=7)

            print(level, l)
            for p in points:
                plt.annotate(xy=p, xytext=(x,y), arrowprops=dict(arrowstyle='-', linewidth=0.3), text='', )

        plt.yscale('log', base=2)
        sns.despine()
        # plt.xticks(range(1440, 1440 + 7*32, 32))
        # plt.xlim(1440 - 5, 1440 + 4*32 + 31)
        plt.xlabel('Ciphertext size (bytes)')
        plt.ylabel('Decryption failure rate')
        ax.annotate(r'$\textrm{DFR}_1$', xy=(1664, 2**-120), xytext=(0, 0.00), fontsize=8, textcoords='offset points', ha='center', va='bottom')
        ax.annotate(r'$\textrm{DFR}_3$', xy=(1664, 2**-136), xytext=(0, 0.00), fontsize=8, textcoords='offset points', ha='center', va='bottom')
        ax.annotate(r'$\textrm{DFR}_5$', xy=(1664, 2**-160), xytext=(0, 0.00), fontsize=8, textcoords='offset points', ha='center', va='bottom')
        # ax.annotate('Saber Level 3', xy=(0, 2**-136), xytext=(1664, 1.00), fontsize=8, xycoords='axes fraction', ha='center', va='bottom')
        # ax.annotate('Saber Level 5', xy=(0, 2**-165), xytext=(1664, 1.00), fontsize=8, xycoords='axes fraction', ha='center', va='bottom')


        # plt.axhline(y=2**-128, lw=0.3, c='k', ls='--')
        # plt.axhline(y=2**-160, lw=0.3, c='k', ls='--')

        plt.axhline(y=2**-120, lw=0.3, c='k', ls='--')
        plt.axhline(y=2**-136, lw=0.3, c='k', ls='--')
        plt.axhline(y=2**-160, lw=0.3, c='k', ls='--')


        plt.xticks(range(672 -32, 1664 + 64, 64), rotation=30)
        xlims = plt.xlim(672, 1664)
        first_level_line_x = 900
        second_level_line_x = 1300

        total_x = (xlims[1] - xlims[0])
        level_1_x = (first_level_line_x - xlims[0])/2 / total_x

        level_3_x = (first_level_line_x - xlims[0] + (second_level_line_x - first_level_line_x)/2) / total_x
        level_5_x = (second_level_line_x - xlims[0] + (xlims[1] - second_level_line_x)/2) / total_x

        ax.annotate('Level 1', xy=(0, 0), xytext=(level_1_x, 1.00), fontsize=8, xycoords='axes fraction', ha='center', va='bottom')
        ax.annotate('Level 3', xy=(0, 0), xytext=(level_3_x, 1.00), fontsize=8, xycoords='axes fraction', ha='center', va='bottom')
        ax.annotate('Level 5', xy=(0, 0), xytext=(level_5_x, 1.00), fontsize=8, xycoords='axes fraction', ha='center', va='bottom')

        print(f'{level_1_x, level_3_x, level_5_x=}')
        plt.axvline(x=first_level_line_x, lw=0.3, c='0.5', ls=':')
        plt.axvline(x=second_level_line_x, lw=0.3, c='0.5', ls=':')

        # plt.plot([0], [1], marker='.', ms=1, c='k', lw=0, label=r'$(d_\mathbf{u}, d_v)$ Kyber code')
        # plt.plot([0], [0], marker='*', c='k', mfc='none', ms=5, lw=0, mew=0.5, label=r'$(d_\mathbf{u}, d_v)$ Best 2D code')

        plt.plot([0], [1], marker='o', mfc='none', mew=0.5, ms=5, c='k', lw=0, label=r'ML-KEM')
        plt.plot([0], [1], marker='.', ms=2, c='k', lw=0, label=r'1D code')
        plt.plot([0], [0], marker='P', mfc='k', mew=0.0, ms=3, c='k', lw=0, label=r'2D Minal code')
        plt.plot([0], [0], marker='x', mfc='none', mew=0.6, ms=4, c='k', lw=0, label=r'4D Minal code')

        legend = ax.legend(fontsize=7, loc='lower left')
        legend.get_frame().set_linewidth(0.3)
        legend.get_frame().set_edgecolor("k")
        # plt.ylim(2**-192, 2**-70)
        # Careful: these were carefully selected to show lines 2**-128 and 2**160
        # plt.ylim(2**-208, 2**-54)
        plt.ylim(2**-208, 2**-72)


        kyber_params = [
            (768, 2**(-139.1)),
            (1088, 2**(-165.2)),
            (1568, 2**(-175.2)),
        ]

        for x, y in kyber_params:
            # sns.scatterplot(x=[x], y=[y], marker='*', mf=None, s=20)
            plt.plot([x], [y], marker='o', mfc='none', mew=0.5, ms=4, c='k')

        # x, y = kyber_params[0]
        # plt.annotate('Kyber512', xy=(x, y), xytext=(x - 150, y*2**-10), fontsize=6,
        #              xycoords='data', arrowprops=dict(arrowstyle='-', linewidth=0.3, color='0.3'), ha='left')

        # x, y = kyber_params[1]
        # plt.annotate('Kyber768', xy=(x, y), xytext=(x - 150, y*2**-10), fontsize=6,
        #              xycoords='data', arrowprops=dict(arrowstyle='-', linewidth=0.3, color='0.3'), ha='left')

        # x, y = kyber_params[2]
        # plt.annotate('Kyber1024', xy=(x, y), xytext=(x - 150, y*2**-10), fontsize=6,
        #              xycoords='data', arrowprops=dict(arrowstyle='-', linewidth=0.3, color='0.3'), ha='left')

        plt.tight_layout()
        plt.savefig(os.path.join(figures_dir, 'compression/balancing_du_and_dv.pdf'))



    def get_relevant_results(datafile='results/kyber2d_dfr/result_kyber2d_dfr.csv'):

        relevant_results_parameters = [
            {'level': 1, 'q': 3329, 'du0': 10, 'n_blocks_for_du0': 2, 'du1': 0, 'dv': 4},
            {'level': 3, 'q': 3329, 'du0': 10, 'n_blocks_for_du0': 3, 'du1': 0, 'dv': 4},
            {'level': 5, 'q': 3329, 'du0': 11, 'n_blocks_for_du0': 4, 'du1': 0, 'dv': 5},
            {'level': 5, 'q': 3329, 'du0': 10, 'n_blocks_for_du0': 4, 'du1': 0, 'dv': 6},
            {'level': 5, 'q': 3329, 'du0': 11, 'n_blocks_for_du0': 2, 'du1': 10, 'dv': 6},
        ]

        df = pd.read_csv(datafile)
        for params in relevant_results_parameters:
            filtered = df[(df.KYBER_LEVEL == params['level']) &
                          (df.KYBER_Q == params['q']) &
                          (df.KYBER_DU0 == params['du0']) &
                          (df.KYBER_N_BLOCKS_FOR_DU0 == params['n_blocks_for_du0']) &
                          (df.KYBER_DU1 == params['du1']) &
                          (df.KYBER_DV == params['dv'])]

            assert len(filtered) > 0
            min_dfr = filtered.DFR.min()
            min_dfr_row = filtered[filtered.DFR == min_dfr].iloc[0]
            results = params.copy()
            results['code_alpha'] = min_dfr_row.CODE_ALPHA
            results['code_beta'] = min_dfr_row.CODE_BETA
            results['dfr'] = min_dfr_row.DFR
            results['log2_dfr'] = math.log2(min_dfr_row.DFR)
            del results['q']
            del results['dfr']
            print(results)

    def get_results_below_dfr_targets(datafile='results/kyber2d_dfr/result_kyber2d_dfr.csv',
                                      datafile4d='results/kyber4d_dfr/result_kyber4d_dfr.csv'):

        df2d = pd.read_csv(datafile)
        df4d = pd.read_csv(datafile4d)


        # Filter data frame for only simple compression parameters (du, dv)
        df2d = df2d[df2d.KYBER_K == df2d.KYBER_N_BLOCKS_FOR_DU0]

        df2d['CODE_DIMENSION'] = 2

        df = pd.concat([df2d, df4d])

        df = df[df.KYBER_Q == 3329]


        DFR_TARGETS = {
            1: 2**(-120),
            3: 2**(-136),
            5: 2**(-160),
        }
        for keys, g in df.groupby(['KYBER_LEVEL', 'KYBER_DU0', 'KYBER_DV', 'KYBER_ETA2']):
            if min(g.DFR) >= DFR_TARGETS[keys[0]]:
                continue

            dfr_1d = min(g[(g.CODE_DIMENSION == 2) & (g.CODE_BETA == 0)].DFR)
            dfr_2d = min(g[(g.CODE_DIMENSION == 2)].DFR)
            dfr_4d = min(g[(g.CODE_DIMENSION == 4)].DFR)

            best_beta2d = g[(g.CODE_DIMENSION == 2) & (g.DFR <= dfr_2d)].iloc[0].CODE_BETA
            best_beta4d = g[(g.CODE_DIMENSION == 4) & (g.DFR <= dfr_4d)].iloc[0].CODE_BETA

            print('KYBER_LEVEL, KYBER_DU0, KYBER_DV, KYBER_ETA2 = ', list(map(int, keys)))
            print('1D DFR = $2^{%.1f}$' % np.log2(dfr_1d))
            print('2D DFR = $2^{%.1f}$' % np.log2(dfr_2d), 'beta = %d' % best_beta2d)
            print('4D DFR = $2^{%.1f}$' % np.log2(dfr_4d), 'beta = %d' % best_beta4d)



        return df


        relevant_results_parameters = [
            {'level': 1, 'q': 3329, 'du0': 10, 'n_blocks_for_du0': 2, 'du1': 0, 'dv': 4},
            {'level': 3, 'q': 3329, 'du0': 10, 'n_blocks_for_du0': 3, 'du1': 0, 'dv': 4},
            {'level': 5, 'q': 3329, 'du0': 11, 'n_blocks_for_du0': 4, 'du1': 0, 'dv': 5},
            {'level': 5, 'q': 3329, 'du0': 10, 'n_blocks_for_du0': 4, 'du1': 0, 'dv': 6},
            {'level': 5, 'q': 3329, 'du0': 11, 'n_blocks_for_du0': 2, 'du1': 10, 'dv': 6},
        ]


        df = pd.read_csv(datafile)
        for params in relevant_results_parameters:
            filtered = df[(df.KYBER_LEVEL == params['level']) &
                          (df.KYBER_Q == params['q']) &
                          (df.KYBER_DU0 == params['du0']) &
                          (df.KYBER_N_BLOCKS_FOR_DU0 == params['n_blocks_for_du0']) &
                          (df.KYBER_DU1 == params['du1']) &
                          (df.KYBER_DV == params['dv'])]

            assert len(filtered) > 0
            min_dfr = filtered.DFR.min()
            min_dfr_row = filtered[filtered.DFR == min_dfr].iloc[0]
            results = params.copy()
            results['code_alpha'] = min_dfr_row.CODE_ALPHA
            results['code_beta'] = min_dfr_row.CODE_BETA
            results['dfr'] = min_dfr_row.DFR
            results['log2_dfr'] = math.log2(min_dfr_row.DFR)
            del results['q']
            del results['dfr']
            print(results)


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

if __name__ == '__main__':
    main()
