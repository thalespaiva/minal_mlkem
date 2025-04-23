#------------------------------------------------------------------------------
# This file was adapted from the implementation
# (available at: Public Domain https://github.com/pq-crystals/kyber)
# of "CRYSTALS - Kyber: a CCA-secure module-lattice-based KEM"
# by: Joppe Bos, Leo Ducas, Eike Kiltz, Tancrede Lepoint,
# Vadim Lyubashevsky, John M. Schanck, Peter Schwabe & Damien Stehle
#
# The authors above write in https://eprint.iacr.org/2017/634.pdf:
# > Availability of software. We place all software described in
# > this paper into the public domain to maximize reusability of
# > our results. It is available for download on GitHub: https:
# > //github.com/pq-crystals/kyber.
#
#------------------------------------------------------------------------------

import operator as op
import sys

from math import factorial as fac
from math import log, ceil, erf, sqrt, log2
from math import factorial as fac

def gaussian_center_weight(sigma, t):
    """ Weight of the gaussian of std deviation s, on the interval [-t, t]
    :param x: (float)
    :param y: (float)
    :returns: erf( t / (sigma*\sqrt 2) )
    """
    return erf(t / (sigma * sqrt(2.)))


def binomial(x, y):
    """ Binomial coefficient
    :param x: (integer)
    :param y: (integer)
    :returns: y choose x
    """
    try:
        binom = fac(x) // fac(y) // fac(x - y)
    except ValueError:
        binom = 0
    return binom


def centered_binomial_pdf(k, x):
    """ Probability density function of the centered binomial law of param k at x
    :param k: (integer)
    :param x: (integer)
    :returns: p_k(x)
    """
    return binomial(2*k, x+k) / 2.**(2*k)


def build_centered_binomial_law(k):
    """ Construct the binomial law as a dictionnary
    :param k: (integer)
    :param x: (integer)
    :returns: A dictionnary {x:p_k(x) for x in {-k..k}}
    """
    D = {}
    for i in range(-k, k+1):
        D[i] = centered_binomial_pdf(k, i)
    return D


def mod_switch(x, q, rq):
    """ Modulus switching (rounding to a different discretization of the Torus)
    :param x: value to round (integer)
    :param q: input modulus (integer)
    :param rq: output modulus (integer)
    """
    return int(round(1.* rq * x / q) % rq)


def mod_centered(x, q):
    """ reduction mod q, centered (ie represented in -q/2 .. q/2)
    :param x: value to round (integer)
    :param q: input modulus (integer)
    """
    a = x % q
    if a < q/2:
        return a
    return a - q


def build_mod_switching_error_law(q, rq):
    """ Construct Error law: law of the difference introduced by switching from and back a uniform value mod q
    :param q: original modulus (integer)
    :param rq: intermediate modulus (integer)
    """
    D = {}
    V = {}
    for x in range(q):
        y = mod_switch(x, q, rq)
        z = mod_switch(y, rq, q)
        d = mod_centered(x - z, q)
        D[d] = D.get(d, 0) + 1./q
        V[y] = V.get(y, 0) + 1

    return D


def law_convolution(A, B):
    """ Construct the convolution of two laws (sum of independent variables from two input laws)
    :param A: first input law (dictionnary)
    :param B: second input law (dictionnary)
    """

    C = {}
    for a in A:
        for b in B:
            c = a+b
            C[c] = C.get(c, 0) + A[a] * B[b]
    return C


def law_product(A, B):
    """ Construct the law of the product of independent variables from two input laws
    :param A: first input law (dictionnary)
    :param B: second input law (dictionnary)
    """
    C = {}
    for a in A:
        for b in B:
            c = a*b
            C[c] = C.get(c, 0) + A[a] * B[b]
    return C


def clean_dist(A):
    """ Clean a distribution to accelerate further computation (drop element of the support with proba less than 2^-300)
    :param A: input law (dictionnary)
    """
    B = {}
    for (x, y) in A.items():
        if y>2**(-300):
            B[x] = y
    return B


def iter_law_convolution(A, i):
    """ compute the -ith forld convolution of a distribution (using double-and-add)
    :param A: first input law (dictionnary)
    :param i: (integer)
    """
    D = {0: 1.0}
    i_bin = bin(i)[2:]  # binary representation of n
    for ch in i_bin:
        D = law_convolution(D, D)
        D = clean_dist(D)
        if ch == '1':
            D = law_convolution(D, A)
            D = clean_dist(D)
    return D


def tail_probability(D, t):
    '''
    Probability that an drawn from D is strictly greater than t in absolute value
    :param D: Law (Dictionnary)
    :param t: tail parameter (integer)
    '''
    s = 0
    ma = max(D.keys())
    if t >= ma:
        return 0
    for i in reversed(range(int(ceil(t)), ma)):  # Summing in reverse for better numerical precision (assuming tails are decreasing)
        s += D.get(i, 0) + D.get(-i, 0)
    return s


class KyberParameterSet:

    def __init__(self, n, k, eta1, eta2, q, du, dv):
        self.n = n
        self.k = k
        self.eta1 = eta1     # binary distribution for the secret key
        self.eta2 = eta2    # binary distribution for the ciphertext errors
        self.q = q
        self.du = du  # 2^(bits in the first ciphertext)
        self.dv = dv  # 2^(bits in the second ciphertext)

        self.rqk = ceil(log2(q))  # 2^(bits in the public key)

    def get_final_error_distribution_1d(self):
        """ construct the final error distribution in our encryption scheme
        :param self: parameter set (ParameterSet)
        """
        chis = build_centered_binomial_law(self.eta1)           # LWE error law for the key
        chie = build_centered_binomial_law(self.eta2)        # LWE error law for the ciphertext
        chie_pk = build_centered_binomial_law(self.eta1)
        Rk = build_mod_switching_error_law(self.q, 2**self.rqk)    # Rounding error public key
        Rc = build_mod_switching_error_law(self.q, 2**self.du)    # rounding error first ciphertext
        chiRs = law_convolution(chis, Rk)                   # LWE+Rounding error key
        chiRe = law_convolution(chie, Rc)                   # LWE + rounding error ciphertext

        B1 = law_product(chie_pk, chiRs)                       # (LWE+Rounding error) * LWE (as in a E*S product)
        B2 = law_product(chis, chiRe)

        C1 = iter_law_convolution(B1, self.k * self.n)
        C2 = iter_law_convolution(B2, self.k * self.n)

        C=law_convolution(C1, C2)

        R2 = build_mod_switching_error_law(self.q, 2**self.dv)    # Rounding2 (in the ciphertext mask part)
        F = law_convolution(R2, chie)                       # LWE+Rounding2 error
        D = law_convolution(C, F)                           # Final error
        return D


    def dfr_using_kyber_code(self):
        F = self.get_final_error_distribution_1d()
        proba = tail_probability(F, self.q/4)
        return F, self.n*proba
