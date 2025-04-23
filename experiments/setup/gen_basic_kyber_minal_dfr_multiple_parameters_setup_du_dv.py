#!/usr/bin/python3
# Copyright 2025 LG Electronics, Inc. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0


from dataclasses import dataclass

@dataclass
class KyberParams:
    q: int
    level: int
    k: int
    eta1: int
    eta2: int
    du: int
    dv: int

    def __str__(self):
        printf()

KYBER_BASE_PARAMS = {
    1: KyberParams(q=3329,level=1,k=2,eta1=3,eta2=2,du=10,dv=4),
    3: KyberParams(q=3329,level=3,k=3,eta1=2,eta2=2,du=10,dv=4),
    5: KyberParams(q=3329,level=5,k=4,eta1=2,eta2=2,du=11,dv=5),
}

def get_ciphertext_size(k, du, dv):
    return (256 // 8) * (k * du + dv)

def get_ciphertext_size_ext(k, n_blocks_for_du0, du0, du1, dv):
    assert(n_blocks_for_du0 <= k)
    return (256 // 8) * (n_blocks_for_du0 * du0 + (k - n_blocks_for_du0) * du1 + dv)

# Only print parameters if they do not increase the ciphertext size
def print_if_compressed(level, k, n_blocks_for_du0, du0, du1, dv):
    base = KYBER_BASE_PARAMS[level]
    if get_ciphertext_size_ext(k, n_blocks_for_du0, du0, du1, dv) <= get_ciphertext_size(base.k, base.du, base.dv):
        print(f'{base.q},{level},{k},{base.eta1},{base.eta2},{n_blocks_for_du0},{du0},{du1},{dv}')

def print_header():
    print('KYBER_Q,KYBER_LEVEL,KYBER_K,KYBER_ETA1,KYBER_ETA2,KYBER_N_BLOCKS_FOR_DU0,KYBER_DU0,KYBER_DU1,KYBER_DV')

def main():
    print_header()

    for l in KYBER_BASE_PARAMS:
        base_params = KYBER_BASE_PARAMS[l]

        n_blocks_for_du0 = base_params.k
        du1 = 0
        for du0 in range(base_params.du - 1, base_params.du + 2):
            for dv in range(base_params.dv - 1, base_params.dv + 2):
                print_if_compressed(base_params.level, base_params.k, n_blocks_for_du0, du0, du1, dv)


if __name__ == '__main__':
    main()