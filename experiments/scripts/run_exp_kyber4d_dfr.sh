#!/usr/bin/bash

# Copyright 2025 LG Electronics, Inc. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# IMPORTANT NOTICE: Scripts should be run from the root directory

# Expected time: About 24 hours

if [[ $# != 2 ]];
then
    echo "[-] Usage $0 <reproduction_directory> <build_directory>"
    exit 1
fi

reproduction_directory=$1
build_directory=$2

echo "[+] Reproduction directory = $reproduction_directory"
echo "[+] Build directory = $build_directory"


set -e # Abort on error

set -x

mkdir -p $reproduction_directory/results/kyber4d_dfr

$build_directory/minal_src/dfr_computation/minal4d_dfr_computation  \
    experiments/setup/exp_kyber4d_dfr.csv | tee \
    $reproduction_directory/results/kyber4d_dfr/result_kyber4d_dfr.csv

$build_directory/minal_src/dfr_computation/minal4d_dfr_computation_with_pk_compression11  \
    experiments/setup/exp_kyber4d_dfr_pk_compression_11bits.csv | tee \
    $reproduction_directory/results/kyber4d_dfr/result_kyber4d_dfr_pk_compression_11bits.csv

set +x
