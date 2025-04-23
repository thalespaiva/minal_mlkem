#!/usr/bin/bash

# Copyright 2025 LG Electronics, Inc. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# IMPORTANT NOTICE: Scripts should be run from the root directory

# Expected time: 5 minutes

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

mkdir -p $reproduction_directory/results/performance/ref

for ref_speed_test in $build_directory/minal_src/kem_implementation/ref/test_speed*; do
    ./$ref_speed_test | tee $reproduction_directory/results/performance/ref/$(basename ${ref_speed_test}.txt)
done

mkdir -p $reproduction_directory/results/performance/avx2

for avx2_speed_test in $build_directory/minal_src/kem_implementation/avx2/test_speed*; do
    ./$avx2_speed_test | tee $reproduction_directory/results/performance/avx2/$(basename ${avx2_speed_test}.txt)
done

./experiments/scripts/raw_performance_to_csv.py $reproduction_directory/results/performance/avx2 | \
    tee $reproduction_directory/results/performance/avx2_performance.csv

./experiments/scripts/raw_performance_to_csv.py $reproduction_directory/results/performance/ref | \
    tee $reproduction_directory/results/performance/ref_performance.csv

set +x
