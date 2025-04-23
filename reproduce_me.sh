#!/bin/bash

# Copyright 2025 LG Electronics, Inc. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

if [[ -z $1 ]];
then
    echo "[-] Usage $0 <reproduction_directory>"
    exit 1
fi

reproduction_directory=$1

echo "Reproduction directory = $reproduction_directory"

if [ -d $reproduction_directory ];
then
    echo "[-] Directory $reproduction_directory already exists. Quitting..."
    exit 1
fi

echo "[+] Creating directory $reproduction_directory"
mkdir -p $reproduction_directory

set -e # Stop on any error
set -x

# Builds C targets
cmake -B ${reproduction_directory}/build
make -C  ${reproduction_directory}/build

build_directory=${reproduction_directory}/build

# Generates decoding parameters to be used by Kyber2D C implementation
./experiments/scripts/generate_decoding_lines.py

# Runs speed tests for Kyber reference and Kyber with 2D/4D codes
./experiments/scripts/run_speed_tests.sh $reproduction_directory $build_directory

# Generates 2D DFR results
./experiments/scripts/run_exp_kyber2d_dfr.sh $reproduction_directory $build_directory

# Generates 4D DFR results
./experiments/scripts/run_exp_kyber4d_dfr.sh $reproduction_directory $build_directory

set +x
