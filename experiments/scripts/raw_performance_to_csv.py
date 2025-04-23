#!/usr/bin/python3
# Copyright 2025 LG Electronics, Inc. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0


import sys
import os

def extract_data_from_file_name(f):
    level = None        # int
    code_dim = None     # int
    null_eta2 = None    # bool

    if '512_' in f:
        level = 1
    elif '768_' in f:
        level = 3
    elif '1024_' in f:
        level = 5
    else:
        raise ValueError(f'Could not find the security level of file {f}')

    if '_1d' in f:
        code_dim = 1
    elif '_2d' in f:
        code_dim = 2
    elif '_4d' in f:
        code_dim = 4
    else:
        code_dim = 1

    if 'null_eta2' in f:
        null_eta2 = True
    else:
        null_eta2 = False

    assert level is not None
    assert code_dim is not None
    assert null_eta2 is not None

    return {
        'level': level,
        'code_dim': code_dim,
        'null_eta2': null_eta2,
    }

def extract_measurements_from_file(filepath):

    f = open(filepath, 'r')

    measurements = {
        'function': [],
        'cycles_avg': [],
        'cycles_median': [],
    }

    for l in f:
        if l.startswith('median'):
            measurements['cycles_median'].append(int(l.split()[1]))
        elif l.startswith('average'):
            measurements['cycles_avg'].append(int(l.split()[1]))
        elif ':' in l:
            # Removes trailing ':' in l:
            measurements['function'].append(l.strip()[:-1])
    f.close()
    assert(len(measurements['function']) == len(measurements['cycles_avg']))
    assert(len(measurements['function']) == len(measurements['cycles_median']))

    return measurements

def main(raw_performance_directory):

    print('level,function,cycles_avg,cycles_median,code_dimension,null_eta2')
    for f in os.listdir(raw_performance_directory):
        if f.startswith('test_speed'):
            base_data = extract_data_from_file_name(f)
            measurements = extract_measurements_from_file(os.path.join(raw_performance_directory, f))

            for function, cycles_avg, cycles_median in (zip(measurements['function'],
                                                            measurements['cycles_avg'],
                                                            measurements['cycles_median'])):
                print(base_data['level'], end=',')
                print(function, end=',')
                print(cycles_avg, end=',')
                print(cycles_median, end=',')
                print(base_data['code_dim'], end=',')
                print(base_data['null_eta2'])

if __name__ == '__main__':

    if len(sys.argv) != 2:
        print(f'Usage: {sys.argv[0]} raw_performance_directory')
        sys.exit(1)

    main(sys.argv[1])