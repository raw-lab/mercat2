# -*- coding: utf-8 -*-
"""mercat2_diversity.py: Module for calculating diversity metrics
"""

import ray
from skbio.diversity import alpha as skbio_alpha


def compute_alpha_beta_diversity(basename, counts_tsv, out_file):

    @ray.remote(num_cpus=1)
    def get(func, count):
        try:
            method = getattr(skbio_alpha, func)
            return func, method(count)
        except:
            return func, 'NA'

    counts = list()
    with open(counts_tsv) as reader:
        reader.readline() # skip header
        for line in reader:
            counts += [int(line.split()[1])]

    jobs = list()
    for func in ['shannon', 'simpson', 'simpson_e', 'goods_coverage', 'fisher_alpha', 'dominance', 'chao1', 'chao1_ci', 'ace']:
        jobs += [get.remote(func, counts)]

    results = dict()
    while jobs:
        ready,jobs = ray.wait(jobs)
        if ready:
            key,value = ray.get(ready[0])
            results[key] = value

    with open(out_file, 'w') as writer:
        print('Metric', basename, sep='\t', file=writer)
        for func in ['shannon', 'simpson', 'simpson_e', 'goods_coverage', 'fisher_alpha', 'dominance', 'chao1', 'chao1_ci', 'ace']:
            if func not in results:
                continue
            if type(results[func]) is not str:
                try:
                    value = round(results[func], 2)
                except:
                    value = [round(x,2) for x in results[func]]
            else:
                value = results[func]
            print(func, value, sep='\t', file=writer)
