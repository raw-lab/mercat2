# -*- coding: utf-8 -*-
"""mercat2_diversity.py: Module for calculating diversity metrics
"""

import ray
from skbio.diversity import alpha as skbio_alpha
from skbio.diversity import beta as skbio_beta


def compute_alpha_beta_diversity(counts_tsv, out_file):

    @ray.remote(num_cpus=1)
    def get(func, count):
        try:
            method = getattr(skbio_alpha, func)
            return func, method(count)
        except:
            return func, 'NONE'

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
            if value is not 'NONE':
                results[key] = value

    with open(out_file, 'w') as writer:
        #writer.write(abmetric + " = " + str(metrics[abmetric]) + "\n")
        print('shannon', round(results['shannon'], 2), sep='\t', file=writer)
        print('simpson', round(results['simpson'], 2), sep='\t', file=writer)
        print('simpson_e', round(results['simpson_e'], 2), sep='\t', file=writer)
        print('goods_coverage', round(results['goods_coverage'], 2), sep='\t', file=writer)
        print('fisher_alpha', round(results['fisher_alpha'], 2), sep='\t', file=writer)
        print('dominance', round(results['dominance'], 2), sep='\t', file=writer)
        print('chao1', round(results['chao1'], 2), sep='\t', file=writer)
        print('chao1_ci', results['chao1_ci'], sep='\t', file=writer)
        print('ace', round(results['ace'], 2), sep='\t', file=writer)
