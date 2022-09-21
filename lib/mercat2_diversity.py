# -*- coding: utf-8 -*-
"""mercat2_diversity.py: Module for calculating diversity metrics
"""

import pandas as pd
from skbio.diversity import alpha as skbio_alpha
from skbio.diversity import beta as skbio_beta


def compute_alpha_beta_diversity(counts_tsv, out_file):

    counts = list()
    with open(counts_tsv) as reader:
        reader.readline() # skip header
        for line in reader:
            counts += [int(line.split()[1])]

    with open(out_file, 'w') as writer:
        #writer.write(abmetric + " = " + str(metrics[abmetric]) + "\n")
        print('shannon', round(skbio_alpha.shannon(counts), 2), sep='\t', file=writer)
        print('simpson', round(skbio_alpha.simpson(counts), 2), sep='\t', file=writer)
        print('simpson_e', round(skbio_alpha.simpson_e(counts), 2), sep='\t', file=writer)
        print('goods_coverage', round(skbio_alpha.goods_coverage(counts), 2), sep='\t', file=writer)
        print('fisher_alpha', round(skbio_alpha.fisher_alpha(counts), 2), sep='\t', file=writer)
        print('dominance', round(skbio_alpha.dominance(counts), 2), sep='\t', file=writer)
        print('chao1', round(skbio_alpha.chao1(counts), 2), sep='\t', file=writer)
        print('chao1_ci', skbio_alpha.chao1_ci(counts), sep='\t', file=writer)
        print('ace', round(skbio_alpha.ace(counts), 2), sep='\t', file=writer)
