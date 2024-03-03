# -*- coding: utf-8 -*-
"""mercat2_diversity.py: Module for calculating diversity metrics
"""


from pathlib import Path
import ray
from skbio.diversity import alpha as skbio_alpha
from skbio.diversity import beta_diversity
import matplotlib.pyplot as pyplot


def compute_alpha_diversity(basename, counts_tsv, out_file):

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
    return


def compute_beta_diversity(basename, counts_tsv:Path, outpath:Path):
    outpath.mkdir(0o777, True, True)

    IDs = list()
    counts = list()
    with open(counts_tsv) as reader:
        reader.readline()
        for line in reader:
            line = line.rstrip('\n').split('\t')
            IDs.append(line[0])
            counts.append([int(x) for x in line[1:]])

    for metric in [
            "euclidean",
            "cityblock",
            "braycurtis",
            "canberra",
            "chebyshev",
            "correlation",
            "cosine",
            "dice",
            "hamming",
            "jaccard",
            "mahalanobis", # FAILED: ValueError: The number of observations (5) is too small; the covariance matrix is singular. For observations with 1024 dimensions, at least 1025 observations are required.
            "manhattan",  # aliases to "cityblock" in beta_diversity
            "matching",
            "minkowski",
            "rogerstanimoto",
            "russellrao",
            "seuclidean",
            "sokalmichener",
            "sokalsneath",
            "sqeuclidean",
            "yule",
        ]:
        try:
            distance = beta_diversity(metric, counts, IDs)

            with open(Path(outpath, f"{metric}-{basename}.tsv"), 'w') as writer:
                print('', *IDs, sep='\t', file=writer)
                for i,row in enumerate(distance.data):
                    print(IDs[i], *row, sep='\t', file=writer)
            fig = distance.plot()
            fig.savefig(Path(outpath, f"{metric}-{basename}.png"))
            pyplot.close(fig)
        except Exception as e:
            print(f'Error with beta metric: {metric.capitalize()}')
            print(e.with_traceback(None))

    return
