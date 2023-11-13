conda environment to create charts: (dumped using `conda env export --from-history`)
```
name: python-charts
channels:
  - conda-forge
dependencies:
  - matplotlib
  - polars
  - numpy
```

To run the tests yourself:
- Install mercat2 and place the jellyfish and kmc binaries in the `progs` folder.
- Download test datasets to the `data` directory, and create list files for KMC (see `helper/kmc.sh`)
- Then, run `run_tests.sh`, `compile_results.py`, and `create_charts.py`.
