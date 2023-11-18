conda environment to create charts: (dumped using `conda env export --from-history`)
```
name: mercat2-charts
channels:
  - conda-forge
dependencies:
  - matplotlib
  - polars
  - numpy
```

Setup:
- Install mercat2 and place the jellyfish and kmc binaries in the `progs` folder.
- Download test datasets from OSF (https://osf.io/4bhp2/) to the `data` directory

Run:
- `run_tests.sh`
- `compile_results.py`
- Activate your conda environment for creating the graphs, then `create_charts.py -i results.csv`.
