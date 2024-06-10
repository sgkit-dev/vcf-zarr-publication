# vcf-zarr-publication

This repo contains the manuscript for the publication describing the vcf-zarr
specification and its compression and query performance on several datasets.
All code required to generate figures and example analyses is in this
repo.

## Layout

- The main text is in the paper.tex/paper.bib files.

- Building the document and creating plots is automated using the
Makefile. Building datasets and running benchmarks are also semi-automated
using the main Makefile (but not entirely, as these benchmarks take a lot
of time to run, and needed to be done bit-by-bit). Run python3
``src/collect_data.py --help`` to see the available commands.

- Code for generating all plots is in ``src/plot.py``, and the data for these
plots is stored in CVS form in ``plot_data``.

- Code for running the data compression analysis on 1000 Genomes data
is in ``src/compression_benchmarks.py``. The dataset is downloaded
and managed in the ``real_data`` directory. See the Makefile there
for details on downloading and converting to the inital Zarr.

- The ``gel_example`` directory contains the Jupyter notebooks used to run the
benchmarks for the Genomics England example.

- Code for running simulation-based scaling benchmarks ``src/collect_data.py``,
  and the actual Zarr-Python benchmarks used in ``src/zarr_afdist.py``. The
``software`` directory contains the ``savvy-afdist`` C++ code used in the
benchmarks, along with downloaded versions of bcftools etc. The Makefile in
this directory should take care of downloading and compiling all of these. The
Makefile is the starting point here, which should take care of downloading
dataset, doing the subsetting and running the various conversion tools. The
simulated dataset and various copies are stored in ``scaling``, along with some
utilities.

To run the simulation based benchmarks:

1. cd to the ``software`` directory and ``make`` (you may need to install some dependencies).
2. cd to the ``scaling`` directory and ``make``. This will take a *long* time and need a
lot of storage space.
3. Run the various benchmarks using ``python src/collect_data.py``

