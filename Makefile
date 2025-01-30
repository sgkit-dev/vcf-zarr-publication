FIGURES=figures/data-scaling.pdf\
	figures/whole-matrix-compute.pdf\
	figures/whole-matrix-compute-zarr-versions.pdf\
	figures/whole-matrix-decode.pdf\
	figures/subset-matrix-compute.pdf\
	figures/subset-matrix-compute-supplemental.pdf\
	figures/s3-throughput.pdf\
	figures/s3-network-throughput.pdf \
	figures/ofh-histograms.pdf

all: paper.pdf response-to-reviewers.pdf

paper.aux: paper.tex
	pdflatex -shell-escape paper.tex

paper.bbl: paper.aux paper.bib
	bibtex paper
	pdflatex -shell-escape paper.tex

paper.pdf: $(FIGURES) paper.bbl
	pdflatex -shell-escape paper.tex

paper.ps: paper.dvi
	dvips paper

paper.dvi: paper.tex paper.bib
	latex paper.tex
	bibtex paper
	latex paper.tex
	latex paper.tex

.PHONY: spellcheck
spellcheck: aspell.conf
	aspell --conf ./aspell.conf --check paper.tex

clean:
	rm -f *.log *.dvi *.aux
	rm -f *.blg *.bbl
	rm -fR _minted*

mrproper: clean
	rm -f *.ps *.pdf



TS_FILES=scaling/data/chr21_10_1.ts\
	scaling/data/chr21_10_2.ts\
	scaling/data/chr21_10_3.ts\
	scaling/data/chr21_10_4.ts\
	scaling/data/chr21_10_5.ts\
	scaling/data/chr21_10_6.ts


# Should probably add rules for other data collection. However, it
# takes a *long* time to do some stuff, and needs to be done in bits
# in practise.
plot_data/data-scaling.csv:
	python3 src/collect_data.py file-size $(TS_FILES) -o $@


# Compression benchmarks:

plot_data/shuffle_benchmarks.csv:
	python3 src/compression_benchmarks.py --input real_data/data/WGS/chr22.zarr \
		--test-config shuffle \
		-o $@

plot_data/compressor_benchmarks.csv:
	python3 src/compression_benchmarks.py --input real_data/data/WGS/chr22.zarr \
		--test-config compressor \
		-o $@

plot_data/chunksize_benchmarks.csv:
	python3 src/compression_benchmarks.py --input real_data/data/WGS/chr22.zarr \
		--test-config chunksize \
		-o $@

plot_data/chunksize_finegrained_benchmarks.csv:
	python3 src/compression_benchmarks.py --input real_data/data/WGS/chr22.zarr \
		--test-config chunksize_finegrained \
		-o $@


# TODO make rule for time-scaling

# TODO make some substitution rules for this later
figures/data-scaling.pdf: plot_data/data-scaling.csv
	python3 src/plot.py data-scaling plot_data/data-scaling.csv  \
		figures/data-scaling.pdf

figures/s3-throughput.pdf: plot_data/gel-s3-throughput.csv
	python3 src/plot.py s3-throughput plot_data/gel-s3-throughput.csv  \
		figures/s3-throughput.pdf

figures/s3-network-throughput.pdf: plot_data/gel-s3-network-throughput.csv
	python3 src/plot.py s3-network-throughput plot_data/gel-s3-network-throughput.csv  \
		figures/s3-network-throughput.pdf

figures/ofh-histograms.pdf:
	python3 src/plot.py ofh-histograms figures/ofh-histograms.pdf

figures/whole-matrix-compute.pdf: plot_data/whole-matrix-compute.csv
	python3 src/plot.py whole-matrix-compute plot_data/whole-matrix-compute.csv  \
		figures/whole-matrix-compute.pdf

figures/whole-matrix-compute-zarr-versions.pdf: plot_data/whole-matrix-compute-zarr-versions.csv
	python3 src/plot.py whole-matrix-compute-zarr-versions plot_data/whole-matrix-compute-zarr-versions.csv  \
		figures/whole-matrix-compute-zarr-versions.pdf

figures/whole-matrix-decode.pdf: plot_data/whole-matrix-decode.csv
	python3 src/plot.py whole-matrix-decode plot_data/whole-matrix-decode.csv  \
		figures/whole-matrix-decode.pdf

figures/column-extract.pdf: plot_data/column-extract.csv
	python3 src/plot.py column-extract plot_data/column-extract.csv  \
		figures/column-extract.pdf

figures/subset-matrix-compute.pdf: plot_data/subset-matrix-compute.csv
	python3 src/plot.py subset-matrix-compute plot_data/subset-matrix-compute.csv  \
		figures/subset-matrix-compute.pdf

figures/subset-matrix-compute-supplemental.pdf: plot_data/subset-matrix-compute.csv
	python3 src/plot.py subset-matrix-compute-supplemental \
		plot_data/subset-matrix-compute.csv  \
		figures/subset-matrix-compute-supplemental.pdf

figures/compression-shuffle.pdf: plot_data/shuffle_benchmarks.csv
	python3 src/plot.py compression-shuffle \
                plot_data/shuffle_benchmarks.csv \
                $@

figures/compression-compressor.pdf: plot_data/compressor_benchmarks.csv
	python3 src/plot.py compression-compressor \
                plot_data/compressor_benchmarks.csv \
                $@

figures/compression-chunksize.pdf: plot_data/chunksize_benchmarks.csv
	python3 src/plot.py compression-chunksize \
		plot_data/chunksize_benchmarks.csv \
		$@

figures/compression-chunksize-finegrained.pdf: plot_data/chunksize_finegrained_benchmarks.csv
	python3 src/plot.py compression-chunksize-finegrained \
		plot_data/chunksize_finegrained_benchmarks.csv \
		$@


review-diff.tex: paper.tex
	latexdiff reviewed-paper.tex paper.tex > review-diff.tex

review-diff.pdf: review-diff.tex
	pdflatex review-diff.tex
	pdflatex review-diff.tex
	bibtex review-diff
	pdflatex review-diff.tex


response-to-reviewers.pdf: response-to-reviewers.tex
	pdflatex $<
