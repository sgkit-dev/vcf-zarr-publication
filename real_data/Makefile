# Script to download & transform 1000G VCF data for benchmarking
# various compression features of vcf2zarr.
#
# Usage:
#       make CHROMOSOME=22 WGS=1
#       make CHROMOSOME=20 WORKERS=1
#
# Explanation of options:
#       * You can specify the chromosome you'd like to process by
#       setting the CHROMOSOME argument.
#       * You can specify whether to use WGS (Whole Genome Sequencing)
#       or genotyping data by setting the flag WGS. If you use WGS=1,
#       the script will download the WGS data. Otherwise, it will download
#       the genotyping data. Make sure you have enough storage space for the
#       WGS data.
# 	 	* You can specify the number of workers to use for the encoding and exploding
# 	 	by setting the WORKERS argument. The default is 8 workers.


# Required software:
BCFTOOLS=../software/bcftools


# Flags / commandline arguments:
CHROMOSOME ?= 22
WGS ?= 1
WORKERS ?= 8

# If WGS is set to 1, use URL for WGS data:
ifeq ($(WGS), 1)
        TGP_URL = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK/CCDG_13607_B01_GRM_WGS_2019-02-19_chr$(CHROMOSOME).recalibrated_variants.vcf.gz"
        OUTPUT_DIR = data/WGS
else
        # Use URL for genotyping data:
        TGP_URL = "http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr$(CHROMOSOME).phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
        OUTPUT_DIR = data/genotypes
endif

all: $(OUTPUT_DIR)/chr$(CHROMOSOME).zarr

%.vcf.gz:
	@echo "> Downloading VCF file from: $(TGP_URL)"
	mkdir -p $(OUTPUT_DIR) # && wget -O "$(OUTPUT_DIR)/chr$(CHROMOSOME).vcf.gz" "$(TGP_URL)"
	${BCFTOOLS} view "$(TGP_URL)" | head -n 100000 | ${BCFTOOLS} view -O z -o "$(OUTPUT_DIR)/chr$(CHROMOSOME).vcf.gz"
	@echo "---------------------- Downloading completed ----------------------"
	@echo "> Indexing the VCF file..."
	${BCFTOOLS} index "$(OUTPUT_DIR)/chr$(CHROMOSOME).vcf.gz"

.PRECIOUS: %.vcf.gz %.icf %.vcf.gz.tbi

%.icf: %.vcf.gz
	@echo "> Exploding the VCF: $< -> $@"
	vcf2zarr explode $< $@ -p$(WORKERS)
	@echo "---------------------- Exploding completed ----------------------"

%.zarr: %.icf
	@echo "> Encoding the exploded VCF: $< to $@"
	vcf2zarr encode $< $@ -p$(WORKERS) --max-memory=10GiB -v
	@echo "---------------------- Encoding completed ----------------------"

