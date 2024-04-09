BCFTOOLS=../software/bcftools
BCFTOOLS_PLUGINS=../software/bcftools-1.18/plugins
SAVVY=../software/savvy/bin/sav
GENOZIP=../software/genozip

simulated_chrom_21.ts.tsz:
	wget https://zenodo.org/record/7702392/files/simulated_chrom_21.ts.tsz


%.ts: %.ts.tsz
	tsunzip -k $<

all: data/chr21_10_1.zarr data/chr21_10_1.sav data/chr21_10_1.genozip data/chr21_10_1.tags.genozip \
	data/chr21_10_2.zarr data/chr21_10_2.sav data/chr21_10_2.genozip data/chr21_10_2.tags.genozip \
	data/chr21_10_3.zarr data/chr21_10_3.sav data/chr21_10_3.genozip data/chr21_10_3.tags.genozip \
	data/chr21_10_4.zarr data/chr21_10_4.sav data/chr21_10_4.genozip data/chr21_10_4.tags.genozip \
	data/chr21_10_5.zarr data/chr21_10_5.sav data/chr21_10_5.genozip data/chr21_10_5.tags.genozip \
	data/chr21_10_6.zarr data/chr21_10_6.sav data/chr21_10_6.genozip data/chr21_10_6.tags.genozip 

.PRECIOUS: %.tags.bcf %.bcf %.vcf.gz %.icf %.ts

data/chr21_10_1.ts: simulated_chrom_21.ts
	python3 subset_trees.py

%.vcf.gz: %.ts
	python3 -m tskit vcf $< | bgzip > $@
	${BCFTOOLS} index $@

%.bcf: %.vcf.gz
	${BCFTOOLS} view -O b $< > $@
	${BCFTOOLS} index $@

%.icf: %.bcf
	vcf2zarr explode $< $@ -p8

%.zarr: %.icf
	vcf2zarr encode $< $@ -p8 --max-memory=200GiB -v

# Compute some minimal INFO summaries to support the af-dist computation
%.tags.bcf: %.bcf
	BCFTOOLS_PLUGINS=${BCFTOOLS_PLUGINS} ${BCFTOOLS} \
			 +fill-tags $< -Ob -o $@ -- -t AN,AC,AF
	${BCFTOOLS} index $@

%.tags.vcf.gz: %.tags.bcf
	${BCFTOOLS} view -O z $< > $@
	${BCFTOOLS} index $@

%.sav: %.bcf
	${SAVVY} import --phasing full $< > $@

%.genozip: %.bcf
	# There's some weird bug in genocat that doesn't allow us 
	# to run --samples on output if it's derived directly from
	# bcf, so we have to filter it through bcftools first.
	${BCFTOOLS} view $< | ${GENOZIP} --no-test -o $@
