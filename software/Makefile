all: bcftools savvy-afdist/savvy-afdist genozip

bcftools:
	wget https://github.com/samtools/bcftools/releases/download/1.18/bcftools-1.18.tar.bz2
	tar -jxvf bcftools-1.18.tar.bz2
	cd bcftools-1.18 && ./configure && make -j 8
	cp bcftools-1.18/bcftools ./

genozip:
	# NOTE: using this specific version of genozip because it started using some 
	# assembly instructions not supported by the CPU I'm using at some point after
	# in the months that.
	wget https://github.com/divonlan/genozip/archive/refs/tags/genozip-15.0.26.tar.gz
	tar -xf genozip-15.0.26.tar.gz
	tar -xf genozip-genozip-15.0.26/installers/genozip-linux-x86_64.tar
	cp genozip-linux-x86_64/* ./

savvy:
	# Assumes cget is installed:
	# python -m pip install cget;
	cget install --prefix savvy statgen/savvy

savvy-afdist/savvy-afdist: savvy
	cmake -DCMAKE_TOOLCHAIN_FILE=../savvy/cget/cget.cmake \
		-DCMAKE_BUILD_TYPE=Release -B savvy-afdist savvy-afdist
	make -C savvy-afdist

clean:
	rm -fR bcftools* savvy* geno*


    
