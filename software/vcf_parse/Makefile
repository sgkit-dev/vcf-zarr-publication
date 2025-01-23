
vcf_parse: vcf_parse.c
	gcc -I../bcftools-1.18/htslib-1.18/htslib \
		-L../bcftools-1.18/htslib-1.18 \
		-o vcf_parse vcf_parse.c -lhts -lpthread -lz -lm -lbz2 -llzma -lcurl -lcrypto
