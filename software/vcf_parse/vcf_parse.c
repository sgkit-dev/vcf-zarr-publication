#include <assert.h>
#include <err.h>
#include <hts.h>
#include <stdio.h>
#include <stdlib.h>
#include <vcf.h>

static void read_vcf(char *path) {
    htsFile *h = hts_open(path, "r");

    if (h == NULL) {
        errx(EXIT_FAILURE, "Error opening VCF");
    }
    bcf_hdr_t *hdr = bcf_hdr_read(h);
    if (hdr == NULL) {
        errx(EXIT_FAILURE, "Error reading header");
    }
    bcf1_t *b = bcf_init();
    int ret = 0;
    int records = 0;
    while ((ret = bcf_read(h, hdr, b)) == 0) {
        /* printf("%d\n", b->pos); */
        records++;
    }
    printf("read = %d records\n", records);

    bcf_destroy(b);
    bcf_hdr_destroy(hdr);
    hts_close(h);
}

int main(int argc, char **argv) {
    if (argc != 2) {
        errx(EXIT_FAILURE, "usage: vcf_parse <filename>\n");
    }
    printf("Reading VCF %s\n", argv[1]);
    read_vcf(argv[1]);
    return 0;
}
