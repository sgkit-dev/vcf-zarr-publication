### Info

Chromosome 22 genotype files taken from [release 9](https://ourfuturehealth.gitbook.io/our-future-health/data/data-releases/2024-data-releases/release-9).

- chr 22
- split across three block files
- 10,221 variants
- 651,050 samples

Format fileds of interest
- LRR: Log R Ratio (7 dp float)
- BAF: B allele frequency (7 dp float)
- GS: GenCall Score (4 dp float)

Machine (Azure VM)
- Standard_D32s_v5
- 32 CPUs
- 128 GB RAM

### Transform genotype VCF to column format 

```
INPUT_NAME=xxx

time vcf2zarr explode --force -p 20 -Q ofh.chr22-b0.vcf.gz ofh.chr22-b1.vcf.gz ofh.chr22-b2.vcf.gz $INPUT_NAME.icf

real    7m53.434s
user    137m11.811s
sys     2m39.746s
```

### Make schema

```
vcf2zarr mkschema $INPUT_NAME.icf > $INPUT_NAME.schema.json
```

### Modify schemas

Create two new schemas with filters applied to largest format fields (`call_LRR`, `call_BAF`, `call_GS`)

- quantize filter

```
{
    "id": "quantize",
    "digits": 5,
    "dtype": "<f4",
    "astype": "<f4",
}
```

- bitround filter

```
{
    "id": "bitround", 
    "keepbits": 5
}
```

### Convert to Zarr

- Create Zarr objects for default and using new filters
- Run `vcf2zarr inspect`

```
time vcf2zarr encode --force -p 20 -Q $INPUT_NAME.icf $INPUT_NAME.vcz

real    6m33.836s
user    57m51.796s
sys     2m49.234s

name                   dtype    stored      size            ratio    nchunks  chunk_size               avg_chunk_stored    shape               chunk_shape       compressor                                                      filters
---------------------  -------  ----------  ---------  ----------  ---------  -----------------------  ------------------  ------------------  ----------------  --------------------------------------------------------------  ------------
/call_LRR              float32  22.1 GiB    24.79 GiB      1.1           726  34.96 MiB                31.17 MiB           (10221, 651050)     (1000, 10000)     Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   None
/call_BAF              float32  15.12 GiB   24.79 GiB      1.6           726  34.96 MiB                21.33 MiB           (10221, 651050)     (1000, 10000)     Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   None
/call_GS               float32  2.09 GiB    24.79 GiB     12             726  34.96 MiB                2.95 MiB            (10221, 651050)     (1000, 10000)     Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   None
/call_genotype         int8     657.21 MiB  12.39 GiB     19             726  17.48 MiB                926.97 KiB          (10221, 651050, 2)  (1000, 10000, 2)  Blosc(cname='zstd', clevel=7, shuffle=BITSHUFFLE, blocksize=0)  None
/call_genotype_mask    bool     79.48 MiB   12.39 GiB    160             726  17.48 MiB                112.1 KiB           (10221, 651050, 2)  (1000, 10000, 2)  Blosc(cname='zstd', clevel=7, shuffle=BITSHUFFLE, blocksize=0)  None
/sample_id             object   2.93 MiB    4.97 MiB       1.7            66  77.07 KiB                45.44 KiB           (651050,)           (10000,)          Blosc(cname='zstd', clevel=7, shuffle=SHUFFLE, blocksize=0)     [VLenUTF8()]
/call_genotype_phased  bool     644.72 KiB  6.2 GiB    10000             726  8.74 MiB                 909 bytes           (10221, 651050)     (1000, 10000)     Blosc(cname='zstd', clevel=7, shuffle=BITSHUFFLE, blocksize=0)  None
/variant_allele        object   57.45 KiB   159.7 KiB      2.8            11  14.52 KiB                5.22 KiB            (10221, 2)          (1000, 2)         Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   [VLenUTF8()]
/variant_filter        bool     48.97 KiB   9.98 KiB       0.2            11  929.1818181818181 bytes  4.45 KiB            (10221, 1)          (1000, 1)         Blosc(cname='zstd', clevel=7, shuffle=BITSHUFFLE, blocksize=0)  None
/variant_id            object   48.11 KiB   79.85 KiB      1.7            11  7.26 KiB                 4.37 KiB            (10221,)            (1000,)           Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   [VLenUTF8()]
/variant_position      int32    36.37 KiB   39.93 KiB      1.1            11  3.63 KiB                 3.31 KiB            (10221,)            (1000,)           Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   None
/region_index          int32    8.62 KiB    264 bytes      0.03            1  264.0 bytes              8.62 KiB            (11, 6)             (11, 6)           Blosc(cname='zstd', clevel=9, shuffle=NOSHUFFLE, blocksize=0)   None
/variant_length        int8     5.17 KiB    9.98 KiB       1.9            11  929.1818181818181 bytes  481 bytes           (10221,)            (1000,)           Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   None
/variant_contig        int8     4.99 KiB    9.98 KiB       2              11  929.1818181818181 bytes  464 bytes           (10221,)            (1000,)           Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   None
/variant_quality       float32  4.93 KiB    39.93 KiB      8.1            11  3.63 KiB                 459 bytes           (10221,)            (1000,)           Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   None
/variant_id_mask       bool     4.87 KiB    9.98 KiB       2              11  929.1818181818181 bytes  453 bytes           (10221,)            (1000,)           Blosc(cname='zstd', clevel=7, shuffle=BITSHUFFLE, blocksize=0)  None
/contig_id             object   4.51 KiB    200 bytes      0.043           1  200.0 bytes              4.51 KiB            (25,)               (25,)             Blosc(cname='zstd', clevel=7, shuffle=SHUFFLE, blocksize=0)     [VLenUTF8()]
/contig_length         int64    4.5 KiB     200 bytes      0.043           1  200.0 bytes              4.5 KiB             (25,)               (25,)             Blosc(cname='zstd', clevel=7, shuffle=SHUFFLE, blocksize=0)     None
/filter_id             object   4.43 KiB    8 bytes        0.0018          1  8.0 bytes                4.43 KiB            (1,)                (1,)              Blosc(cname='zstd', clevel=7, shuffle=SHUFFLE, blocksize=0)     [VLenUTF8()]
```

```
time vcf2zarr encode --force -p 20 -Qs schema_bitround.json $INPUT_NAME.icf $INPUT_NAME.bitround.vcz

real    11m21.202s
user    105m17.498s
sys     2m20.037s

name                   dtype    stored      size            ratio    nchunks  chunk_size               avg_chunk_stored    shape               chunk_shape       compressor                                                      filters
---------------------  -------  ----------  ---------  ----------  ---------  -----------------------  ------------------  ------------------  ----------------  --------------------------------------------------------------  ----------------------
/call_LRR              float32  8.35 GiB    24.79 GiB      3             726  34.96 MiB                11.78 MiB           (10221, 651050)     (1000, 10000)     Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   [BitRound(keepbits=5)]
/call_BAF              float32  6.02 GiB    24.79 GiB      4.1           726  34.96 MiB                8.49 MiB            (10221, 651050)     (1000, 10000)     Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   [BitRound(keepbits=5)]
/call_GS               float32  961.37 MiB  24.79 GiB     26             726  34.96 MiB                1.32 MiB            (10221, 651050)     (1000, 10000)     Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   [BitRound(keepbits=5)]
/call_genotype         int8     657.21 MiB  12.39 GiB     19             726  17.48 MiB                926.97 KiB          (10221, 651050, 2)  (1000, 10000, 2)  Blosc(cname='zstd', clevel=7, shuffle=BITSHUFFLE, blocksize=0)  None
/call_genotype_mask    bool     79.48 MiB   12.39 GiB    160             726  17.48 MiB                112.1 KiB           (10221, 651050, 2)  (1000, 10000, 2)  Blosc(cname='zstd', clevel=7, shuffle=BITSHUFFLE, blocksize=0)  None
/sample_id             object   2.93 MiB    4.97 MiB       1.7            66  77.07 KiB                45.44 KiB           (651050,)           (10000,)          Blosc(cname='zstd',
clevel=7, shuffle=SHUFFLE, blocksize=0)     [VLenUTF8()]
/call_genotype_phased  bool     644.72 KiB  6.2 GiB    10000             726  8.74 MiB                 909 bytes           (10221, 651050)     (1000, 10000)     Blosc(cname='zstd',
clevel=7, shuffle=BITSHUFFLE, blocksize=0)  None
/variant_allele        object   57.45 KiB   159.7 KiB      2.8            11  14.52 KiB                5.22 KiB            (10221, 2)          (1000, 2)         Blosc(cname='zstd',
clevel=7, shuffle=NOSHUFFLE, blocksize=0)   [VLenUTF8()]
/variant_filter        bool     48.97 KiB   9.98 KiB       0.2            11  929.1818181818181 bytes  4.45 KiB            (10221, 1)          (1000, 1)         Blosc(cname='zstd',
clevel=7, shuffle=BITSHUFFLE, blocksize=0)  None
/variant_id            object   48.11 KiB   79.85 KiB      1.7            11  7.26 KiB                 4.37 KiB            (10221,)            (1000,)           Blosc(cname='zstd',
clevel=7, shuffle=NOSHUFFLE, blocksize=0)   [VLenUTF8()]
/variant_position      int32    36.37 KiB   39.93 KiB      1.1            11  3.63 KiB                 3.31 KiB            (10221,)            (1000,)           Blosc(cname='zstd',
clevel=7, shuffle=NOSHUFFLE, blocksize=0)   None
/region_index          int32    8.62 KiB    264 bytes      0.03            1  264.0 bytes              8.62 KiB            (11, 6)             (11, 6)           Blosc(cname='zstd',
clevel=9, shuffle=NOSHUFFLE, blocksize=0)   None
/variant_length        int8     5.17 KiB    9.98 KiB       1.9            11  929.1818181818181 bytes  481 bytes           (10221,)            (1000,)           Blosc(cname='zstd',
clevel=7, shuffle=NOSHUFFLE, blocksize=0)   None
/variant_contig        int8     4.99 KiB    9.98 KiB       2              11  929.1818181818181 bytes  464 bytes           (10221,)            (1000,)           Blosc(cname='zstd',
clevel=7, shuffle=NOSHUFFLE, blocksize=0)   None
/variant_quality       float32  4.93 KiB    39.93 KiB      8.1            11  3.63 KiB                 459 bytes           (10221,)            (1000,)           Blosc(cname='zstd',
clevel=7, shuffle=NOSHUFFLE, blocksize=0)   None
/variant_id_mask       bool     4.87 KiB    9.98 KiB       2              11  929.1818181818181 bytes  453 bytes           (10221,)            (1000,)           Blosc(cname='zstd',
clevel=7, shuffle=BITSHUFFLE, blocksize=0)  None
/contig_id             object   4.51 KiB    200 bytes      0.043           1  200.0 bytes              4.51 KiB            (25,)               (25,)             Blosc(cname='zstd',
clevel=7, shuffle=SHUFFLE, blocksize=0)     [VLenUTF8()]
/contig_length         int64    4.5 KiB     200 bytes      0.043           1  200.0 bytes              4.5 KiB             (25,)               (25,)             Blosc(cname='zstd',
clevel=7, shuffle=SHUFFLE, blocksize=0)     None
/filter_id             object   4.43 KiB    8 bytes        0.0018          1  8.0 bytes                4.43 KiB            (1,)                (1,)              Blosc(cname='zstd',
clevel=7, shuffle=SHUFFLE, blocksize=0)     [VLenUTF8()]
```

```
time vcf2zarr encode --force -p 20 -Qs schema_quantize.json $INPUT_NAME.icf $INPUT_NAME.quantize.vcz

real    8m17.282s
user    75m46.997s
sys     2m30.519s

name                   dtype    stored      size            ratio    nchunks  chunk_size               avg_chunk_stored    shape               chunk_shape       compressor
                                            filters
---------------------  -------  ----------  ---------  ----------  ---------  -----------------------  ------------------  ------------------  ----------------  --------------------
------------------------------------------  ---------------------------------
/call_LRR              float32  16.31 GiB   24.79 GiB      1.5           726  34.96 MiB                23.01 MiB           (10221, 651050)     (1000, 10000)     Blosc(cname='zstd',
clevel=7, shuffle=NOSHUFFLE, blocksize=0)   [Quantize(digits=5, dtype='<f4')]
/call_BAF              float32  11.34 GiB   24.79 GiB      2.2           726  34.96 MiB                15.99 MiB           (10221, 651050)     (1000, 10000)     Blosc(cname='zstd',
clevel=7, shuffle=NOSHUFFLE, blocksize=0)   [Quantize(digits=5, dtype='<f4')]
/call_GS               float32  1.81 GiB    24.79 GiB     14             726  34.96 MiB                2.56 MiB            (10221, 651050)     (1000, 10000)     Blosc(cname='zstd',
clevel=7, shuffle=NOSHUFFLE, blocksize=0)   [Quantize(digits=5, dtype='<f4')]
/call_genotype         int8     657.21 MiB  12.39 GiB     19             726  17.48 MiB                926.97 KiB          (10221, 651050, 2)  (1000, 10000, 2)  Blosc(cname='zstd',
clevel=7, shuffle=BITSHUFFLE, blocksize=0)  None
/call_genotype_mask    bool     79.48 MiB   12.39 GiB    160             726  17.48 MiB                112.1 KiB           (10221, 651050, 2)  (1000, 10000, 2)  Blosc(cname='zstd',
clevel=7, shuffle=BITSHUFFLE, blocksize=0)  None
/sample_id             object   2.93 MiB    4.97 MiB       1.7            66  77.07 KiB                45.44 KiB           (651050,)           (10000,)          Blosc(cname='zstd',
clevel=7, shuffle=SHUFFLE, blocksize=0)     [VLenUTF8()]
/call_genotype_phased  bool     644.72 KiB  6.2 GiB    10000             726  8.74 MiB                 909 bytes           (10221, 651050)     (1000, 10000)     Blosc(cname='zstd',
clevel=7, shuffle=BITSHUFFLE, blocksize=0)  None
/variant_allele        object   57.45 KiB   159.7 KiB      2.8            11  14.52 KiB                5.22 KiB            (10221, 2)          (1000, 2)         Blosc(cname='zstd',
clevel=7, shuffle=NOSHUFFLE, blocksize=0)   [VLenUTF8()]
/variant_filter        bool     48.97 KiB   9.98 KiB       0.2            11  929.1818181818181 bytes  4.45 KiB            (10221, 1)          (1000, 1)         Blosc(cname='zstd',
clevel=7, shuffle=BITSHUFFLE, blocksize=0)  None
/variant_id            object   48.11 KiB   79.85 KiB      1.7            11  7.26 KiB                 4.37 KiB            (10221,)            (1000,)           Blosc(cname='zstd',
clevel=7, shuffle=NOSHUFFLE, blocksize=0)   [VLenUTF8()]
/variant_position      int32    36.37 KiB   39.93 KiB      1.1            11  3.63 KiB                 3.31 KiB            (10221,)            (1000,)           Blosc(cname='zstd',
clevel=7, shuffle=NOSHUFFLE, blocksize=0)   None
/region_index          int32    8.62 KiB    264 bytes      0.03            1  264.0 bytes              8.62 KiB            (11, 6)             (11, 6)           Blosc(cname='zstd',
clevel=9, shuffle=NOSHUFFLE, blocksize=0)   None
/variant_length        int8     5.17 KiB    9.98 KiB       1.9            11  929.1818181818181 bytes  481 bytes           (10221,)            (1000,)           Blosc(cname='zstd',
clevel=7, shuffle=NOSHUFFLE, blocksize=0)   None
/variant_contig        int8     4.99 KiB    9.98 KiB       2              11  929.1818181818181 bytes  464 bytes           (10221,)            (1000,)           Blosc(cname='zstd',
clevel=7, shuffle=NOSHUFFLE, blocksize=0)   None
/variant_quality       float32  4.93 KiB    39.93 KiB      8.1            11  3.63 KiB                 459 bytes           (10221,)            (1000,)           Blosc(cname='zstd',
clevel=7, shuffle=NOSHUFFLE, blocksize=0)   None
/variant_id_mask       bool     4.87 KiB    9.98 KiB       2              11  929.1818181818181 bytes  453 bytes           (10221,)            (1000,)           Blosc(cname='zstd',
clevel=7, shuffle=BITSHUFFLE, blocksize=0)  None
/contig_id             object   4.51 KiB    200 bytes      0.043           1  200.0 bytes              4.51 KiB            (25,)               (25,)             Blosc(cname='zstd',
clevel=7, shuffle=SHUFFLE, blocksize=0)     [VLenUTF8()]
/contig_length         int64    4.5 KiB     200 bytes      0.043           1  200.0 bytes              4.5 KiB             (25,)               (25,)             Blosc(cname='zstd',
clevel=7, shuffle=SHUFFLE, blocksize=0)     None
/filter_id             object   4.43 KiB    8 bytes        0.0018          1  8.0 bytes                4.43 KiB            (1,)                (1,)              Blosc(cname='zstd',
clevel=7, shuffle=SHUFFLE, blocksize=0)     [VLenUTF8()]
```

- each `.vcz` also has a CSV file created via `vcf2zarr inspect` in the `csvs` directory

### Compare data sizes

```
20G	    ofh.chr22-b0.vcf.gz
19G	    ofh.chr22-b1.vcf.gz
3.3G    ofh.chr22-b2.vcf.gz
41G	    ofh.chr22.icf
41G	    ofh.chr22.vcz
17G	    ofh.chr22.bitround.vcz
31G     ofh.chr22.quantize.vcz
```

### Compare distributions 

Compare default and filteed values via histograms of each format field

- call_LRR

![title](plots/default_call_LRR.png)
![title](plots/bitround_call_LRR.png)
![title](plots/quantize_call_LRR.png)

- call_BAF

![title](plots/default_call_BAF.png)
![title](plots/bitround_call_BAF.png)
![title](plots/quantize_call_BAF.png)

- call_GS

![title](plots/default_call_GS.png)
![title](plots/bitround_call_GS.png)
![title](plots/quantize_call_GS.png)

### Calculate Mean Absolute Error

Look at the mean absolute error by comparing default vs filtered for each format field

```
import sgkit
import dask.array as da

format_fields = ["call_LRR", "call_BAF", "call_GS"]

ds_1 = sgkit.load_dataset(zarr_path1)
ds_2 = sgkit.load_dataset(zarr_path2)

for format_field in format_fields:
    z1 = ds_1[format_field]
    z2 = ds_2[format_field]
    mae = da.nanmean(da.abs(z1 - z2)).compute()
```

| format field | filter type | mae | 
| --- | --- |  --- |
| call_LRR | bitround | 0.0006934664561413229 |
| call_BAF | bitround | 0.0007864997605793178 |
| call_GS | bitround | 0.003377074608579278 |
| call_LRR | quantize | 1.907108071463881e-06 |
| call_BAF | quantize | 1.1768814829338226e-06 |
| call_GS | quantize | 1.8836773278962937e-06 |
