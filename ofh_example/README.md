### Info

One genotype file taken from [release 9](https://ourfuturehealth.gitbook.io/our-future-health/data/data-releases/2024-data-releases/release-9).

- chr 22 (block 3)
- 814 variants
- 651050 samples

Format fileds of interest
- LRR: Log R Ratio (7 dp float)
- BAF: B allele frequency (7 dp float)
- GS: GenCall Score (4 dp float)


### Transform genotype VCF to column format 

```
INPUT_NAME=xxx

time vcf2zarr explode --force -p 20 -Q $INPUT_NAME.vcf.gz $INPUT_NAME.icf

real	1m27.363s
user	12m45.271s
sys	    0m37.283s
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

real	6m32.662s
user	6m32.180s
sys	    0m12.640s

name                   dtype    stored      size            ratio    nchunks  chunk_size    avg_chunk_stored    shape             chunk_shape       compressor                                                      filters
---------------------  -------  ----------  -----------  --------  ---------  ------------  ------------------  ----------------  ----------------  --------------------------------------------------------------  ------------
/call_LRR              float32  1.75 GiB    1.97 GiB       1.1           652  3.1 MiB       2.74 MiB            (814, 651050)     (10000, 1000)     Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   None
/call_BAF              float32  1.17 GiB    1.97 GiB       1.7           652  3.1 MiB       1.83 MiB            (814, 651050)     (10000, 1000)     Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   None
/call_GS               float32  190.35 MiB  1.97 GiB      11             652  3.1 MiB       298.95 KiB          (814, 651050)     (10000, 1000)     Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   None
/call_genotype         int8     58.47 MiB   1010.81 MiB   17             652  1.55 MiB      91.83 KiB           (814, 651050, 2)  (10000, 1000, 2)  Blosc(cname='zstd', clevel=7, shuffle=BITSHUFFLE, blocksize=0)  None
/call_genotype_mask    bool     10.76 MiB   1010.81 MiB   94             652  1.55 MiB      16.89 KiB           (814, 651050, 2)  (10000, 1000, 2)  Blosc(cname='zstd', clevel=7, shuffle=BITSHUFFLE, blocksize=0)  None
/sample_id             object   3.19 MiB    4.97 MiB       1.6           652  7.8 KiB       5.01 KiB            (651050,)         (1000,)           Blosc(cname='zstd', clevel=7, shuffle=SHUFFLE, blocksize=0)     [VLenUTF8()]
/call_genotype_phased  bool     555.95 KiB  505.4 MiB    930             652  793.76 KiB    873 bytes           (814, 651050)     (10000, 1000)     Blosc(cname='zstd', clevel=7, shuffle=BITSHUFFLE, blocksize=0)  None
/variant_allele        object   9.25 KiB    12.72 KiB      1.4             1  12.72 KiB     9.25 KiB            (814, 2)          (10000, 2)        Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   [VLenUTF8()]
/variant_filter        bool     8.5 KiB     814 bytes      0.093           1  814.0 bytes   8.5 KiB             (814, 1)          (10000, 1)        Blosc(cname='zstd', clevel=7, shuffle=BITSHUFFLE, blocksize=0)  None
/variant_id            object   7.96 KiB    6.36 KiB       0.8             1  6.36 KiB      7.96 KiB            (814,)            (10000,)          Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   [VLenUTF8()]
/variant_position      int32    7.11 KiB    3.18 KiB       0.45            1  3.18 KiB      7.11 KiB            (814,)            (10000,)          Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   None
/contig_id             object   4.51 KiB    200 bytes      0.043           1  200.0 bytes   4.51 KiB            (25,)             (25,)             Blosc(cname='zstd', clevel=7, shuffle=SHUFFLE, blocksize=0)     [VLenUTF8()]
/contig_length         int64    4.5 KiB     200 bytes      0.043           1  200.0 bytes   4.5 KiB             (25,)             (25,)             Blosc(cname='zstd', clevel=7, shuffle=SHUFFLE, blocksize=0)     None
/variant_quality       float32  4.47 KiB    3.18 KiB       0.71            1  3.18 KiB      4.47 KiB            (814,)            (10000,)          Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   None
/variant_contig        int8     4.46 KiB    814 bytes      0.18            1  814.0 bytes   4.46 KiB            (814,)            (10000,)          Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   None
/variant_id_mask       bool     4.46 KiB    814 bytes      0.18            1  814.0 bytes   4.46 KiB            (814,)            (10000,)          Blosc(cname='zstd', clevel=7, shuffle=BITSHUFFLE, blocksize=0)  None
/filter_id             object   4.43 KiB    8 bytes        0.0018          1  8.0 bytes     4.43 KiB            (1,)              (1,)              Blosc(cname='zstd', clevel=7, shuffle=SHUFFLE, blocksize=0)     [VLenUTF8()]

time vcf2zarr encode --force -p 20 -Qs schema_bitround.json $INPUT_NAME.icf $INPUT_NAME.bitround.vcz

real	11m46.174s
user	11m29.312s
sys	    0m29.800s

name                   dtype    stored      size            ratio    nchunks  chunk_size    avg_chunk_stored    shape             chunk_shape       compressor                                                      filters
---------------------  -------  ----------  -----------  --------  ---------  ------------  ------------------  ----------------  ----------------  --------------------------------------------------------------  ----------------------
/call_LRR              float32  682.05 MiB  1.97 GiB       3             652  3.1 MiB       1.05 MiB            (814, 651050)     (10000, 1000)     Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   [BitRound(keepbits=5)]
/call_BAF              float32  488.38 MiB  1.97 GiB       4.1           652  3.1 MiB       767.02 KiB          (814, 651050)     (10000, 1000)     Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   [BitRound(keepbits=5)]
/call_GS               float32  87.58 MiB   1.97 GiB      23             652  3.1 MiB       137.55 KiB          (814, 651050)     (10000, 1000)     Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   [BitRound(keepbits=5)]
/call_genotype         int8     58.47 MiB   1010.81 MiB   17             652  1.55 MiB      91.83 KiB           (814, 651050, 2)  (10000, 1000, 2)  Blosc(cname='zstd', clevel=7, shuffle=BITSHUFFLE, blocksize=0)  None
/call_genotype_mask    bool     10.76 MiB   1010.81 MiB   94             652  1.55 MiB      16.89 KiB           (814, 651050, 2)  (10000, 1000, 2)  Blosc(cname='zstd', clevel=7, shuffle=BITSHUFFLE, blocksize=0)  None
/sample_id             object   3.19 MiB    4.97 MiB       1.6           652  7.8 KiB       5.01 KiB            (651050,)         (1000,)           Blosc(cname='zstd', clevel=7, shuffle=SHUFFLE, blocksize=0)     [VLenUTF8()]
/call_genotype_phased  bool     555.95 KiB  505.4 MiB    930             652  793.76 KiB    873 bytes           (814, 651050)     (10000, 1000)     Blosc(cname='zstd', clevel=7, shuffle=BITSHUFFLE, blocksize=0)  None
/variant_allele        object   9.25 KiB    12.72 KiB      1.4             1  12.72 KiB     9.25 KiB            (814, 2)          (10000, 2)        Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   [VLenUTF8()]
/variant_filter        bool     8.5 KiB     814 bytes      0.093           1  814.0 bytes   8.5 KiB             (814, 1)          (10000, 1)        Blosc(cname='zstd', clevel=7, shuffle=BITSHUFFLE, blocksize=0)  None
/variant_id            object   7.96 KiB    6.36 KiB       0.8             1  6.36 KiB      7.96 KiB            (814,)            (10000,)          Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   [VLenUTF8()]
/variant_position      int32    7.11 KiB    3.18 KiB       0.45            1  3.18 KiB      7.11 KiB            (814,)            (10000,)          Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   None
/contig_id             object   4.51 KiB    200 bytes      0.043           1  200.0 bytes   4.51 KiB            (25,)             (25,)             Blosc(cname='zstd', clevel=7, shuffle=SHUFFLE, blocksize=0)     [VLenUTF8()]
/contig_length         int64    4.5 KiB     200 bytes      0.043           1  200.0 bytes   4.5 KiB             (25,)             (25,)             Blosc(cname='zstd', clevel=7, shuffle=SHUFFLE, blocksize=0)     None
/variant_quality       float32  4.47 KiB    3.18 KiB       0.71            1  3.18 KiB      4.47 KiB            (814,)            (10000,)          Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   None
/variant_contig        int8     4.46 KiB    814 bytes      0.18            1  814.0 bytes   4.46 KiB            (814,)            (10000,)          Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   None
/variant_id_mask       bool     4.46 KiB    814 bytes      0.18            1  814.0 bytes   4.46 KiB            (814,)            (10000,)          Blosc(cname='zstd', clevel=7, shuffle=BITSHUFFLE, blocksize=0)  None
/filter_id             object   4.43 KiB    8 bytes        0.0018          1  8.0 bytes     4.43 KiB            (1,)              (1,)              Blosc(cname='zstd', clevel=7, shuffle=SHUFFLE, blocksize=0)     [VLenUTF8()]

time vcf2zarr encode --force -p 20 -Qs schema_quantize.json $INPUT_NAME.icf $INPUT_NAME.quantize.vcz

real	8m47.324s
user	8m9.471s
sys	    0m49.021s

name                   dtype    stored      size            ratio    nchunks  chunk_size    avg_chunk_stored    shape             chunk_shape       compressor                                                      filters
---------------------  -------  ----------  -----------  --------  ---------  ------------  ------------------  ----------------  ----------------  --------------------------------------------------------------  ---------------------------------
/call_LRR              float32  1.3 GiB     1.97 GiB       1.5           652  3.1 MiB       2.04 MiB            (814, 651050)     (10000, 1000)     Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   [Quantize(digits=5, dtype='<f4')]
/call_BAF              float32  911.2 MiB   1.97 GiB       2.2           652  3.1 MiB       1.4 MiB             (814, 651050)     (10000, 1000)     Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   [Quantize(digits=5, dtype='<f4')]
/call_GS               float32  165.62 MiB  1.97 GiB      12             652  3.1 MiB       260.12 KiB          (814, 651050)     (10000, 1000)     Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   [Quantize(digits=5, dtype='<f4')]
/call_genotype         int8     58.47 MiB   1010.81 MiB   17             652  1.55 MiB      91.83 KiB           (814, 651050, 2)  (10000, 1000, 2)  Blosc(cname='zstd', clevel=7, shuffle=BITSHUFFLE, blocksize=0)  None
/call_genotype_mask    bool     10.76 MiB   1010.81 MiB   94             652  1.55 MiB      16.89 KiB           (814, 651050, 2)  (10000, 1000, 2)  Blosc(cname='zstd', clevel=7, shuffle=BITSHUFFLE, blocksize=0)  None
/sample_id             object   3.19 MiB    4.97 MiB       1.6           652  7.8 KiB       5.01 KiB            (651050,)         (1000,)           Blosc(cname='zstd', clevel=7, shuffle=SHUFFLE, blocksize=0)     [VLenUTF8()]
/call_genotype_phased  bool     555.95 KiB  505.4 MiB    930             652  793.76 KiB    873 bytes           (814, 651050)     (10000, 1000)     Blosc(cname='zstd', clevel=7, shuffle=BITSHUFFLE, blocksize=0)  None
/variant_allele        object   9.25 KiB    12.72 KiB      1.4             1  12.72 KiB     9.25 KiB            (814, 2)          (10000, 2)        Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   [VLenUTF8()]
/variant_filter        bool     8.5 KiB     814 bytes      0.093           1  814.0 bytes   8.5 KiB             (814, 1)          (10000, 1)        Blosc(cname='zstd', clevel=7, shuffle=BITSHUFFLE, blocksize=0)  None
/variant_id            object   7.96 KiB    6.36 KiB       0.8             1  6.36 KiB      7.96 KiB            (814,)            (10000,)          Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   [VLenUTF8()]
/variant_position      int32    7.11 KiB    3.18 KiB       0.45            1  3.18 KiB      7.11 KiB            (814,)            (10000,)          Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   None
/contig_id             object   4.51 KiB    200 bytes      0.043           1  200.0 bytes   4.51 KiB            (25,)             (25,)             Blosc(cname='zstd', clevel=7, shuffle=SHUFFLE, blocksize=0)     [VLenUTF8()]
/contig_length         int64    4.5 KiB     200 bytes      0.043           1  200.0 bytes   4.5 KiB             (25,)             (25,)             Blosc(cname='zstd', clevel=7, shuffle=SHUFFLE, blocksize=0)     None
/variant_quality       float32  4.47 KiB    3.18 KiB       0.71            1  3.18 KiB      4.47 KiB            (814,)            (10000,)          Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   None
/variant_contig        int8     4.46 KiB    814 bytes      0.18            1  814.0 bytes   4.46 KiB            (814,)            (10000,)          Blosc(cname='zstd', clevel=7, shuffle=NOSHUFFLE, blocksize=0)   None
/variant_id_mask       bool     4.46 KiB    814 bytes      0.18            1  814.0 bytes   4.46 KiB            (814,)            (10000,)          Blosc(cname='zstd', clevel=7, shuffle=BITSHUFFLE, blocksize=0)  None
/filter_id             object   4.43 KiB    8 bytes        0.0018          1  8.0 bytes     4.43 KiB            (1,)              (1,)              Blosc(cname='zstd', clevel=7, shuffle=SHUFFLE, blocksize=0)     [VLenUTF8()]
```

- each `.vcz` also has a CSV file created via `vcf2zarr inspect` in the `csvs` directory

### Compare data sizes

```
3.3G	xxx.vcf.gz
3.2G	xxx.icf
3.2G	xxx.vcz
1.4G	xxx.bitround.vcz
2.5G	xxx.quantize.vcz
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
| call_LRR | bitround | 0.0007233350770547986 |
| call_BAF | bitround | 0.0008083758875727654 |
| call_GS | bitround | 0.0033126999624073505 |
| call_LRR | quantize | 1.9064399339185911e-06 |
| call_BAF | quantize | 1.1348728321536328e-06 |
| call_GS | quantize | 1.925799779201043e-06 |
