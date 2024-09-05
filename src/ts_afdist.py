import numpy as np
import pandas as pd
import tensorstore as ts

from zarr_afdist import count_genotypes_chunk, GenotypeCounts


def ts_afdist(path):
    path = str(path)
    store = ts.open(
        {
            "driver": "zarr",
            "kvstore": {
                "driver": "file",
                "path": path,
            },
            "context": {
                "cache_pool": {"total_bytes_limit": 0},
                "data_copy_concurrency": {"limit": 1},
                "file_io_concurrency": {"limit": 1},
            },
        },
        write=False,
    ).result()

    variant_count = store.shape[0]
    sample_count = store.shape[1]
    chunk_shape = store.chunk_layout.read_chunk.shape
    variant_chunk_size = chunk_shape[0]
    sample_chunk_size = chunk_shape[1]

    het = np.zeros(variant_count, dtype=np.int32)
    hom_alt = np.zeros(variant_count, dtype=np.int32)
    hom_ref = np.zeros(variant_count, dtype=np.int32)
    ref_count = np.zeros(variant_count, dtype=np.int32)

    for variant_chunk_start in range(0, variant_count, variant_chunk_size):
        variant_chunk_end = min(variant_count, variant_chunk_start + variant_chunk_size)

        for sample_chunk_start in range(0, sample_count, sample_chunk_size):
            sample_chunk_end = min(sample_count, sample_chunk_start + sample_chunk_size)

            G = (
                store[
                    variant_chunk_start:variant_chunk_end,
                    sample_chunk_start:sample_chunk_end,
                ]
                .read()
                .result()
            )
            count_genotypes_chunk(
                variant_chunk_start, G, hom_ref, hom_alt, het, ref_count
            )

    counts = GenotypeCounts(hom_ref, hom_alt, het, ref_count)

    num_bins = 10
    n = sample_count
    alt_count = 2 * n - counts.ref_count
    af = alt_count / (n * 2)
    bins = np.linspace(0, 1.0, num_bins + 1)
    bins[-1] += 0.0125
    pRA = 2 * af * (1 - af)
    pAA = af * af
    a = np.bincount(np.digitize(pRA, bins), weights=counts.het, minlength=num_bins + 1)
    b = np.bincount(
        np.digitize(pAA, bins), weights=counts.hom_alt, minlength=num_bins + 1
    )
    count = (a + b).astype(int)

    return pd.DataFrame({"start": bins[:-1], "stop": bins[1:], "prob_dist": count[1:]})
