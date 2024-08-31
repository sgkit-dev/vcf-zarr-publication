import numpy as np
import pandas as pd
import tensorstore as ts

def ts_afdist(path):
    path = str(path)
    store = ts.open({
            "driver": "zarr",
            "kvstore": {
                "driver": "file",
                "path": path,
            },
    }, write=False).result()
    variant_count = store.shape[0]
    sample_count = store.shape[1]
    chunk_shape = store.chunk_layout.read_chunk.shape
    variant_chunk_size = chunk_shape[0]
    sample_chunk_size = chunk_shape[1]
    bin_counts = np.zeros((11,), dtype=int)

    for variant_chunk_start in range(0, variant_count, variant_chunk_size):
        variant_chunk_end = min(variant_count, variant_chunk_start + variant_chunk_size)
        variant_chunk_len = variant_chunk_end - variant_chunk_start
        ref_counts = np.zeros((variant_chunk_len,), dtype=int)
        het_counts = np.zeros((variant_chunk_len,), dtype=int)
        hom_alt_counts = np.zeros((variant_chunk_len,), dtype=int)

        for sample_chunk_start in range(0, sample_count, sample_chunk_size):
            sample_chunk_end = min(sample_count, sample_chunk_start + sample_chunk_size)

            chunk = store[variant_chunk_start:variant_chunk_end, sample_chunk_start:sample_chunk_end].read().result()
            a = chunk[:, :, 0]
            b = chunk[:, :, 1]

            chunk_ref_counts = ((a == 0).astype(int) + (b == 0).astype(int)).sum(axis=1)
            chunk_het_counts = (a != b).sum(axis=1)
            chunk_hom_alt_counts = np.logical_and(a == b, a > 0).sum(axis=1)

            np.add(ref_counts, chunk_ref_counts, out=ref_counts)
            np.add(het_counts, chunk_het_counts, out=het_counts)
            np.add(hom_alt_counts, chunk_hom_alt_counts, out=hom_alt_counts)

        alt_count = 2 * sample_count - ref_counts
        alt_freq = alt_count / (2 * sample_count)
        het_ref_freq = 2 * alt_freq * (1 - alt_freq)
        hom_alt_freq = alt_freq * alt_freq

        bins = np.linspace(0, 1.0, len(bin_counts))
        bins[-1] += 0.0125
        a = np.bincount(np.digitize(het_ref_freq, bins), weights=het_counts, minlength=len(bins)).astype(int)
        b = np.bincount(np.digitize(hom_alt_freq, bins), weights=hom_alt_counts, minlength=len(bins)).astype(int)
        np.add(bin_counts, a, out=bin_counts)
        np.add(bin_counts, b, out=bin_counts)

    return pd.DataFrame({"start": bins[:-1], "stop": bins[1:], "prob dist": bin_counts[1:]})
