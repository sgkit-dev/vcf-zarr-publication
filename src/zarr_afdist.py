import dataclasses

import numpy as np
import pandas as pd
import numba
import zarr


@numba.njit("void(int64, int8[:], int32[:], int32[:], int32[:], int32[:])")
def count_genotypes(index, g, hom_ref, hom_alt, het, ref_count):
    n = g.shape[0] // 2
    # NB Assuming no missing data!
    for i in range(n):
        j = 2 * i
        a = g[j]
        b = g[j + 1]
        if a == b:
            if a == 0:
                hom_ref[index] += 1
            else:
                hom_alt[index] += 1
        else:
            het[index] += 1
        ref_count[index] += (a == 0) + (b == 0)


@numba.njit("void(int64, int8[:, :, :], int32[:], int32[:], int32[:], int32[:])")
def count_genotypes_chunk(offset, G, hom_ref, hom_alt, het, ref_count):
    # NB Assuming diploids and no missing data!
    index = offset
    for j in range(G.shape[0]):
        for k in range(G.shape[1]):
            a = G[j, k, 0]
            b = G[j, k, 1]
            if a == b:
                if a == 0:
                    hom_ref[index] += 1
                else:
                    hom_alt[index] += 1
            else:
                het[index] += 1
            ref_count[index] += (a == 0) + (b == 0)
        index += 1


@numba.njit(
    "void(int64, int8[:, :, :], int8[:], int8[:], int32[:], int32[:], int32[:], int32[:])"
)
def count_genotypes_chunk_subset(
    offset, G, variant_mask, sample_mask, hom_ref, hom_alt, het, ref_count
):
    # NB Assuming diploids and no missing data!
    index = offset
    for j in range(G.shape[0]):
        if variant_mask[j]:
            for k in range(G.shape[1]):
                if sample_mask[k]:
                    a = G[j, k, 0]
                    b = G[j, k, 1]
                    if a == b:
                        if a == 0:
                            hom_ref[index] += 1
                        else:
                            hom_alt[index] += 1
                    else:
                        het[index] += 1
                    ref_count[index] += (a == 0) + (b == 0)
            index += 1


@dataclasses.dataclass
class GenotypeCounts:
    hom_ref: list
    hom_alt: list
    het: list
    ref_count: list


def classify_genotypes_variant_wise(call_genotype):
    m = call_genotype.shape[0]
    n = call_genotype.shape[1]

    het = np.zeros(m, dtype=np.int32)
    hom_alt = np.zeros(m, dtype=np.int32)
    hom_ref = np.zeros(m, dtype=np.int32)
    ref_count = np.zeros(m, dtype=np.int32)
    # This way is quite a bit slower, leading to substantially higher sys time.
    # Not clear why, since the IO should be done synchronously, but hey.
    for j, genotypes in enumerate(call_genotype):
        count_genotypes(j, genotypes.reshape(2 * n), hom_ref, hom_alt, het, ref_count)
    return GenotypeCounts(hom_ref, hom_alt, het, ref_count)


def classify_genotypes(call_genotype):
    m = call_genotype.shape[0]

    het = np.zeros(m, dtype=np.int32)
    hom_alt = np.zeros(m, dtype=np.int32)
    hom_ref = np.zeros(m, dtype=np.int32)
    ref_count = np.zeros(m, dtype=np.int32)
    j = 0
    for v_chunk in range(call_genotype.cdata_shape[0]):
        for s_chunk in range(call_genotype.cdata_shape[1]):
            G = call_genotype.blocks[v_chunk, s_chunk]
            count_genotypes_chunk(j, G, hom_ref, hom_alt, het, ref_count)
        j += G.shape[0]
    return GenotypeCounts(hom_ref, hom_alt, het, ref_count)


def classify_genotypes_subset(call_genotype, variant_mask, sample_mask):
    m = np.sum(variant_mask)

    # Use zarr arrays to get mask chunks aligned with the main data
    # for convenience.
    z_variant_mask = zarr.array(
        variant_mask, chunks=call_genotype.chunks[0], dtype=np.int8
    )
    z_sample_mask = zarr.array(
        sample_mask, chunks=call_genotype.chunks[1], dtype=np.int8
    )

    het = np.zeros(m, dtype=np.int32)
    hom_alt = np.zeros(m, dtype=np.int32)
    hom_ref = np.zeros(m, dtype=np.int32)
    ref_count = np.zeros(m, dtype=np.int32)
    j = 0
    # We should probably skip to the first non-zero chunk, but there probably
    # isn't much difference unless we have a huge number of chunks, and we're
    # only selecting a tiny subset
    for v_chunk in range(call_genotype.cdata_shape[0]):
        variant_mask_chunk = z_variant_mask.blocks[v_chunk]
        count = np.sum(variant_mask_chunk)
        if count > 0:
            for s_chunk in range(call_genotype.cdata_shape[1]):
                sample_mask_chunk = z_sample_mask.blocks[s_chunk]
                if np.sum(sample_mask_chunk) > 0:
                    G = call_genotype.blocks[v_chunk, s_chunk]
                    count_genotypes_chunk_subset(
                        j,
                        G,
                        variant_mask_chunk,
                        sample_mask_chunk,
                        hom_ref,
                        hom_alt,
                        het,
                        ref_count,
                    )
            j += count
    return GenotypeCounts(hom_ref, hom_alt, het, ref_count)


def zarr_afdist(path, num_bins=10):
    root = zarr.open(path)
    call_genotype = root["call_genotype"]
    m = call_genotype.shape[0]
    n = call_genotype.shape[1]

    # Using the more general code is slightly slower, 35s vs 30s on one
    # of the intermediate sized benchmarks.
    # counts = classify_genotypes_subset(
    #     call_genotype, np.ones(m, dtype=bool), np.ones(n, dtype=bool)
    # )
    counts = classify_genotypes(call_genotype)
    alt_count = 2 * n - counts.ref_count
    af = alt_count / (n * 2)
    # print("af", af)

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
