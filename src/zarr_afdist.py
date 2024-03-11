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
    n = call_genotype.shape[1]

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


def zarr_afdist(path, num_bins=10):
    root = zarr.open(path)
    call_genotype = root["call_genotype"]
    n = call_genotype.shape[1]

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
