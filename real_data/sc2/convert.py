import sys

sys.path.append("/home/jk/work/github/sc2ts")

# import pyfaidx
import sc2ts
import tqdm
import numpy as np
import zarr
import numcodecs

DEFAULT_ZARR_COMPRESSOR = numcodecs.Blosc(cname="zstd", clevel=7)


def encode_alignment(h):
    # FIXME Map anything that's not ACGT- to N.
    a = np.full(h.shape, -1, dtype=np.int8)
    for code, char in enumerate("ACGT-"):
        a[h == char] = code
    return a


def fa2zarr(fasta_path, zarr_path):
    root = zarr.open(zarr_path, mode="w")
    L = 29903
    with sc2ts.AlignmentStore(fasta_path) as store:
        N = len(store)
        sample_chunk_size = 1000
        shape = ((L, N),)
        gt_array = root.empty(
            "call_genotype",
            shape=(L, N),
            chunks=(L, sample_chunk_size),
            dtype=np.int8,
            compressor=numcodecs.Blosc(
                cname="zstd", clevel=7, shuffle=numcodecs.Blosc.BITSHUFFLE
            ),
            # filters=[numcodecs.get_codec(filt) for filt in array_spec.filters],
            # object_codec=object_codec,
            dimension_separator="/",
        )

        chunk = np.zeros_like(gt_array, shape=(L, sample_chunk_size)).T
        k = 0
        j = 0
        bar = tqdm.tqdm(total=N)
        for key, alignment in store.items():
            bar.update()
            # print(key, len(alignment), "\t", alignment)
            a = encode_alignment(alignment[1:])
            chunk[j] = a
            j += 1
            if j == sample_chunk_size:
                gt_array[:, k : k + sample_chunk_size] = chunk.T
                j = 0
                k += sample_chunk_size


if __name__ == "__main__":
    fa2zarr(sys.argv[1], sys.argv[2])
