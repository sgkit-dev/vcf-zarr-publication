import itertools
import zarr
import numcodecs
import os.path as osp
import pandas as pd
import numpy as np
from tqdm import tqdm
import argparse


def generate_chunksize_variations(z_arr):
    default_variant_chunksize = z_arr.chunks[0]
    default_sample_chunksize = z_arr.chunks[1]
    variant_chunksizes = np.logspace(2, np.log10(z_arr.shape[0]), num=6).astype(int)
    sample_chunksizes = np.linspace(100, z_arr.shape[1], num=6).astype(int)

    final_chunksizes = []
    for vcs in variant_chunksizes:
        final_chunksizes.append((vcs, default_sample_chunksize))

    for scs in sample_chunksizes:
        final_chunksizes.append((default_variant_chunksize, scs))

    return final_chunksizes


def test_compression_ratio_vs_chunksize(z_arr, dry_run=False):
    """
    Test the impact of chunksize on compression ratio

    :param z_arr: A Zarr array.
    :param dry_run: If True, generate the variations table without copying any data.
    """

    chunksize_var = generate_chunksize_variations(z_arr)
    comp_results_table = []

    if not dry_run:
        # Extract the data
        arr = z_arr[:]

    for cs_var in tqdm(chunksize_var, total=len(chunksize_var), desc="Chunksizes"):
        cs_var = (cs_var[0] or z_arr.shape[0], cs_var[1] or z_arr.shape[1])

        variant_cs, sample_cs = cs_var

        if len(z_arr.shape) == 3:
            cs_var += (z_arr.shape[2],)

        if not dry_run:
            if z_arr.filters is not None:
                object_codec = z_arr.filters[0]
            else:
                object_codec = None

            z2 = zarr.empty_like(
                arr,
                chunks=cs_var,
                compressor=z_arr.compressor,
                dtype=z_arr.dtype,
                object_codec=object_codec,
            )

            z2[:] = arr

            n_chunks = z2.nchunks
            compress_ratio = float(dict(z2.info_items())["Storage ratio"])
        else:
            n_chunks = None
            compress_ratio = None

        comp_results_table.append(
            {
                "chunksize": cs_var,
                "variant_chunksize": variant_cs,
                "sample_chunksize": sample_cs,
                "nchunks": n_chunks,
                "CompressionRatio": compress_ratio,
            }
        )

    return pd.DataFrame(comp_results_table).sort_values(
        "CompressionRatio", ascending=False
    )


def test_compression_ratio_vs_shuffle(z_arr, dry_run=False):
    shuffle_var = range(3)
    shuffle_names = ["No Shuffle", "Byte Shuffle", "Bit Shuffle"]

    comp_results_table = []

    if not dry_run:
        # Extract the data
        arr = z_arr[:]

    for sv, sn in tqdm(
        zip(shuffle_var, shuffle_names), desc="Testing shuffles", total=3
    ):
        if not dry_run:
            if z_arr.filters is not None:
                object_codec = z_arr.filters[0]
            else:
                object_codec = None

            new_compressor = numcodecs.Blosc(cname="zstd", shuffle=sv, clevel=7)

            z2 = zarr.empty_like(
                arr,
                chunks=z_arr.chunks,
                compressor=new_compressor,
                dtype=z_arr.dtype,
                object_codec=object_codec,
            )

            z2[:] = arr

            n_chunks = z2.nchunks
            compress_ratio = float(dict(z2.info_items())["Storage ratio"])
        else:
            n_chunks = None
            compress_ratio = None

        comp_results_table.append(
            {"Shuffle": sn, "nchunks": n_chunks, "CompressionRatio": compress_ratio}
        )

    return pd.DataFrame(comp_results_table).sort_values(
        "CompressionRatio", ascending=False
    )


def test_vcf2zarr_compression_variations(
    z_group,
    keys=("call_GQ", "call_DP", "call_AD", "call_AB", "call_genotype"),
    test_config="chunksize",
    dry_run=False,
):
    """
    Test the impact of compressor options + chunksizes on compression ratios
    of vcf2zarr zarr arrays.

    :param z_group: A Zarr group
    :param keys: A list of keys to test the compression variations on. If None,
    it tests on all the Zarr array in the hierarchy.
    :param dry_run: If True, generate the variations table without copying any data.


    """

    if keys is None:
        group_keys = list(z_group.keys())

        if not test_all:
            group_keys = [k for k in group_keys if "call_" in k]

        keys = group_keys

    tabs = []

    for k in keys:
        # Skip call_PL for now
        if k == "call_PL":
            continue

        print("Testing:", k)

        if test_config == "chunksize":
            res = test_compression_ratio_vs_chunksize(z_group[k], dry_run=dry_run)
        else:
            res = test_compression_ratio_vs_shuffle(z_group[k], dry_run=dry_run)

        tabs.append(res)
        tabs[-1]["ArrayName"] = k

    return pd.concat(tabs)


if __name__ == "__main__":
    # Home directory for the vcf-zarr-publication repo:
    home_dir = osp.dirname(osp.dirname(osp.abspath(__file__)))

    parser = argparse.ArgumentParser(description="vcf2zarr compression benchmarks")
    parser.add_argument(
        "-i",
        "--input",
        dest="input_zarr_group",
        type=str,
        required=True,
        help="The path to the input Zarr group.",
    )
    parser.add_argument(
        "--test-config",
        dest="test_config",
        type=str,
        default="chunksize",
        choices={"chunksize", "shuffle"},
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output_file",
        type=str,
        default=osp.join(home_dir, "plot_data/compression_benchmarks.csv"),
        help="The path to the output CSV file.",
    )
    parser.add_argument(
        "--dry-run",
        dest="dry_run",
        action="store_true",
        default=False,
        help="If True, generate the variations table without copying any data.",
    )

    args = parser.parse_args()

    results_df = test_vcf2zarr_compression_variations(
        zarr.open(args.input_zarr_group),
        test_config=args.test_config,
        dry_run=args.dry_run,
    )

    results_df["input_data"] = args.input_zarr_group

    results_df.to_csv(args.output_file, index=False)
