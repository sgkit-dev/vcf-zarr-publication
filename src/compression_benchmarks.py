import itertools
import zarr
import numcodecs
import os.path as osp
import pandas as pd
import numpy as np
from tqdm import tqdm
import argparse


def generate_shape_permutations(shape):
    return list(itertools.permutations(range(len(shape))))


def generate_chunksize_variations(shape,
                                  min_chunksize=(1000, 100),
                                  max_chunksize=(10_000, None)):
    """
    Generate variations for the chunk sizes
    """

    if isinstance(min_chunksize, int):
        min_chunksize = list(np.repeat(min_chunksize, len(shape)))

    if isinstance(max_chunksize, int):
        max_chunksize = list(np.repeat(max_chunksize, len(shape)))

    dim_chunksizes = []

    for sdim, min_cs, max_cs in itertools.zip_longest(shape, min_chunksize, max_chunksize):

        if min_cs is None:
            min_cs = sdim

        if max_cs is None or max_cs > sdim:
            max_cs = sdim

        if sdim <= min_cs:
            dim_chunksizes.append([sdim])

        else:

            dim_chunksizes.append([min_cs])

            while True:
                if dim_chunksizes[-1][-1] * 10 <= max_cs:
                    dim_chunksizes[-1].append(dim_chunksizes[-1][-1] * 10)
                else:
                    break

    return list(itertools.product(*dim_chunksizes))


def generate_blosc_variations():
    """
    Generate variations for the Blosc compressor
    """

    from numcodecs.blosc import list_compressors

    compressors = ['zstd']  # list_compressors()
    bit_shuffle = list(range(3))  # AUTOSHUFFLE is discarded here
    c_levels = [7] #[5, 7, 9]

    return [
        {'cname': cn, 'shuffle': bits, 'clevel': cl}
        for cn, bits, cl in itertools.product(compressors, bit_shuffle, c_levels)
    ]


def test_compression_variations(z_arr,
                                max_chunks=10,
                                min_chunksize=(1000, 100),
                                max_chunksize=(10_000, None),
                                dry_run=False):
    """
    Test the impact of compressor options + chunksizes on compression ratios
    of Zarr arrays.

    :param z_arr: A Zarr array.
    :param max_chunks: Only extract up to `max_chunks` from the original array for processing
    (to save memory)
    :param max_chunksize: The maximum chunksize (along any dimension) to test.
    :param dry_run: If True, generate the variations table without copying any data.
    """

    shape_perm = generate_shape_permutations(z_arr.shape)
    chunksize_var = generate_chunksize_variations(z_arr.shape,
                                                  min_chunksize=min_chunksize,
                                                  max_chunksize=max_chunksize)
    blosc_var = generate_blosc_variations()

    comp_results_table = []

    # Determine the shape to extract based on the chunksize variations and
    # the max_chunks parameter:

    max_idx = np.minimum(np.array(chunksize_var).max(axis=0) * max_chunks,
                         np.array(z_arr.shape))
    print(max_idx)
    slices = slices = tuple(slice(0, n) for n in max_idx)

    if not dry_run:
        # Extract the data
        arr = z_arr[slices]

    for cs_var, bcomp_var, dim_order in tqdm(itertools.product(chunksize_var, blosc_var, shape_perm),
                                             total=len(chunksize_var) * len(blosc_var) * len(shape_perm),
                                             desc='Compression variations'):

        reshaped_chunks = tuple(np.array(cs_var)[np.array(dim_order)])
        
        variant_cs, sample_cs = cs_var[:2]

        named_chunksize = np.array([f'{k}{v}' for k, v in zip(['v', 's', 'p'], cs_var)])
        named_chunksize = named_chunksize[np.array(dim_order)]
        named_chunksize = ','.join(list(named_chunksize))

        if not dry_run:
            arr_reshaped = arr.transpose(dim_order)

            new_compressor = numcodecs.Blosc(**bcomp_var)

            if z_arr.filters is not None:
                object_codec = z_arr.filters[0]
            else:
                object_codec = None

            z2 = zarr.empty_like(arr_reshaped,
                                 chunks=reshaped_chunks,
                                 compressor=new_compressor,
                                 dtype=z_arr.dtype,
                                 object_codec=object_codec)

            z2[:] = arr_reshaped

            n_chunks = z2.nchunks
            compress_ratio = float(dict(z2.info_items())['Storage ratio'])
        else:
            n_chunks = None
            compress_ratio = None

        comp_results_table.append({**bcomp_var,
                                   **{'chunksize': reshaped_chunks,
                                      'named_chunksize': named_chunksize,
                                      'variant_chunksize': variant_cs,
                                      'sample_chunksize': sample_cs,
                                      'nchunks': n_chunks,
                                      'dim_order': dim_order,
                                      'CompressionRatio': compress_ratio}})

        # Test the combination of filters + compressor:
        if not dry_run:
            if arr.dtype == bool:
            
                z2 = zarr.empty_like(arr_reshaped,
                                 chunks=reshaped_chunks,
                                 compressor=new_compressor,
                                 filters=[numcodecs.PackBits()],
                                 dtype=z_arr.dtype,
                                 object_codec=object_codec)

                z2[:] = arr_reshaped

                n_chunks = z2.nchunks
                compress_ratio = float(dict(z2.info_items())['Storage ratio'])

                comp_results_table.append({**bcomp_var,
                                   **{'chunksize': reshaped_chunks,
                                      'named_chunksize': named_chunksize,
                                      'variant_chunksize': variant_cs,
                                      'sample_chunksize': sample_cs,
                                      'nchunks': n_chunks,
                                      'dim_order': dim_order,
                                      'CompressionRatio': compress_ratio}})
                comp_results_table[-1]['cname'] += '+PackBits'

            elif np.issubdtype(arr.dtype, np.floating):
            
                z2 = zarr.empty_like(arr_reshaped,
                                 chunks=reshaped_chunks,
                                 compressor=new_compressor,
                                 filters=[numcodecs.Quantize(1, arr.dtype)],
                                 dtype=z_arr.dtype,
                                 object_codec=object_codec)

                z2[:] = arr_reshaped

                n_chunks = z2.nchunks
                compress_ratio = float(dict(z2.info_items())['Storage ratio'])

                comp_results_table.append({**bcomp_var,
                                   **{'chunksize': reshaped_chunks,
                                      'named_chunksize': named_chunksize,
                                      'variant_chunksize': variant_cs,
                                      'sample_chunksize': sample_cs,
                                      'nchunks': n_chunks,
                                      'dim_order': dim_order,
                                      'CompressionRatio': compress_ratio}})
                comp_results_table[-1]['cname'] += '+Quantize'



    return pd.DataFrame(comp_results_table).sort_values('CompressionRatio', ascending=False)


def test_vcf2zarr_compression_variations(z_group,
                                         keys=None,
                                         max_chunks=10,
                                         min_chunksize=(1000, 100),
                                         max_chunksize=(10_000, None),
                                         test_all=False,
                                         dry_run=False):
    """
    Test the impact of compressor options + chunksizes on compression ratios
    of vcf2zarr zarr arrays.

    :param z_group: A Zarr group
    :param keys: A list of keys to test the compression variations on. If None,
    it tests on all the Zarr array in the hierarchy.
    :param max_chunks: Only extract up to `max_chunks` from the original array for processing
    (to save memory)
    :param max_chunksize: The maximum chunksize (along any dimension) to test.
    :param test_all: By default, we only test the call_* arrays. To test all arrays, set 
    this flag to True.
    :param dry_run: If True, generate the variations table without copying any data.


    """

    group_keys = list(z_group.keys())
    
    if not test_all:
        group_keys = [k for k in group_keys if 'call_' in k]

    keys = keys or group_keys

    tabs = []

    for k in keys:
        
        # Skip call_PL for now
        if k == 'call_PL':
            continue

        print("Testing:", k)
        tabs.append(test_compression_variations(z_group[k],
                                                max_chunks=max_chunks,
                                                min_chunksize=min_chunksize,
                                                max_chunksize=max_chunksize,
                                                dry_run=dry_run))
        tabs[-1]['ArrayName'] = k

    return pd.concat(tabs)


if __name__ == '__main__':

    # Home directory for the vcf-zarr-publication repo:
    home_dir = osp.dirname(osp.dirname(osp.abspath(__file__)))

    parser = argparse.ArgumentParser(description='vcf2zarr compression benchmarks')
    parser.add_argument('-i', '--input', dest='input_zarr_group', type=str, required=True,
                        help='The path to the input Zarr group.')
    parser.add_argument('-o', '--output', dest='output_file', type=str, 
                        default=osp.join(home_dir, 'plot_data/compression_benchmarks.csv'),
                        help='The path to the output CSV file.')
    parser.add_argument('--keys', dest='keys', type=str, nargs='+', default=None,
                        help='The keys to test the compression variations on. If not provided, '
                             'it tests on all the Zarr arrays in the hierarchy.')
    parser.add_argument('--test-all-arrays', dest='test_all', default=False, action='store_true',
                        help='By default, tests only `call_*` arrays. To test all arrays, set this flag.')
    parser.add_argument('--max-chunks', dest='max_chunks', type=int, default=10,
                        help='Only extract up to `max_chunks` from the original array for processing.')
    parser.add_argument('--dry-run', dest='dry_run', action='store_true', default=False,
                        help='If True, generate the variations table without copying any data.')

    args = parser.parse_args()

    results_df = test_vcf2zarr_compression_variations(zarr.open(args.input_zarr_group),
                                                      keys=args.keys,
                                                      max_chunks=args.max_chunks,
                                                      dry_run=args.dry_run,
                                                      test_all=args.test_all)
    results_df['input_data'] = args.input_zarr_group

    results_df.to_csv(args.output_file, index=False)

