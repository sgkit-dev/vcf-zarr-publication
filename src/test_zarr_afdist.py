import pytest
import zarr
import numpy as np
import pandas as pd
import numpy.testing as nt
import msprime
import sgkit as sg
import bio2zarr.vcf
import pysam


from . import zarr_afdist


# def _sgkit_afdist_subset_worker(
#     ds_path, variant_slice, sample_slice, num_threads, debug, conn
# ):
#     before = time.time()
#     with dask.distributed.Client(
#         processes=False, threads_per_worker=num_threads
#     ) as client:
#         ds = sg.load_dataset(ds_path)
#         if sample_slice is not None:
#             assert variant_slice is not None
#             ds = ds.isel(variants=variant_slice, samples=sample_slice)
#         df = get_prob_dist(ds)
#     wall_time = time.time() - before
#     cpu_times = psutil.Process().cpu_times()
#     if debug:
#         print(df)
#     conn.send(f"{wall_time} {cpu_times.user} {cpu_times.system}")


def sgkit_afdist(ds, num_bins=10):
    ds = sg.variant_stats(ds, merge=False).compute()
    bins = np.linspace(0, 1.0, num_bins + 1)
    bins[-1] += 0.0125
    af = 1 - ds.variant_allele_count[:, 0].values / ds.variant_allele_total.values
    pRA = 2 * af * (1 - af)
    pAA = af * af
    hets = ds.variant_n_het.values
    homs = ds.variant_n_hom_alt.values
    a = np.bincount(np.digitize(pRA, bins), weights=hets, minlength=num_bins + 1)
    b = np.bincount(np.digitize(pAA, bins), weights=homs, minlength=num_bins + 1)

    count = (a + b).astype(int)
    return pd.DataFrame({"start": bins[:-1], "stop": bins[1:], "prob_dist": count[1:]})


class TestSimulations:
    def run_simulation(
        self, tmp_path, n, L, seed, variant_chunk_size=None, sample_chunk_size=None
    ):
        ts = msprime.sim_ancestry(
            n, population_size=10_000, sequence_length=L, random_seed=seed
        )
        ts = msprime.sim_mutations(ts, rate=1e-8, random_seed=seed)
        assert ts.num_sites > 5

        vcf_file = tmp_path / "sim.vcf"
        with open(vcf_file, "w") as f:
            ts.write_vcf(f)
        # This also compresses the input file
        pysam.tabix_index(str(vcf_file), preset="vcf")
        out = tmp_path / "example.vcf.zarr"
        bio2zarr.vcf.convert(
            [tmp_path / "sim.vcf.gz"],
            out,
            chunk_length=variant_chunk_size,
            chunk_width=sample_chunk_size,
        )
        return out

    @pytest.mark.parametrize(
        ["n", "L"],
        [
            (8, 10**5),
            (8, 10**6),
        ],
    )
    @pytest.mark.parametrize("seed", range(1, 10))
    def test_afdist(self, tmp_path, n, L, seed):
        zarr_path = self.run_simulation(tmp_path, n, L, seed)
        afdist1 = zarr_afdist.zarr_afdist(zarr_path, 8)
        ds = sg.load_dataset(zarr_path)
        afdist2 = sgkit_afdist(ds, 8)
        nt.assert_array_equal(afdist1.start, afdist2.start)
        nt.assert_array_equal(afdist1.stop, afdist2.stop)
        nt.assert_array_equal(afdist1.prob_dist, afdist2.prob_dist)

    @pytest.mark.parametrize(
        ["n", "L"],
        [
            (10, 10**5),
            (10, 10**6),
            (100, 10**5),
            (1000, 10**5),
        ],
    )
    @pytest.mark.parametrize("seed", range(1, 10))
    def test_classify_genotypes(self, tmp_path, n, L, seed):
        zarr_path = self.run_simulation(
            tmp_path, n, L, seed, variant_chunk_size=200, sample_chunk_size=103
        )
        root = zarr.open(zarr_path)
        counts = zarr_afdist.classify_genotypes(root["call_genotype"])
        ds = sg.load_dataset(zarr_path)
        sg_counts = sg.variant_stats(ds, merge=False).compute()
        nt.assert_array_equal(sg_counts.variant_n_het, counts.het)
        nt.assert_array_equal(sg_counts.variant_n_hom_alt, counts.hom_alt)
        nt.assert_array_equal(sg_counts.variant_n_hom_ref, counts.hom_ref)

        counts = zarr_afdist.classify_genotypes_variant_wise(root["call_genotype"])
        nt.assert_array_equal(sg_counts.variant_n_het, counts.het)
        nt.assert_array_equal(sg_counts.variant_n_hom_alt, counts.hom_alt)
        nt.assert_array_equal(sg_counts.variant_n_hom_ref, counts.hom_ref)

    @pytest.mark.parametrize(
        ["n", "L"],
        [
            (10, 10**5),
            (10, 10**6),
            (100, 10**5),
            (1000, 10**5),
        ],
    )
    @pytest.mark.parametrize("seed", range(1, 10))
    def test_classify_genotypes_subset(self, tmp_path, n, L, seed):
        zarr_path = self.run_simulation(
            tmp_path, n, L, seed, variant_chunk_size=200, sample_chunk_size=103
        )
        root = zarr.open(zarr_path)
        G = root["call_genotype"]
        m = G.shape[0]
        n = G.shape[1]
        variants_mask = np.ones(m, dtype=bool)
        variants_mask[0 : m // 2] = 0
        samples_mask = np.ones(n, dtype=bool)
        samples_mask[n // 2 :] = 0
        counts = zarr_afdist.classify_genotypes_subset(G, variants_mask, samples_mask)
        ds = sg.load_dataset(zarr_path)
        ds = ds.isel(variants=variants_mask, samples=samples_mask)
        sg_counts = sg.variant_stats(ds, merge=False).compute()
        nt.assert_array_equal(sg_counts.variant_n_het, counts.het)
        nt.assert_array_equal(sg_counts.variant_n_hom_alt, counts.hom_alt)
        nt.assert_array_equal(sg_counts.variant_n_hom_ref, counts.hom_ref)

    @pytest.mark.parametrize(
        ["n", "L"],
        [
            (10, 10**5),
            (10, 10**6),
            (100, 10**5),
            (1000, 10**5),
        ],
    )
    @pytest.mark.parametrize("seed", range(11, 15))
    def test_af(self, tmp_path, n, L, seed):
        zarr_path = self.run_simulation(tmp_path, n, L, seed)
        counts = zarr_afdist.classify_genotypes(zarr.open(zarr_path)["call_genotype"])
        alt_count = 2 * n - counts.ref_count
        af = alt_count / (n * 2)
        ds = sg.load_dataset(zarr_path)
        sg_counts = sg.variant_stats(ds, merge=False).compute()
        sg_af = (
            1 - sg_counts.variant_allele_count[:, 0] / sg_counts.variant_allele_total
        )
        nt.assert_allclose(af, sg_af)
