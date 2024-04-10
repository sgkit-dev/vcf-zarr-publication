import os
import sys
import pathlib
import subprocess
import time
import dataclasses
import multiprocessing
import tempfile

import numba
import zarr
import psutil
import humanize
import numpy as np
import pandas as pd
import tskit
import click
import sgkit as sg


# yuck - but simplest way to avoid changing directory structure
sys.path.insert(0, "src")
from zarr_afdist import zarr_afdist, classify_genotypes_subset_filter, zarr_decode


def du(path):
    """
    Return the total disk usage of path in bytes.
    """
    out = subprocess.run(f"du -sb {path}", shell=True, check=True, capture_output=True)
    return int(out.stdout.decode().split()[0])


@dataclasses.dataclass
class ProcessTimeResult:
    wall: float
    system: float
    user: float


def time_cli_command(cmd, debug):
    # FIXME this doesn't look like it's capturing time
    # time spent by threads correctly

    # S: Total number of CPU-seconds used by the system on behalf of
    #    the process (in kernel mode), in seconds.
    # U: Total number of CPU-seconds that the process used directly
    #    (in user mode), in seconds.
    before = time.time()
    if debug:
        print(cmd)
    out = subprocess.run(cmd, shell=True, check=True, capture_output=True)
    wall_time = time.time() - before
    time_str = out.stderr.decode()
    sys_time, user_time = map(float, time_str.split())
    if debug:
        print(out.stdout.decode())
    return ProcessTimeResult(wall_time, sys_time, user_time)


def savvy_version():
    cmd = "./software/savvy/bin/sav --version"
    out = subprocess.run(cmd, shell=True, check=True, capture_output=True)
    return out.stdout.decode()


def bcftools_version():
    cmd = "./software/bcftools --version"
    out = subprocess.run(cmd, shell=True, check=True, capture_output=True)
    return out.stdout.decode()


def genozip_version():
    cmd = "./software/genozip --version"
    out = subprocess.run(cmd, shell=True, check=True, capture_output=True)
    return out.stdout.decode()


def variant_slice_coords(ds, variant_slice):
    pos = ds.variant_position[variant_slice].values
    return ds.contig_id.values[0], pos[0], pos[-1]


def get_variant_slice_region(ds, variant_slice):
    pos = ds.variant_position[variant_slice].values
    return f"{ds.contig_id.values[0]}:{pos[0]}-{pos[-1]}"


def write_sample_names(ds, sample_slice, f):
    sample_names = ds.sample_id[sample_slice].values
    for s in sample_names:
        print(s, file=f)
    f.flush()


def run_bcftools_afdist_subset(
    path, ds, variant_slice, sample_slice, *, num_threads=1, debug=False
):
    region = get_variant_slice_region(ds, variant_slice)
    with tempfile.NamedTemporaryFile("w") as f:
        write_sample_names(ds, sample_slice, f.file)
        cmd = (
            "export BCFTOOLS_PLUGINS=software/bcftools-1.18/plugins; "
            "/usr/bin/time -f'%S %U' "
            # Need to run the pipeline in a subshell to make sure we're
            # timing correctly
            "sh -c '"
            # bcftools view updates the AC and AN fields by default, but
            # it doesn't update AF (which +af-dist uses). So, we use
            # -I to prevent it updating, and then call fill-tags
            f"software/bcftools view -I -r {region} -S {f.name} "
            # Output uncompressed BCF to make pipeline more efficient
            f"-Ou --threads {num_threads} {path} | "
            f"software/bcftools +fill-tags -Ou --threads {num_threads} | "
            f"software/bcftools +af-dist --threads {num_threads}"
            "'"
        )
        return time_cli_command(cmd, debug)


def run_bcftools_afdist_filter(path, *, num_threads=1, debug=False):
    cmd = (
        "export BCFTOOLS_PLUGINS=software/bcftools-1.18/plugins; "
        "/usr/bin/time -f'%S %U' "
        # Need to run the pipeline in a subshell to make sure we're
        # timing correctly
        "sh -c '"
        # bcftools view updates the AC and AN fields by default, but
        # it doesn't update AF (which +af-dist uses). So, we use
        # -I to prevent it updating, and then call fill-tags
        f'software/bcftools view -I --include "FORMAT/DP>10 & FORMAT/GQ>20" '
        # Output uncompressed BCF to make pipeline more efficient
        f"-Ou --threads {num_threads} {path} | "
        f"software/bcftools +fill-tags -Ou --threads {num_threads} | "
        f"software/bcftools +af-dist --threads {num_threads}"
        "'"
    )
    return time_cli_command(cmd, debug)


# NOTE! These must be called on a file that has had fill-tags run on it.
def run_bcftools_afdist(path, *, num_threads=1, debug=False):
    cmd = (
        "BCFTOOLS_PLUGINS=software/bcftools-1.18/plugins "
        "/usr/bin/time -f'%S %U' "
        f"./software/bcftools +af-dist --threads {num_threads} {path}"
    )
    return time_cli_command(cmd, debug)


def run_genozip_afdist(path, *, num_threads=1, debug=False):
    cmd = (
        "export BCFTOOLS_PLUGINS=software/bcftools-1.18/plugins; "
        "/usr/bin/time -f'%S %U' "
        # Need to run the pipeline in a subshell to make sure we're
        # timing correctly
        "sh -c '"
        f"software/genocat --threads {num_threads} {path} | "
        f"software/bcftools +af-dist --threads {num_threads}"
        "'"
    )
    return time_cli_command(cmd, debug)


def run_genozip_afdist_subset(
    path, ds, variant_slice, sample_slice, *, num_threads=1, debug=False
):
    region = get_variant_slice_region(ds, variant_slice)
    # There's no "file" option for specifying samples with genozip,
    # so have to put on the command line. Hopefully this won't
    # break things
    samples = ",".join(ds.sample_id[sample_slice].values)
    cmd = (
        "export BCFTOOLS_PLUGINS=software/bcftools-1.18/plugins; "
        "/usr/bin/time -f'%S %U' "
        # Need to run the pipeline in a subshell to make sure we're
        # timing correctly
        "sh -c '"
        # f"software/bcftools view -I -r {region} -S {f.name} "
        f"software/genocat -r {region} -s {samples} "
        f"--threads {num_threads} {path} | "
        f"software/bcftools +fill-tags -Ou --threads {num_threads} | "
        f"software/bcftools +af-dist --threads {num_threads}"
        "'"
    )
    return time_cli_command(cmd, debug)


def run_savvy_afdist(path, *, debug=False):
    cmd = "/usr/bin/time -f'%S %U' " f"software/savvy-afdist/sav-afdist {path}"
    return time_cli_command(cmd, debug)


def run_savvy_decode(path, *, debug=False):
    cmd = (
        "/usr/bin/time -f'%S %U' "
        f"software/savvy-afdist/sav-afdist {path} --decode-only"
    )
    return time_cli_command(cmd, debug)


def run_savvy_afdist_subset(
    path, ds, variant_slice, sample_slice, *, num_threads=1, debug=False
):
    region_coords = variant_slice_coords(ds, variant_slice)
    assert region_coords[0] == "1"
    start, end = region_coords[1:]
    with tempfile.NamedTemporaryFile("w") as f:
        write_sample_names(ds, sample_slice, f.file)
        cmd = (
            "/usr/bin/time -f'%S %U' "
            f"software/savvy-afdist/sav-afdist --samples-file {f.name} "
            f"--start {start} --end {end} {path}"
        )
        return time_cli_command(cmd, debug)


def _zarr_afdist_subset_worker(ds_path, variant_slice, sample_slice, debug, conn):
    before = time.time()
    df = zarr_afdist(ds_path, variant_slice=variant_slice, sample_slice=sample_slice)
    wall_time = time.time() - before
    cpu_times = psutil.Process().cpu_times()
    if debug:
        print(df)
    conn.send(f"{wall_time} {cpu_times.user} {cpu_times.system}")


def zarr_afdist_worker(ds_path, debug, conn):
    return _zarr_afdist_subset_worker(ds_path, None, None, debug, conn)


def zarr_decode_worker(ds_path, debug, conn):
    before = time.time()
    bytes_decoded = zarr_decode(ds_path)
    wall_time = time.time() - before
    cpu_times = psutil.Process().cpu_times()
    if debug:
        print(bytes_decoded)
    conn.send(f"{wall_time} {cpu_times.user} {cpu_times.system}")


def zarr_afdist_subset_worker(ds_path, variant_slice, sample_slice, debug, conn):
    return _zarr_afdist_subset_worker(ds_path, variant_slice, sample_slice, debug, conn)


def run_zarr_afdist(ds_path, *, debug=False):
    conn1, conn2 = multiprocessing.Pipe()
    p = multiprocessing.Process(target=zarr_afdist_worker, args=(ds_path, debug, conn2))
    p.start()
    value = conn1.recv()
    wall_time, user_time, sys_time = map(float, value.split())
    p.join()
    if p.exitcode != 0:
        raise ValueError()
    p.close()
    return ProcessTimeResult(wall_time, sys_time, user_time)


def run_zarr_decode(ds_path, *, debug=False):
    conn1, conn2 = multiprocessing.Pipe()
    p = multiprocessing.Process(target=zarr_decode_worker, args=(ds_path, debug, conn2))
    p.start()
    value = conn1.recv()
    wall_time, user_time, sys_time = map(float, value.split())
    p.join()
    if p.exitcode != 0:
        raise ValueError()
    p.close()
    return ProcessTimeResult(wall_time, sys_time, user_time)


def run_zarr_afdist_subset(ds_path, ds, variant_slice, sample_slice, *, debug=False):
    conn1, conn2 = multiprocessing.Pipe()
    p = multiprocessing.Process(
        target=zarr_afdist_subset_worker,
        args=(ds_path, variant_slice, sample_slice, debug, conn2),
    )
    p.start()
    value = conn1.recv()
    wall_time, user_time, sys_time = map(float, value.split())
    p.join()
    if p.exitcode != 0:
        raise ValueError()
    p.close()
    return ProcessTimeResult(wall_time, sys_time, user_time)


@dataclasses.dataclass
class Tool:
    name: str
    suffix: str
    afdist_func: None
    afdist_subset_func: None
    version_func: None
    decode_func: None = None


all_tools = [
    Tool(
        "savvy",
        ".sav",
        run_savvy_afdist,
        run_savvy_afdist_subset,
        savvy_version,
        run_savvy_decode,
    ),
    Tool(
        "zarr", ".zarr", run_zarr_afdist, run_zarr_afdist_subset, None, run_zarr_decode
    ),
    # Making sure we run on the output of bcftools fill-tags
    Tool(
        "bcftools",
        ".tags.bcf",
        run_bcftools_afdist,
        run_bcftools_afdist_subset,
        bcftools_version,
    ),
    Tool(
        "genozip",
        ".tags.genozip",
        run_genozip_afdist,
        run_genozip_afdist_subset,
        genozip_version,
    ),
]


@click.command()
@click.argument("src", type=click.Path(), nargs=-1)
@click.option("-o", "--output", type=click.Path(), default=None)
@click.option("-t", "--tool", multiple=True, default=[t.name for t in all_tools])
@click.option("-s", "--storage", default="hdd")
@click.option("--debug", is_flag=True)
def whole_matrix_compute(src, output, tool, storage, debug):
    if len(src) == 0:
        raise ValueError("Need at least one input file!")
    tool_map = {t.name: t for t in all_tools}
    tools = [tool_map[tool_name] for tool_name in tool]

    data = []
    paths = [pathlib.Path(p) for p in sorted(src)]
    for ts_path in paths:
        ts = tskit.load(ts_path)
        click.echo(f"{ts_path} n={ts.num_samples // 2}, m={ts.num_sites}")
        sg_path = ts_path.with_suffix(".zarr")

        if not sg_path.exists:
            print("Skipping missing", sg_path)
            continue
        ds = sg.load_dataset(sg_path)
        num_sites = ts.num_sites
        assert ts.num_samples // 2 == ds.samples.shape[0]
        assert num_sites == ds.variant_position.shape[0]
        assert np.array_equal(ds.variant_position, ts.tables.sites.position.astype(int))

        for tool in tools:
            tool_path = ts_path.with_suffix(tool.suffix)
            if debug:
                print("Running:", tool)
            result = tool.afdist_func(tool_path, debug=debug)
            data.append(
                {
                    "num_samples": ds.samples.shape[0],
                    "num_sites": ts.num_sites,
                    "tool": tool.name,
                    "user_time": result.user,
                    "sys_time": result.system,
                    "wall_time": result.wall,
                    "storage": storage,
                }
            )
            df = pd.DataFrame(data).sort_values(["num_samples", "tool"])
            if output is not None:
                df.to_csv(output, index=False)
            print(df)


@click.command()
@click.argument("src", type=click.Path(), nargs=-1, required=True)
@click.option("-o", "--output", type=click.Path(), default=None)
@click.option("--debug", is_flag=True)
def file_size(src, output, debug):
    paths = [pathlib.Path(p) for p in sorted(src)]
    data = []
    for ts_path in paths:
        ts = tskit.load(ts_path)
        click.echo(f"{ts_path} n={ts.num_samples // 2}, m={ts.num_sites}")
        bcf_path = ts_path.with_suffix(".bcf")
        vcf_path = ts_path.with_suffix(".vcf.gz")
        # FIXME
        zarr_path = ts_path.with_suffix(".zarr")
        sav_path = ts_path.with_suffix(".sav")
        genozip_path = ts_path.with_suffix(".genozip")
        if not zarr_path.exists:
            print("Skipping missing", zarr_path)
            continue
        ds = sg.load_dataset(zarr_path)
        assert ts.num_samples // 2 == ds.samples.shape[0]
        assert ts.num_sites == ds.variant_position.shape[0]
        assert np.array_equal(ds.variant_position, ts.tables.sites.position.astype(int))
        tmap = {
            "tsk": ts_path,
            "vcf": vcf_path,
            "bcf": bcf_path,
            "zarr": zarr_path,
            "sav": sav_path,
            "genozip": genozip_path,
        }
        for tool, path in tmap.items():
            size = du(path)
            data.append(
                {
                    "sequence_length": int(ts.sequence_length),
                    "num_samples": ds.samples.shape[0],
                    "num_sites": ts.num_sites,
                    "tool": tool,
                    "size": size,
                }
            )
            df = pd.DataFrame(data).sort_values(["num_samples", "tool"])
            df.to_csv(output, index=False)
        print(df)


@click.command()
@click.argument("src", type=click.Path(), nargs=-1)
@click.option("-o", "--output", type=click.Path(), default=None)
@click.option("-t", "--tool", multiple=True, default=[t.name for t in all_tools[:2]])
@click.option("-s", "--storage", default="hdd")
@click.option("--debug", is_flag=True)
def whole_matrix_decode(src, output, tool, storage, debug):
    if len(src) == 0:
        raise ValueError("Need at least one input file!")
    tool_map = {t.name: t for t in all_tools[:2]}
    tools = [tool_map[tool_name] for tool_name in tool]

    data = []
    paths = [pathlib.Path(p) for p in sorted(src)]
    for ts_path in paths:
        ts = tskit.load(ts_path)
        click.echo(f"{ts_path} n={ts.num_samples // 2}, m={ts.num_sites}")
        sg_path = ts_path.with_suffix(".zarr")

        if not sg_path.exists:
            print("Skipping missing", sg_path)
            continue
        ds = sg.load_dataset(sg_path)
        num_sites = ts.num_sites
        num_samples = ds.samples.shape[0]
        assert ts.num_samples // 2 == ds.samples.shape[0]
        assert num_sites == ds.variant_position.shape[0]
        assert np.array_equal(ds.variant_position, ts.tables.sites.position.astype(int))

        for tool in tools:
            tool_path = ts_path.with_suffix(tool.suffix)
            if debug:
                print("Running:", tool)
            result = tool.decode_func(tool_path, debug=debug)
            data.append(
                {
                    "num_samples": ds.samples.shape[0],
                    "num_sites": ts.num_sites,
                    "tool": tool.name,
                    "user_time": result.user,
                    "sys_time": result.system,
                    "wall_time": result.wall,
                    "storage": storage,
                    "total_genotypes": num_sites * num_samples * 2,
                }
            )
            print(
                "rate = ",
                humanize.naturalsize(
                    data[-1]["total_genotypes"] / data[-1]["wall_time"], format="%.1f"
                ),
                "/ s",
            )
            df = pd.DataFrame(data).sort_values(["num_samples", "tool"])
            df.to_csv(output, index=False)
            print(df)


def midslice(n, k):
    """
    Return a slice of size k from the middle of an array of size n.
    """
    n2 = n // 2
    k2 = k // 2
    return slice(n2 - k2, n2 + k2)


@click.command()
@click.argument("src", type=click.Path(), nargs=-1)
@click.option("-o", "--output", type=click.Path(), default=None)
@click.option("-t", "--tool", multiple=True, default=[t.name for t in all_tools])
@click.option("-s", "--slice-id", multiple=True, default=["n10", "n/2"])
@click.option("--debug", is_flag=True)
def subset_processing_time(src, output, tool, slice_id, debug):
    if len(src) == 0:
        raise ValueError("Need at least one input file!")
    tool_map = {t.name: t for t in all_tools}
    tools = [tool_map[tool_name] for tool_name in tool]

    data = []
    paths = [pathlib.Path(p) for p in sorted(src)]
    for ts_path in paths:
        ts = tskit.load(ts_path)
        click.echo(f"{ts_path} n={ts.num_individuals}, m={ts.num_sites}")
        sg_path = ts_path.with_suffix(".zarr")

        if not sg_path.exists:
            print("Skipping missing", sg_path)
            continue
        ds = sg.load_dataset(sg_path)
        num_sites = ts.num_sites
        num_samples = ds.samples.shape[0]
        assert ts.num_samples // 2 == num_samples
        assert num_sites == ds.variant_position.shape[0]
        assert np.array_equal(ds.variant_position, ts.tables.sites.position.astype(int))

        slices = {
            "n10": (midslice(num_sites, 10000), midslice(num_samples, 10)),
            "n/2": (
                midslice(num_sites, 10000),
                midslice(num_samples, num_samples // 2),
            ),
        }

        for sid in slice_id:
            variant_slice, sample_slice = slices[sid]
            for tool in tools:
                tool_path = ts_path.with_suffix(tool.suffix)
                if debug:
                    print("Running:", tool)
                result = tool.afdist_subset_func(
                    tool_path,
                    ds,
                    variant_slice,
                    sample_slice,
                    debug=debug,
                )
                data.append(
                    {
                        "num_samples": ds.samples.shape[0],
                        "num_sites": ts.num_sites,
                        "slice": sid,
                        "tool": tool.name,
                        "user_time": result.user,
                        "sys_time": result.system,
                        "wall_time": result.wall,
                    }
                )
                df = pd.DataFrame(data).sort_values(["num_samples", "tool"])
                df.to_csv(output, index=False)
                print(df)


@click.command()
@click.argument("prefix", type=click.Path())
def genotype_filtering_processing_time(prefix):
    # run_bcftools_afdist_filter(prefix + ".vcf.gz", debug=True)

    # run_zarr_afdist_filter(prefix + ".zarr", debug=True)
    zarr_ds = zarr.open(prefix + ".zarr")
    counts = classify_genotypes_subset_filter(zarr_ds)
    print(counts)


@click.command()
@click.argument("prefix", type=click.Path())
def site_filtering_processing_time(prefix):
    # bcftools command:

    # bcftools query -i 'FILTER="PASS"' -f "%CHROM,%POS,%INFO/AN_EUR\n" tmp/1kg_chr20_all.vcf.gz > tmp1.csv
    zarr_ds = zarr.open(prefix + ".zarr")
    pass_filter = zarr_ds["variant_filter"][:, 0]
    contig_id = zarr_ds.contig_id[:]
    df = pd.DataFrame(
        {
            "CHROM": contig_id[zarr_ds["variant_contig"].vindex[pass_filter]],
            "POS": zarr_ds["variant_position"].vindex[pass_filter],
            "INFO/AN_EUR": zarr_ds["variant_AN_EUR"].vindex[pass_filter],
        }
    )
    before = time.perf_counter()
    df.to_csv("tmp2.csv", header=False, index=False)
    duration = time.perf_counter() - before
    print("csv writing took", duration, "for", len(df), " variants")


@click.command()
def report_versions():
    for tool in all_tools:
        print(tool.name)
        print(tool.version_func())


@click.group()
def cli():
    pass


cli.add_command(file_size)
cli.add_command(whole_matrix_compute)
cli.add_command(whole_matrix_decode)
cli.add_command(subset_processing_time)
cli.add_command(genotype_filtering_processing_time)
cli.add_command(site_filtering_processing_time)
cli.add_command(report_versions)


if __name__ == "__main__":
    cli()
