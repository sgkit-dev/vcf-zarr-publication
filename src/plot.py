import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import humanize
import pandas as pd
import click
import numpy as np
import scipy.optimize as optimize

# Main text size is 9pt
plt.rcParams.update({"font.size": 7})
plt.rcParams.update({"legend.fontsize": 6})
plt.rcParams.update({"lines.markersize": 4})


bcf_colour = "tab:orange"
vcf_colour = "tab:green"
sav_colour = "tab:red"
genozip_colour = "tab:purple"
zarr_colour = "tab:blue"
zarr_nshf_colour = "tab:cyan"
zarr_lz4_nshf_colour = "tab:purple"
zarr_zip_colour = "tab:green"
two_bit_colour = "tab:pink"
zarr_java_colour = "tab:olive"
ts_py_colour = "tab:brown"
ts_cpp_colour = "tab:gray"


def one_panel_fig(**kwargs):
    # The columnwidth of the format is ~250pt, which is
    # 3 15/32 inch, = 3.46
    width = 3.46
    fig, ax = plt.subplots(1, 1, figsize=(width, 2 * width / 3), **kwargs)
    return fig, ax


def two_panel_fig(**kwargs):
    # The columnwidth of the format is ~250pt, which is
    # 3 15/32 inch, = 3.46
    width = 3.46
    fig, ax = plt.subplots(1, 2, figsize=(width, 2 * width / 3), **kwargs)
    return fig, ax


def plot_size(ax, df, label_y_offset=None, order=None):
    colour_map = {
        "2bit": two_bit_colour,
        "vcf": vcf_colour,
        "bcf": bcf_colour,
        "zarr": zarr_colour,
        # "zarr_nshf": zarr_nshf_colour,
        "sav": sav_colour,
        "genozip": genozip_colour,
    }

    GB = 2**30
    label_y_offset = {} if label_y_offset is None else label_y_offset

    for tool, colour in colour_map.items():
        dfs = df[df.tool == tool]
        dfs = dfs.sort_values("num_samples")
        ax.loglog(
            dfs["num_samples"].values,
            dfs["size"].values,
            ".-",
            color=colour,
            label="vcf.gz" if tool == "vcf" else tool,
        )
        row = dfs.iloc[-1]
        size = row["size"] / GB
        label = f"{size:.0f}G"
        if size > 100:
            size /= 1024
            label = f"{size:.1f}T"
        ax.annotate(
            label,
            textcoords="offset points",
            xytext=(15, label_y_offset.get(tool, 0)),
            xy=(row.num_samples, row["size"]),
            xycoords="data",
        )

    df_large = df[df.num_samples == 10**6].copy()
    df_large["size"] /= GB
    print(df_large)

    ax.legend()
    add_number_of_variants(df, ax)
    ax.set_xlabel("Number of samples")
    ax.set_ylabel("Storage size (bytes)")
    plt.tight_layout()


def plot_s3_network_throughput(ax, df):

    GB = 2**30

    x = df["processes"].values
    for col in ["c5.9xlarge", "c5n.9xlarge"]:
        y = df[col].values
        ax.loglog(
            x ,
            y, 
            marker=".",
            label=col,
            base=2,
        )
        argmax = np.argmax(y)
        ax.plot(x[argmax], y[argmax], marker="+", color="black")
        ax.annotate(
            f"{y[argmax] / GB:.1f} GiB/s",
            textcoords="offset points",
            xytext=(-8, -10),
            xy=(x[argmax], y[argmax]),
            xycoords="data",
        )

    ax.set_xticks(x)
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

    ax.legend()
    ax.set_xlabel("Number of processes")
    ax.set_ylabel("Data download rate (bytes/s)")
    plt.tight_layout()


def plot_s3_throughput(ax, df):

    GB = 2**30

    argmax = 6
    ax.loglog(
        df["processes"].values,
        df["throughput_decompress"].values,
        marker=".",
        label="Decompression only",
        base=2,
    )
    ax.plot(
        df["processes"][argmax],
        df["throughput_decompress"][argmax],
        marker="+",
        color="black",
    )
    maxval = df["throughput_decompress"][argmax]
    ax.annotate(
        f"{maxval / GB:.1f} GiB/s",
        textcoords="offset points",
        xytext=(-8, -10),
        xy=(df["processes"][argmax], df["throughput_decompress"][argmax]),
        xycoords="data",
    )
    ax.loglog(
        df["processes"].values,
        df["throughput_afdist"].values,
        marker=".",
        label="Compute af-dist",
        base=2,
    )
    argmax = 7
    ax.plot(
        df["processes"][argmax],
        df["throughput_afdist"][argmax],
        marker="+",
        color="black",
    )
    maxval = df["throughput_afdist"][argmax]
    ax.annotate(
        f"{maxval / GB:.1f} GiB/s",
        textcoords="offset points",
        xytext=(-8, -10),
        xy=(df["processes"][argmax], df["throughput_afdist"][argmax]),
        xycoords="data",
    )

    ax.set_xticks(df["processes"].values)
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    # ax.set_xticklabels(labels)

    ax.legend()
    ax.set_xlabel("Number of processes")
    ax.set_ylabel("Genotype data throughput (bytes/s)")
    plt.tight_layout()


def plot_total_cpu(
    ax,
    df,
    toolname=None,
    colours=None,
    time_units="h",
    extrapolate=None,
    order=None,
    label_y_offset=None,
    side_digits=None,
):
    if colours is None:
        colours = {
            "bcftools+vcf": vcf_colour,
            "bcftools": bcf_colour,
            "genozip": genozip_colour,
            "zarr": zarr_colour,
            "zarr_nshf": zarr_nshf_colour,
            "zarr_lz4_nshf": zarr_lz4_nshf_colour,
            "zarr.zip": zarr_zip_colour,
            "savvy": sav_colour,
            "zarr_java": zarr_java_colour,
            "ts_py": ts_py_colour,
            "ts_cpp": ts_cpp_colour,
        }
    have_genozip = False
    toolname = {} if toolname is None else toolname
    divisors = {"s": 1, "h": 3600, "m": 60}
    extrapolate = [] if extrapolate is None else extrapolate

    label_y_offset = {} if label_y_offset is None else label_y_offset
    if order is None:
        order = colours.keys()

    # for tool in df.tool.unique():
    for tool in order:
        dfs = df[(df.tool == tool)]
        if dfs.empty:
            continue
        total_cpu = dfs["user_time"].values + dfs["sys_time"].values
        ax.loglog(
            dfs["num_samples"].values,
            total_cpu,
            label=f"{toolname.get(tool, tool)}",
            # linestyle=ls,
            marker=".",
            color=colours[tool],
        )

        # Show wall-time too. Pipeline nature of the bcftools and genozip
        # commands means that it automatically threads, even if we don't
        # really want it to.
        ax.loglog(
            dfs["num_samples"].values,
            dfs["wall_time"].values,
            linestyle=":",
            # marker=".",
            color=colours[tool],
        )
        row = dfs.iloc[-1]

        if tool not in extrapolate:
            time = total_cpu[-1] / divisors[time_units]
            if side_digits is None:
                digits = 0 if time > 1 else 1
            else:
                digits = side_digits

            side_label = f"{time:.{digits}f}{time_units}"
            ax.annotate(
                side_label,
                textcoords="offset points",
                xytext=(15, label_y_offset.get(tool, 0)),
                xy=(row.num_samples, total_cpu[-1]),
                xycoords="data",
            )

    df_large = df[df.num_samples == 10**6].copy()
    df_large["total_time"] = (df["user_time"] + df["sys_time"]) / 3600
    print(df_large)

    def multiplicative_model(n, a, b):
        # Fit a simple exponential function.
        return a * np.power(n, b)

    for tool in extrapolate:
        dfs = df[(df.tool == tool)]
        fit_params, _ = optimize.curve_fit(
            multiplicative_model, dfs.num_samples[2:], dfs.wall_time[2:]
        )
        num_samples = df[(df.tool == "zarr")].num_samples.values
        fit = multiplicative_model(num_samples, *fit_params)
        # print(fit)

        ax.loglog(num_samples[3:], fit[3:], linestyle=":", color="lightgray")
        time = fit[-1] / divisors[time_units]
        ax.annotate(
            f"{time:.0f}{time_units}*",
            textcoords="offset points",
            xytext=(15, label_y_offset.get(tool, 0)),
            xy=(num_samples[-1], fit[-1]),
            xycoords="data",
        )
    ax.legend()
    add_number_of_variants(df, ax)
    ax.set_xlabel("Number of samples")
    ax.set_ylabel("Time (seconds)")
    plt.tight_layout()


def add_number_of_variants(df, ax):
    dfs = df[df["tool"] == "zarr"]
    num_samples = dfs["num_samples"].values
    num_sites = dfs["num_sites"].values

    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xscale("log")
    ax2.set_xlabel("Number of variants")
    ax2.set_xticks(num_samples)
    ax2.set_xticklabels([humanize.metric(m) for m in num_sites])


@click.command()
@click.argument("size_data", type=click.File("r"))
@click.argument("output", type=click.Path())
def data_scaling(size_data, output):
    """
    Plot the figure showing file size.
    """
    df1 = pd.read_csv(size_data, index_col=None).sort_values("num_samples")

    sav = df1[df1.tool == "sav"]
    zarr = df1[df1.tool == "zarr"]
    ratio = sav["size"].values / zarr["size"].values
    print("sav / zarr ratio:", ratio)

    plink_ish = []
    for _, row in zarr.iterrows():
        d = dict(row)
        d["tool"] = "2bit"
        d["size"] = row.num_samples * row.num_sites / 4
        plink_ish.append(d)

    df1 = pd.concat([df1, pd.DataFrame(plink_ish)])

    fig, ax1 = one_panel_fig()
    plot_size(ax1, df1, label_y_offset={"vcf": 4, "bcf": 1, "sav": -5.5, "genozip": -7})

    # I tried putting an inset axis showing the ratio, but it was too small.
    # ax_inset = ax1.inset_axes([0.70, 0.1, 0.25, 0.25])
    # ax_inset.semilogx(sav["num_samples"], ratio)

    plt.savefig(output)


@click.command()
@click.argument("data", type=click.File("r"))
@click.argument("output", type=click.Path())
def s3_throughput(data, output):
    """
    Plot the figure showing compute performance on whole-matrix afdist.
    """
    df = pd.read_csv(data, index_col=False).sort_values("processes")
    print(df)
    # Data in here is power-of-10 MB/s, so convert back to bytes per
    # second for consistency
    df["throughput_decompress"] *= 1000**2
    df["throughput_afdist"] *= 1000**2

    fig, ax1 = one_panel_fig()
    plot_s3_throughput(ax1, df)
    plt.savefig(output)


@click.command()
@click.argument("data", type=click.File("r"))
@click.argument("output", type=click.Path())
def s3_network_throughput(data, output):
    """
    Plot the figure showing compute performance on whole-matrix afdist.
    """

    # Data derived from S3-1-DownloadChunks-MP.ipynb notebook. Numbers there]
    # are in power-of-10 MB

    # x = [1, 2, 4, 8, 16, 32, 64]
    # data1 = np.array([
    #     159.0341786762204,
    #     349.88301727773677,
    #     660.1621287063576,
    #     1009.4825168651081,
    #     1253.8112536680967,
    #     1280.6974940263988,
    #     1274.494689412922,
    # ])
    # data2 = np.array([
    #     162.80624306237462,
    #     345.7732721919105,
    #     658.9358771498762,
    #     1062.419542658184,
    #     1972.057722720383,
    #     2541.4617122839895,
    #     2523.501053909922,
    # ])
    # MB = 1000 * 1000
    # df = pd.DataFrame({"processes": x, "c5.9xlarge": data1 * MB, "c5n.9xlarge":
    #     data2 * MB})
    # print(df)
    # df.to_csv("plot_data/gel-s3-network-throughput.csv", index=False)

    df = pd.read_csv(data, index_col=False).sort_values("processes")
    print(df)

    fig, ax1 = one_panel_fig()
    plot_s3_network_throughput(ax1, df)
    plt.savefig(output)


@click.command()
@click.argument("time_data", type=click.File("r"))
@click.argument("output", type=click.Path())
def whole_matrix_compute(time_data, output):
    """
    Plot the figure showing compute performance on whole-matrix afdist.
    """
    df = pd.read_csv(time_data, index_col=False).sort_values("num_samples")
    df = df[df.storage == "hdd"]

    fig, ax1 = one_panel_fig()
    name_map = {
        "bcftools": "bcftools +af-dist <BCF_FILE>",
        "bcftools+vcf": "bcftools +af-dist <VCF_FILE>",
        "genozip": "genocat <FILE> | bcftools +af-dist",
        "zarr": "zarr-python API",
        "savvy": "savvy C++ API",
    }
    plot_total_cpu(
        ax1,
        df,
        name_map,
        extrapolate=["genozip", "bcftools+vcf"],
        order=["genozip", "bcftools+vcf", "bcftools", "zarr", "savvy"],
    )

    plt.savefig(output)


@click.command()
@click.argument("time_data", type=click.File("r"))
@click.argument("output", type=click.Path())
def whole_matrix_compute_zarr_versions(time_data, output):
    """
    Plot the figure showing compute performance on whole-matrix afdist for the zarr implementations.
    """
    df = pd.read_csv(time_data, index_col=False).sort_values("num_samples")
    df = df[df.storage == "hdd"]

    fig, ax1 = one_panel_fig()
    name_map = {
        "zarr_java": "JZarr",
        "ts_py": "Tensorstore Python",
        "zarr": "Zarr-python",
        "ts_cpp": "Tensorstore C++",
    }
    plot_total_cpu(
        ax1,
        df,
        name_map,
        order=["zarr_java", "ts_py", "zarr", "ts_cpp"],
        label_y_offset={"ts_cpp": -9, "zarr": -5},
        side_digits=2,
    )
    plt.savefig(output)
    cpp = df[(df["tool"] == "ts_cpp") & (df["num_samples"] == 10**6)].user_time
    zp = df[(df["tool"] == "zarr") & (df["num_samples"] == 10**6)].user_time
    print("Cpp / Zarr python", cpp.values / zp.values)


@click.command()
@click.argument("time_data", type=click.File("r"))
@click.argument("output", type=click.Path())
def whole_matrix_decode(time_data, output):
    """
    Plot the figure showing raw decode performance on whole-matrix afdist.
    """
    df = pd.read_csv(time_data, index_col=False).sort_values("num_samples")
    df = df[df.storage == "hdd"]

    fig, ax1 = one_panel_fig()

    name_map = {
        "zarr": "Zarr (Zstd + BitShuffle)",
        "savvy": "Savvy",
        "zarr_nshf": "Zarr (Zstd)",
        "zarr_lz4_nshf": "Zarr (LZ4)",
    }
    plot_total_cpu(
        ax1,
        df,
        toolname=name_map,
        time_units="m",
        order=["zarr", "zarr_nshf", "zarr_lz4_nshf", "savvy"],
        label_y_offset={"savvy": -6, "zarr_lz4_nshf": -4},
    )
    df["genotypes_per_second"] = df["total_genotypes"] / df["user_time"]
    for tool in name_map.keys():
        max_rate = df[df.tool == tool]["genotypes_per_second"].max()
        print(tool, humanize.naturalsize(max_rate, binary=True))

    plt.savefig(output)

    df["cpu_time"] = df["user_time"] + df["sys_time"]
    df = df[((df.tool == "zarr") | (df.tool == "zarr.zip")) & (df.num_samples == 10**6)]
    dfs = df[["tool", "user_time", "sys_time", "cpu_time", "wall_time"]].copy()
    dfs["wall_time"] /= 60
    dfs["cpu_time"] /= 60
    print(dfs)


@click.command()
@click.argument("time_data", type=click.File("r"))
@click.argument("output", type=click.Path())
def column_extract(time_data, output):
    """
    Plot the figure showing time to extract the POS column
    """
    df = pd.read_csv(time_data, index_col=False).sort_values("num_samples")
    df_mem = df[df.destination == "memory"]
    df = df[df.destination == "file"]

    toolname = {
        "bcftools": "bcftools query",
        "savvy": "Savvy C++",
        "zarr": "Zarr + pandas to_csv",
    }
    fig, ax1 = one_panel_fig()
    plot_total_cpu(
        ax1,
        df,
        toolname=toolname,
        time_units="s",
        extrapolate=["bcftools"],
        order=["bcftools", "savvy", "zarr"],
    )
    plot_total_cpu(
        ax1,
        df_mem,
        colours={"zarr": "black"},
        toolname={"zarr": "Zarr (memory)"},
        time_units="s",
    )
    plt.savefig(output)


def run_subset_matrix_plot(data, output, subset, extrapolate, **kwargs):
    df = pd.read_csv(data, index_col=False).sort_values("num_samples")
    fig, ax1 = one_panel_fig()

    label_map = {
        "bcftools": "bcftools pipeline",
        "genozip": "genozip + bcftools pipeline",
        "zarr": "zarr-python API",
        "savvy": "savvy C++ API",
    }

    plot_total_cpu(
        ax1,
        df[df.slice == subset],
        toolname=label_map,
        time_units="s",
        extrapolate=extrapolate,
        order=["genozip", "bcftools", "savvy", "zarr"],
        **kwargs,
    )
    plt.savefig(output)


@click.command()
@click.argument("data", type=click.File("r"))
@click.argument("output", type=click.Path())
def subset_matrix_compute(data, output):
    """
    Plot the figure showing compute performance on subsets of matrix afdist.
    """
    run_subset_matrix_plot(
        data, output, "n10", extrapolate=[], label_y_offset={"savvy": -1}
    )


@click.command()
@click.argument("data", type=click.File("r"))
@click.argument("output", type=click.Path())
def subset_matrix_compute_supplemental(data, output):
    """
    Plot the figure showing compute performance on subsets of matrix afdist.
    """
    run_subset_matrix_plot(
        data,
        output,
        "n/2",
        extrapolate=["genozip"],
        label_y_offset={"savvy": 2, "zarr": -2.5},
    )


@click.command()
@click.argument("data", type=click.File("r"))
@click.argument("output", type=click.Path())
def compression_shuffle(data, output):
    """
    Plot figure showing the effect of shuffle settings on compression ratio.
    """
    df = pd.read_csv(data)

    # Note this is ordered by best-to-worst compression for viz

    arrays = [
        "call_GQ",
        "call_DP",
        "call_AD",
        "call_AB",
        "call_genotype",
    ]

    fig, ax = one_panel_fig()
    sns.barplot(
        df,
        orient="h",
        order=arrays,
        y="ArrayName",
        x="CompressionRatio",
        hue="Shuffle",
        ax=ax,
    )
    ax.set_ylabel("")
    ax.get_legend().set_title("")

    plt.tight_layout()
    plt.savefig(output)


@click.command()
@click.argument("data", type=click.File("r"))
@click.argument("output", type=click.Path())
def compression_compressor(data, output):
    """
    Plot figure showing the effect of compressor codec on compression ratio.
    """
    df = pd.read_csv(data)

    # Note this is ordered by best-to-worst compression for viz

    arrays = [
        "call_GQ",
        "call_DP",
        "call_AD",
        "call_AB",
        "call_genotype",
    ]

    fig, ax = one_panel_fig()
    sns.barplot(
        df,
        orient="h",
        order=arrays,
        y="ArrayName",
        x="CompressionRatio",
        hue="Compressor",
        ax=ax,
    )
    ax.set_ylabel("")
    ax.get_legend().set_title("")

    plt.tight_layout()
    plt.savefig(output)


@click.command()
@click.argument("data", type=click.File("r"))
@click.argument("output", type=click.Path())
def compression_chunksize(data, output):
    """
    Plot figure showing the effect of chunksize settings on compression ratio.
    """

    df = pd.read_csv(data)
    sample_df = df.loc[df.variant_chunksize == 10000]
    variant_df = df.loc[df.sample_chunksize == 1000]

    fig, axes = two_panel_fig()

    for arr in df.ArrayName.unique():
        arr_sdf = sample_df.loc[sample_df.ArrayName == arr].sort_values(
            "sample_chunksize"
        )
        arr_vdf = variant_df.loc[variant_df.ArrayName == arr].sort_values(
            "variant_chunksize"
        )
        axes[0].plot(
            arr_sdf.sample_chunksize, arr_sdf.CompressionRatio, label=arr, marker="o"
        )
        axes[1].semilogx(
            arr_vdf.variant_chunksize, arr_vdf.CompressionRatio, label=arr, marker="o"
        )

    plt.legend()
    axes[0].set_title("(A)")
    axes[1].set_title("(B)")
    axes[0].set_xlabel("Sample chunk size")
    axes[1].set_xlabel("Variant chunk size")
    axes[0].set_ylabel("Compression ratio")

    plt.tight_layout()
    plt.savefig(output)


@click.command()
@click.argument("data", type=click.File("r"))
@click.argument("output", type=click.Path())
def compression_chunksize_finegrained(data, output):
    """
    Plot figure showing the effect of finegrained chunksize settings on compression ratio.
    """

    df = pd.read_csv(data)
    fig, ax = one_panel_fig()

    markers = {"Odd": ".", "Even": "x", "Multiple of 4": "+"}

    remain_df = df.copy()

    norm = plt.Normalize(df.sample_chunksize.min(), df.sample_chunksize.max())

    # Plot each category with a different marker
    for category in ["Multiple of 4", "Even", "Odd"]:
        if category == "Even":
            cond = remain_df.sample_chunksize % 2 == 0
            sub_df = remain_df.loc[cond]
            remain_df = remain_df.loc[~cond]
        elif category == "Multiple of 4":
            cond = remain_df.sample_chunksize % 4 == 0
            sub_df = remain_df.loc[cond]
            remain_df = remain_df.loc[~cond]
        else:
            sub_df = remain_df

        s = ax.scatter(
            sub_df.sample_chunksize,
            sub_df.CompressionRatio,
            c=sub_df.smallest_chunk_size,
            marker=markers[category],
            label=category,
            cmap="plasma",
        )

    ax.set_xlabel("Sample Chunksize")
    ax.set_ylabel("Compression Ratio")
    plt.legend()
    leg = ax.get_legend()
    for handle in leg.legend_handles:
        handle.set_color("black")
    cbar = fig.colorbar(s)  # , cax=ax)
    cbar.ax.set_ylabel("Size of last chunk", rotation=270, labelpad=12)

    plt.tight_layout()
    plt.savefig(output)


@click.group()
def cli():
    pass


cli.add_command(data_scaling)
cli.add_command(s3_throughput)
cli.add_command(s3_network_throughput)
cli.add_command(whole_matrix_compute)
cli.add_command(whole_matrix_compute_zarr_versions)
cli.add_command(whole_matrix_decode)
cli.add_command(column_extract)
cli.add_command(subset_matrix_compute)
cli.add_command(subset_matrix_compute_supplemental)
cli.add_command(compression_shuffle)
cli.add_command(compression_chunksize)
cli.add_command(compression_compressor)
cli.add_command(compression_chunksize_finegrained)


if __name__ == "__main__":
    cli()
