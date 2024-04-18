import matplotlib
import matplotlib.pyplot as plt
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
vcf_colour = "tab:pink"
sav_colour = "tab:red"
genozip_colour = "tab:purple"
zarr_colour = "tab:blue"
zarr_nshf_colour = "tab:cyan"


def one_panel_fig(**kwargs):
    # The columnwidth of the format is ~250pt, which is
    # 3 15/32 inch, = 3.46
    width = 3.46
    fig, ax = plt.subplots(1, 1, figsize=(width, 2 * width / 3), **kwargs)
    return fig, ax


def plot_size(ax, df, label_y_offset=None):
    colour_map = {
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
        ax.annotate(
            f"{size:.0f}G",
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
    plt.tight_layout()


def plot_total_cpu(ax, df, toolname=None, time_units="h", extrapolate=None):
    colours = {
        "bcftools+vcf": vcf_colour,
        "bcftools": bcf_colour,
        "genozip": genozip_colour,
        "zarr": zarr_colour,
        "zarr_nshf": zarr_nshf_colour,
        "savvy": sav_colour,
    }
    have_genozip = False
    toolname = {} if toolname is None else toolname
    divisors = {"s": 1, "h": 3600, "m": 60}
    extrapolate = [] if extrapolate is None else extrapolate

    # for tool in df.tool.unique():
    for tool in colours.keys():
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
            ax.annotate(
                f"{time:.0f}{time_units}",
                textcoords="offset points",
                xytext=(15, 0),
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
            xytext=(15, 0),
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

    fig, ax1 = one_panel_fig()
    plot_size(ax1, df1, label_y_offset={"vcf": 4, "sav": -5.5, "genozip": -7})

    sav = df1[df1.tool == "sav"]
    zarr = df1[df1.tool == "zarr"]
    ratio = sav["size"].values / zarr["size"].values
    print("sav / zarr ratio:", ratio)
    # I tried putting an inset axis showing the ratio, but it was too small.
    # ax_inset = ax1.inset_axes([0.70, 0.1, 0.25, 0.25])
    # ax_inset.semilogx(sav["num_samples"], ratio)

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
        "genozip": "genozip <FILE> | bcftools +af-dist",
        "zarr": "zarr-python API",
        "savvy": "savvy C++ API",
    }
    plot_total_cpu(ax1, df, name_map, extrapolate=["genozip", "bcftools+vcf"])

    plt.savefig(output)


@click.command()
@click.argument("time_data", type=click.File("r"))
@click.argument("output", type=click.Path())
def whole_matrix_decode(time_data, output):
    """
    Plot the figure showing raw decode performance on whole-matrix afdist.
    """
    df = pd.read_csv(time_data, index_col=False).sort_values("num_samples")
    print(df)
    df = df[df.storage == "hdd"]

    fig, ax1 = one_panel_fig()

    name_map = {
        "zarr": "Zarr (Zstd + bit shuffle)",
        "savvy": "Savvy",
        "zarr_nshf": "Zarr (Zstd)",
    }
    plot_total_cpu(ax1, df, toolname=name_map, time_units="m")
    df["genotypes_per_second"] = df["total_genotypes"] / df["user_time"]
    for tool in name_map.keys():
        max_rate = df[df.tool == tool]["genotypes_per_second"].max()
        print(tool, humanize.naturalsize(max_rate, binary=True))

    plt.savefig(output)


@click.command()
@click.argument("time_data", type=click.File("r"))
@click.argument("output", type=click.Path())
def column_extract(time_data, output):
    """
    Plot the figure showing time to extract the POS column
    """
    df = pd.read_csv(time_data, index_col=False).sort_values("num_samples")
    df = df[df.storage == "hdd"]

    toolname = {
        "bcftools": "bcftools query",
        "savvy": "savvy C++",
        "zarr": "zarr + pandas to_csv",
    }
    fig, ax1 = one_panel_fig()
    plot_total_cpu(ax1, df, toolname=toolname, time_units="s", extrapolate=["bcftools"])
    plt.savefig(output)


def run_subset_matrix_plot(data, output, subset, extrapolate):
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
    )
    plt.savefig(output)


@click.command()
@click.argument("data", type=click.File("r"))
@click.argument("output", type=click.Path())
def subset_matrix_compute(data, output):
    """
    Plot the figure showing compute performance on subsets of matrix afdist.
    """
    run_subset_matrix_plot(data, output, "n10", extrapolate=[])


@click.command()
@click.argument("data", type=click.File("r"))
@click.argument("output", type=click.Path())
def subset_matrix_compute_supplemental(data, output):
    """
    Plot the figure showing compute performance on subsets of matrix afdist.
    """
    run_subset_matrix_plot(data, output, "n/2", extrapolate=["genozip"])


@click.group()
def cli():
    pass


cli.add_command(data_scaling)
cli.add_command(whole_matrix_compute)
cli.add_command(whole_matrix_decode)
cli.add_command(column_extract)
cli.add_command(subset_matrix_compute)
cli.add_command(subset_matrix_compute_supplemental)


if __name__ == "__main__":
    cli()
