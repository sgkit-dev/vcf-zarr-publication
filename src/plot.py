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


sgkit_colour = "tab:blue"
bcf_colour = "tab:orange"
sav_colour = "tab:red"
genozip_colour = "tab:purple"
zarr_colour = "tab:blue"


def one_panel_fig(**kwargs):
    # The columnwidth of the format is ~250pt, which is
    # 3 15/32 inch, = 3.46
    width = 3.46
    fig, ax = plt.subplots(1, 1, figsize=(width, 2 * width / 3), **kwargs)
    return fig, ax


def plot_size(ax, df, label_y_offset=None):
    colour_map = {
        "vcf": "tab:pink",
        "bcf": bcf_colour,
        "zarr": zarr_colour,
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
        print(f"{tool} : {size:.2f}")
        ax.annotate(
            f"{size:.0f}G",
            textcoords="offset points",
            xytext=(15, label_y_offset.get(tool, 0)),
            xy=(row.num_samples, row["size"]),
            xycoords="data",
        )


    ax.legend()
    add_number_of_variants(df, ax)


def plot_total_cpu(ax, df, toolname=None, time_units="h"):
    colours = {
        "bcftools": bcf_colour,
        "genozip": genozip_colour,
        "zarr": zarr_colour,
        "savvy": sav_colour,
    }
    have_genozip = False
    toolname = {} if toolname is None else toolname
    divisors = {"h": 3600, "m": 60}

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

        hours = total_cpu[-1] // divisors[time_units]

        if tool != "genozip":
            ax.annotate(
                f"{hours:.0f}{time_units}",
                textcoords="offset points",
                xytext=(15, 0),
                xy=(row.num_samples, total_cpu[-1]),
                xycoords="data",
            )
        else:
            have_genozip = True

    # Extrapolate for genozip
    def mulplicative_model(n, a, b):
        # Fit a simple exponential function.
        return a * np.power(n, b)

    if have_genozip:
        dfs = df[(df.tool == "genozip")]
        fit_params, _ = optimize.curve_fit(
            mulplicative_model, dfs.num_samples[2:], dfs.wall_time[2:]
        )
        num_samples = df[(df.tool == "bcftools")].num_samples.values
        fit = mulplicative_model(num_samples, *fit_params)
        # print(fit)

        ax.loglog(num_samples[4:], fit[4:], linestyle=":", color="lightgray")
        hours = fit[-1] // 3600
        ax.annotate(
            f"{hours:.0f}h*",
            textcoords="offset points",
            xytext=(15, 0),
            xy=(num_samples[-1], fit[-1]),
            xycoords="data",
        )
    ax.legend()
    add_number_of_variants(df, ax)


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

    ax1.set_xlabel("Number of samples")
    ax1.set_ylabel("File size (bytes)")

    sav = df1[df1.tool == "sav"]
    zarr = df1[df1.tool == "zarr"]
    ratio = sav["size"].values / zarr["size"].values
    print("sav / zarr ratio:", ratio)
    # I tried putting an inset axis showing the ratio, but it was too small.
    # ax_inset = ax1.inset_axes([0.70, 0.1, 0.25, 0.25])
    # ax_inset.semilogx(sav["num_samples"], ratio)
    plt.tight_layout()
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
        "genozip": "genozip <FILE> | bcftools +af-dist",
        "zarr": "zarr-python API",
        "savvy": "savvy C++ API",
    }
    plot_total_cpu(ax1, df, name_map)

    ax1.set_xlabel("Number of samples")
    ax1.set_ylabel("Time (seconds)")

    # ax1.set_title(f"af-dist CPU time")
    ax1.legend()

    plt.tight_layout()
    plt.savefig(output)


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

    plot_total_cpu(ax1, df, time_units="m")
    df["genotypes_per_second"] = df["total_genotypes"] / df["wall_time"]
    for tool in ["zarr", "savvy"]:
        max_rate = df[df.tool == tool]["genotypes_per_second"].max()
        print(tool, humanize.naturalsize(max_rate, binary=True))

    ax1.set_xlabel("Number of samples")
    ax1.set_ylabel("Time (seconds)")
    ax1.legend()

    plt.tight_layout()
    plt.savefig(output)


def plot_subset_time(ax, df, extrapolate_genozip=False):
    colours = {
        "bcftools": bcf_colour,
        "genozip": genozip_colour,
        # "savvy": sav_colour,
        "zarr": zarr_colour,
    }

    label_map = {
        "bcftools": "bcftools pipeline",
        "genozip": "genozip + bcftools pipeline",
        "zarr": "zarr-python API",
    }

    for tool in colours.keys():
        dfs = df[(df.threads == 1) & (df.tool == tool)]
        total_cpu = dfs["user_time"].values + dfs["sys_time"].values
        n = dfs["num_samples"].values
        ax.loglog(
            n,
            total_cpu,
            label=label_map[tool],
            # linestyle=ls,
            marker=".",
            color=colours[tool],
        )

        # Show wall-time too. Pipeline nature of the bcftools and genozip
        # commands means that it automatically threads, even if we don't
        # really want it to.
        ax.loglog(
            n,
            dfs["wall_time"].values,
            # label=f"{tool}",
            linestyle=":",
            # marker=".",
            color=colours[tool],
        )
        row = dfs.iloc[-1]

        hours = total_cpu[-1]  # // 60

        ax.annotate(
            f"{hours:.0f}s",
            textcoords="offset points",
            xytext=(15, 0),
            xy=(row.num_samples, total_cpu[-1]),
            xycoords="data",
        )

    if extrapolate_genozip:
        # Extrapolate for genozip
        def mulplicative_model(n, a, b):
            # Fit a simple exponential function.
            return a * np.power(n, b)

        dfs = df[(df.threads == 1) & (df.tool == "genozip")]
        fit_params, _ = optimize.curve_fit(
            mulplicative_model, dfs.num_samples[2:], dfs.wall_time[2:]
        )
        num_samples = df[(df.threads == 1) & (df.tool == "bcftools")].num_samples.values
        fit = mulplicative_model(num_samples, *fit_params)
        # print(fit)

        ax.loglog(num_samples[3:], fit[3:], linestyle=":", color="lightgray")
        t = fit[-1]
        ax.annotate(
            f"{t:.0f}s*",
            textcoords="offset points",
            xytext=(15, 0),
            xy=(num_samples[-1], fit[-1]),
            xycoords="data",
        )

    add_number_of_variants(df, ax)


def run_subset_matrix_plot(data, output, subset, extrapolate_genozip=False):
    df = pd.read_csv(data, index_col=False).sort_values("num_samples")
    fig, ax1 = one_panel_fig()
    plot_subset_time(
        ax1, df[df.slice == subset], extrapolate_genozip=extrapolate_genozip
    )

    ax1.set_xlabel("Number of samples")
    ax1.set_ylabel("Time (seconds)")
    ax1.legend()

    plt.tight_layout()
    plt.savefig(output)


@click.command()
@click.argument("data", type=click.File("r"))
@click.argument("output", type=click.Path())
def subset_matrix_compute(data, output):
    """
    Plot the figure showing compute performance on subsets of matrix afdist.
    """
    run_subset_matrix_plot(data, output, "n10")


@click.command()
@click.argument("data", type=click.File("r"))
@click.argument("output", type=click.Path())
def subset_matrix_compute_supplemental(data, output):
    """
    Plot the figure showing compute performance on subsets of matrix afdist.
    """
    run_subset_matrix_plot(data, output, "n/2", True)


@click.group()
def cli():
    pass


cli.add_command(data_scaling)
cli.add_command(whole_matrix_compute)
cli.add_command(whole_matrix_decode)
cli.add_command(subset_matrix_compute)
cli.add_command(subset_matrix_compute_supplemental)


if __name__ == "__main__":
    cli()
