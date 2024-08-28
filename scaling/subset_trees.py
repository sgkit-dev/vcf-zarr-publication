import tskit
import os

ts_full = tskit.load("simulated_chrom_21.ts")
for k in range(1, 7):
    num_individuals = 10**k
    # Assuming diploids
    ts = ts_full.simplify(ts_full.samples()[: 2 * num_individuals])
    outfile = f"data/chr21_10_{k}.ts"

    os.makedirs(os.path.dirname(outfile), exist_ok=True)

    print("Writing", num_individuals, outfile)
    ts.dump(outfile)

