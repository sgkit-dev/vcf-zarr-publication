import tskit

ts = tskit.load("simulated_chrom_21.ts")
for k in range(1, 7):
    num_individuals = 10**k
    # Assuming diploids
    ts = ts.simplify(ts.samples()[: 2 * num_individuals])
    outfile = f"data/chr21_10_{k}.ts"
    print("Writing", num_individuals, outfile)
    ts.dump(outfile)

