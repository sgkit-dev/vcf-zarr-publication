package dev.tskit.afdist;
import com.bc.zarr.*;
import java.util.*;

public class AFDist {
    private static final int VARIANT_DIM = 0;
    private static final int SAMPLE_DIM = 1;

    public static void main(String[] args) {
        if (args.length < 1) {
            System.err.println("Please provide the path to the Zarr file.");
            System.exit(1);
        }
        String zarrPath = args[0];
        
        try {
            run(zarrPath);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void run(String zarrPath) throws Exception {
        ZarrArray zarrArray = openZarrArray(zarrPath);
        int[] shape = zarrArray.getShape();
        int variantCount = shape[VARIANT_DIM];
        int sampleCount = shape[SAMPLE_DIM];

        int[] chunkShape = zarrArray.getChunks();
        int variantChunkSize = Math.min(chunkShape[VARIANT_DIM], variantCount);
        int sampleChunkSize = Math.min(chunkShape[SAMPLE_DIM], sampleCount);

        long[] binCounts = new long[11];
        processGenotypes(zarrArray, variantCount, sampleCount, variantChunkSize, sampleChunkSize, binCounts);
        
        printProbabilityDistribution(binCounts);
    }

    private static ZarrArray openZarrArray(String path) throws Exception {
        return ZarrArray.open(path);
    }

    private static void processGenotypes(ZarrArray store, int variantCount, int sampleCount,
                                         int variantChunkSize, int sampleChunkSize, long[] binCounts) throws Exception {
        byte[] dataBuffer = new byte[variantChunkSize * sampleChunkSize * 2];
        long[] refCounts = new long[variantChunkSize];
        long[] hetCounts = new long[variantChunkSize];
        long[] homAltCounts = new long[variantChunkSize];
        int[] chunk_shape = {variantChunkSize, sampleChunkSize, 2};

        for (int variantChunkStart = 0; variantChunkStart < variantCount; variantChunkStart += variantChunkSize) {
            int variantChunkEnd = Math.min(variantCount, variantChunkStart + variantChunkSize);
            int variantChunkLen = variantChunkEnd - variantChunkStart;

            Arrays.fill(refCounts, 0, variantChunkLen, 0);
            Arrays.fill(hetCounts, 0, variantChunkLen, 0);
            Arrays.fill(homAltCounts, 0, variantChunkLen, 0);

            for (int sampleChunkStart = 0; sampleChunkStart < sampleCount; sampleChunkStart += sampleChunkSize) {
                int sampleChunkEnd = Math.min(sampleCount, sampleChunkStart + sampleChunkSize);
                int sampleChunkLen = sampleChunkEnd - sampleChunkStart;

                int[] offset = {variantChunkStart, sampleChunkStart, 0};
                store.read(dataBuffer, chunk_shape, offset);

                for (int variantIndex = 0; variantIndex < variantChunkLen; variantIndex++) {
                    for (int sampleIndex = 0; sampleIndex < sampleChunkLen; sampleIndex++) {
                        int callIndex = 2 * (sampleChunkSize * variantIndex + sampleIndex);
                        byte a = dataBuffer[callIndex];
                        byte b = dataBuffer[callIndex + 1];

                        refCounts[variantIndex] += (a == 0 ? 1 : 0) + (b == 0 ? 1 : 0);
                        hetCounts[variantIndex] += (a != b ? 1 : 0);
                        homAltCounts[variantIndex] += (a == b && a > 0 ? 1 : 0);
                    }
                }
            }

            for (int variantIndex = 0; variantIndex < variantChunkLen; variantIndex++) {
                long refCount = refCounts[variantIndex];
                long hetCount = hetCounts[variantIndex];
                long homAltCount = homAltCounts[variantIndex];

                long altCount = 2L * sampleCount - refCount;
                double altFreq = altCount / (2.0 * sampleCount);
                double hetRefFreq = 2 * altFreq * (1 - altFreq);
                double homAltFreq = altFreq * altFreq;

                int binIndex = (int) (10 * hetRefFreq);
                binCounts[binIndex] += hetCount;

                binIndex = (int) (10 * homAltFreq);
                binCounts[binIndex] += homAltCount;
            }
        }
        binCounts[9] += binCounts[10];
        binCounts = Arrays.copyOf(binCounts, 10);        
    }

    private static void printProbabilityDistribution(long[] binCounts) {
        System.out.println("# PROB_DIST, genotype probability distribution, assumes HWE");
        for (int binIndex = 0; binIndex < binCounts.length; binIndex++) {
            double binStart = binIndex / 10.0;
            double binEnd = binStart + 0.1;
            System.out.printf("PROB_DIST\t%.1f\t%.1f\t%d%n", binStart, binEnd, binCounts[binIndex]);
        }
    }
}