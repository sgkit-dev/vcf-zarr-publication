/* To any C++ programmers reading this file: apologies! It's obviously the work
 * of non C++ programmers, but hopefully the nastiness won't make it
 * slower. If so, please do let us know and we'll fix. */

#include <fstream>
#include <iostream>
#include <savvy/reader.hpp>
#include <string>
#include <thread>
#include <vector>

void classify_genotypes(savvy::reader& f, std::vector<double>& af,
                        std::vector<int>& hets, std::vector<int>& homs) {
    savvy::variant var;
    std::vector<int> geno;

    while (f >> var) {
        var.get_format("GT", geno);
        int total_alleles = 0;
        int alt_alleles = 0;
        int het_count = 0;
        int hom_count = 0;

        /* Note: this assumes binary alleles, and is therefore slightly
         * different to the Zarr version. This shouldn't make any difference
         * to the gross timings of the operations, though.
         */
        for (size_t i = 0; i < geno.size(); i += 2) {
            int allele1 = geno[i];
            int allele2 = geno[i + 1];
            total_alleles += 2;
            alt_alleles += allele1 + allele2;

            if (allele1 != allele2) {
                het_count++;
            } else if (allele1 == 1 && allele2 == 1) {
                hom_count++;
            }
        }

        double allele_freq =
            total_alleles > 0 ? static_cast<double>(alt_alleles) / total_alleles
                              : 0.0;
        af.push_back(allele_freq);
        hets.push_back(het_count);
        homs.push_back(hom_count);
    }
}

void decode(const std::string& filename) {
    savvy::reader f(filename);
    savvy::variant var;
    std::vector<int> geno;
    int count = 0;

    while (f >> var) {
        var.get_format("GT", geno);
        /* This is just to make sure the compiler doesn't somehow ignore the
         * get_format bit·
         */
        if (geno[0] >= 0) {
            count++;
        }
    }

    std::cout << "Decoded variants: " << count << std::endl;
}

void output_pos(const std::string& filename, bool mem_only) {
    savvy::reader f(filename);
    savvy::variant var;
    std::vector<int> pos;

    while (f >> var) {
        pos.push_back(var.position());
    }
    if (! mem_only) {
        for(unsigned int j = 0; j < pos.size(); j++) {
            std::cout << pos[j] << std::endl;
        }
    }
}

// Use a lambda expression for finding the bin index
auto findBinIndex = [](const std::vector<double>& bins, double value) -> int {
    if (value < bins.front()) return 0;
    if (value >= bins.back()) return static_cast<int>(bins.size()) - 2;

    auto it = std::lower_bound(bins.begin(), bins.end(), value);
    return static_cast<int>((it == bins.end() || *it != value)
                                ? it - bins.begin() - 1
                                : it - bins.begin());
};

// Function to check if a file exists
bool fileExists(const std::string& filename) {
    std::ifstream file(filename);
    return file.good();
}

int main(int argc, char* argv[]) {
    std::string filename;
    std::string samples_file;
    bool decode_only = false;
    bool pos_only = false;
    bool mem_only = false;
    int start = -1;
    int end = -1;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--decode-only") {
            decode_only = true;
        } else if (arg == "--pos-only") {
            pos_only = true;
        } else if (arg == "--mem-only") {
            mem_only = true;
        } else if (arg == "--start") {
            i++;
            start = std::stoi(argv[i]);
        } else if (arg == "--end") {
            i++;
            end = std::stoi(argv[i]);
        } else if (arg == "--samples-file") {
            i++;
            samples_file = argv[i];
        } else {
            filename = arg;
        }
    }

    if (filename.empty()) {
        std::cerr << "Error: filename argument is required" << std::endl;
        return 1;
    }

    if (!fileExists(filename)) {
        std::cerr << "Error: File '" << filename << "' does not exist"
                  << std::endl;
        return 1;
    }

    if ((start == -1) != (end == -1)) {
        std::cerr << "Error: --start and --end must be supplied together"
                  << std::endl;
        return 1;
    }

    if (decode_only) {
        decode(filename);
    } else if (pos_only) {
        output_pos(filename, mem_only);
    } else {
        std::vector<double> af;
        std::vector<int> hets;
        std::vector<int> homs;
        savvy::reader f(filename);

        if (start != -1) {
            // Hardcoding a CHROM to "1" here for our bencharks for simplicity
            f.reset_bounds({"1", (uint64_t)start, (uint64_t)end});
        }
        if (!samples_file.empty()) {
            if (!fileExists(samples_file)) {
                std::cerr << "Error: File '" << samples_file
                          << "' does not exist" << std::endl;
                return 1;
            }
            std::vector<std::string> samples;
            std::ifstream input;

            input.open(samples_file);
            for (std::string line; std::getline(input, line);) {
                samples.push_back(line);
            }
            input.close();

            std::cout << "Read " << samples.size() << " samples from file "
                      << samples_file << std::endl;

            f.subset_samples({samples.begin(), samples.end()});
        }

        classify_genotypes(f, std::ref(af), std::ref(hets), std::ref(homs));

        std::cout << "Number of variants: " << af.size() << std::endl;

        const int num_bins = 10;
        std::vector<double> bins(num_bins + 1);
        for (int i = 0; i <= num_bins; ++i) {
            bins[i] = i * (1.0 / num_bins);
        }
        bins.back() += 0.01;  // Adjust the last bin

        std::vector<int> het_bins(num_bins, 0), hom_bins(num_bins, 0),
            total_counts(num_bins, 0);

        int index = 0;
        for (auto& freq : af) {
            double pRA = 2 * freq * (1 - freq);
            double pAA = freq * freq;

            int binIndex = findBinIndex(bins, pRA);
            het_bins[binIndex] += hets[index];

            binIndex = findBinIndex(bins, pAA);
            hom_bins[binIndex] += homs[index];

            index++;
        }

        // Summing het and hom counts
        std::transform(het_bins.begin(), het_bins.end(), hom_bins.begin(),
                       total_counts.begin(), std::plus<int>());
        // Output the results
        for (int i = 0; i < num_bins; ++i) {
            std::cout << "Bin " << i << " (" << bins[i] << " - " << bins[i + 1]
                      << "): " << total_counts[i] << std::endl;
        }
    }
    return 0;
}
