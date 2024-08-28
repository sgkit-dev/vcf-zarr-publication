
#include <iostream>

#include "absl/status/status.h"
#include "tensorstore/array.h"
#include "tensorstore/chunk_layout.h"
#include "tensorstore/index_space/dim_expression.h"
#include "tensorstore/open.h"
#include "tensorstore/spec.h"
#include "tensorstore/tensorstore.h"

using ::tensorstore::Index;
using ::tensorstore::TensorStore;

absl::Status Run(tensorstore::Spec &input_spec) {
    TensorStore<int8_t, tensorstore::dynamic_rank,
                tensorstore::ReadWriteMode::dynamic>
        store;
    TENSORSTORE_ASSIGN_OR_RETURN(
        store,
        tensorstore::Open<int8_t>(input_spec, tensorstore::OpenMode::open,
                                  tensorstore::ReadWriteMode::read)
            .result());

    tensorstore::IndexDomainView<tensorstore::dynamic_rank> domain =
        store.domain();
    tensorstore::span<const Index, tensorstore::dynamic_rank> shape =
        domain.shape();
    const auto variant_count = shape[0];
    const auto sample_count = shape[1];

    tensorstore::ChunkLayout chunk_layout;
    TENSORSTORE_ASSIGN_OR_RETURN(chunk_layout, store.chunk_layout());

    tensorstore::ChunkLayout::ReadChunkShape chunk_shape =
        chunk_layout.read_chunk_shape();
    const auto variant_chunk_size = std::min(chunk_shape[0], variant_count);
    const auto sample_chunk_size = std::min(chunk_shape[1], sample_count);
    std::cout << "Chunk size: " << variant_chunk_size << ", " << sample_chunk_size  << std::endl;
    std::vector<uint64_t> bin_counts(11);
    std::vector<int8_t> data_vector(variant_chunk_size * sample_chunk_size * 2);

    for (auto variant_chunk_start = 0; variant_chunk_start < variant_count;
         variant_chunk_start += variant_chunk_size) {
        const auto variant_chunk_end =
            std::min(variant_count, variant_chunk_start + variant_chunk_size);
        const auto variant_chunk_len = variant_chunk_end - variant_chunk_start;
        /* std::cout << "Chunk : " << variant_chunk_start << " len = " << variant_chunk_len << std::endl; */
        std::vector<uint64_t> ref_counts(variant_chunk_len);
        std::vector<uint64_t> het_counts(variant_chunk_len);
        std::vector<uint64_t> hom_alt_counts(variant_chunk_len);

        for (auto sample_chunk_start = 0; sample_chunk_start < sample_count;
             sample_chunk_start += sample_chunk_size) {
            const auto sample_chunk_end =
                std::min(sample_count, sample_chunk_start + sample_chunk_size);
            const auto sample_chunk_len = sample_chunk_end - sample_chunk_start;

            tensorstore::SharedArray<int8_t, 3L,
                                     tensorstore::ArrayOriginKind::zero>
                array = tensorstore::UnownedToShared(tensorstore::Array(
                    data_vector.data(),
                    {variant_chunk_size, sample_chunk_size, 2},
                    tensorstore::c_order));
            tensorstore::Read(
                store | tensorstore::Dims(0, 1).HalfOpenInterval(
                            {variant_chunk_start, sample_chunk_start},
                            {variant_chunk_end, sample_chunk_end}),
                array)
                .Wait();

            for (auto variant_index = 0; variant_index < variant_chunk_len; variant_index++) {
                for (auto sample_index = 0; sample_index < sample_chunk_len; sample_index++) {
                    const auto call_index =
                        2 * (sample_chunk_size * variant_index + sample_index);
                    const int8_t a = data_vector[call_index];
                    const int8_t b = data_vector[call_index + 1];

                    ref_counts[variant_index] += (a == 0) + (b == 0);
                    het_counts[variant_index] += a != b;
                    hom_alt_counts[variant_index] += a == b && a > 0;
                }
            }
        }

        for (auto variant_index = 0; variant_index < variant_chunk_len; variant_index++) {
            const uint64_t ref_count = ref_counts[variant_index];
            const uint64_t het_count = het_counts[variant_index];
            const uint64_t hom_alt_count = hom_alt_counts[variant_index];

            const uint64_t alt_count = 2 * sample_count - ref_count;
            const double alt_freq = alt_count / (2.0 * sample_count);
            const double het_ref_freq = 2 * alt_freq * (1 - alt_freq);
            const double hom_alt_freq = alt_freq * alt_freq;
            int bin_index = 10 * het_ref_freq;

            bin_counts[bin_index] += het_count;

            bin_index = 10 * hom_alt_freq;

            bin_counts[bin_index] += hom_alt_count;
        }
    }

    bin_counts[9] += bin_counts[10];
    bin_counts.pop_back();

    std::cout << "# PROB_DIST, genotype probability distribution, assumes HWE"
              << std::endl;

    for (auto bin_index = 0; bin_index < bin_counts.size(); bin_index++) {
        const double bin_start = bin_index / 10.0;
        const double bin_end = bin_start + 0.1;

        std::cout << "PROB_DIST\t" << bin_start << "\t" << bin_end << "\t"
                  << bin_counts[bin_index] << std::endl;
    }

    return absl::OkStatus();
}

int main(int argc, char *argv[]) {
    char *path = argv[1];
    tensorstore::Spec input_spec =
        tensorstore::Spec::FromJson(
            {
                {"driver", "zarr"},
                {"kvstore", {{"driver", "file"}, {"path", path}}},
            })
            .value();
    absl::Status status = Run(input_spec);

    if (!status.ok()) {
        std::cerr << status << std::endl;
    }

    return 0;
}
