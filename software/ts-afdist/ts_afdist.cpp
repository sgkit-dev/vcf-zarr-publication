
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

absl::Status Run(tensorstore::Spec input_spec) {
    TensorStore<int8_t, tensorstore::dynamic_rank, tensorstore::ReadWriteMode::dynamic> store;
    TENSORSTORE_ASSIGN_OR_RETURN(
        store,
        tensorstore::Open<int8_t>(input_spec, tensorstore::OpenMode::open, tensorstore::ReadWriteMode::read).result());

    tensorstore::IndexDomainView<tensorstore::dynamic_rank> domain = store.domain();
    tensorstore::span<const Index, tensorstore::dynamic_rank> shape = domain.shape();
    const Index variant_count = shape[0];
    const Index sample_count = shape[1];

    tensorstore::ChunkLayout chunk_layout;
    TENSORSTORE_ASSIGN_OR_RETURN(
        chunk_layout, store.chunk_layout());

    tensorstore::ChunkLayout::ReadChunkShape chunk_shape = chunk_layout.read_chunk_shape();
    const Index variant_chunk_size = chunk_shape[0];
    const Index sample_chunk_size = chunk_shape[1];
    std::vector<uint64_t> bin_counts(11);

    for (Index variant_chunk_start = 0; variant_chunk_start < variant_count; variant_chunk_start += variant_chunk_size) {
        for (Index sample_chunk_start = 0; sample_chunk_start < sample_count; sample_chunk_start += sample_chunk_size) {
            const Index variant_chunk_end = std::min(variant_count, variant_chunk_start + variant_chunk_size);
            const Index sample_chunk_end = std::min(sample_count, sample_chunk_start + sample_chunk_size);
            tensorstore::SharedArray<int8_t, tensorstore::dynamic_rank, tensorstore::offset_origin> array;

            TENSORSTORE_ASSIGN_OR_RETURN(array, tensorstore::Read(store | tensorstore::Dims(0, 1).HalfOpenInterval(
                                                                              {variant_chunk_start, sample_chunk_start},
                                                                              {variant_chunk_end, sample_chunk_end}))
                                                    .result());

            for (Index variant_index = variant_chunk_start; variant_index < variant_chunk_end; variant_index++) {
                uint64_t hom_ref_count = 0;
                uint64_t hom_alt_count = 0;
                uint64_t het_count = 0;
                uint64_t ref_count = 0;

                for (Index sample_index = sample_chunk_start; sample_index < sample_chunk_end; sample_index++) {
                    const int8_t a = array(variant_index, sample_index, 0);
                    const int8_t b = array(variant_index, sample_index, 1);

                    if (std::min(a, b) < 0) continue;

                    if (a == b)
                        if (a == 0)
                            hom_ref_count += 1;
                        else
                            hom_alt_count += 1;
                    else
                        het_count += 1;

                    if (a == 0) ref_count += 1;
                    if (b == 0) ref_count += 1;
                }

                const uint64_t alt_count = 2 * sample_count - ref_count;
                const double alt_freq = alt_count / (2.0 * sample_count);
                const double het_ref_freq = 2 * alt_freq * (1 - alt_freq);
                const double hom_alt_freq = alt_freq * alt_freq;
                Index bin_index = 10 * het_ref_freq;

                bin_counts[bin_index] += het_count;

                bin_index = 10 * hom_alt_freq;

                bin_counts[bin_index] += hom_alt_count;
            }
        }
    }

    bin_counts[9] += bin_counts[10];
    bin_counts.pop_back();

    std::cout << "# PROB_DIST, genotype probability distribution, assumes HWE" << std::endl;

    for (Index bin_index = 0; bin_index < bin_counts.size(); bin_index++) {
        const double bin_start = bin_index / 10.0;
        const double bin_end = bin_start + 0.1;

        std::cout << "PROB_DIST\t" << bin_start << "\t" << bin_end << "\t" << bin_counts[bin_index] << std::endl;
    }

    return absl::OkStatus();
}

int main(int argc, char *argv[]) {
    char *path = argv[1];
    tensorstore::Spec input_spec = tensorstore::Spec::FromJson({
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
