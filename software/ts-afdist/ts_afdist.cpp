
#include <iostream>
#include "absl/status/status.h"
#include "tensorstore/array.h"
#include "tensorstore/open.h"
#include "tensorstore/spec.h"

/* #include "tensorstore/index.h" */
#include "tensorstore/index_space/dim_expression.h"
/* #include "tensorstore/index_space/transformed_array.h" */
/* #include "tensorstore/util/iterate_over_index_range.h" */
/* #include "tensorstore/util/status.h" */

using ::tensorstore::Context;
using ::tensorstore::Index;

/* void count_alleles(tensorstore::Array array) */
/* { */

/* } */

absl::Status Run(tensorstore::Spec input_spec, std::string output_filename) {

    auto context = Context::Default();
    TENSORSTORE_ASSIGN_OR_RETURN(
        auto input,
        tensorstore::Open(input_spec, context, tensorstore::OpenMode::open,
                          tensorstore::ReadWriteMode::read)
            .result());

    auto shape = input.domain().shape();
    auto num_variants = shape[0];
    auto num_samples = shape[1];
    auto ploidy = shape[2];
    auto chunk_shape = input.chunk_layout().value().read_chunk_shape();
    auto v_chunk_size = chunk_shape[0];
    auto s_chunk_size = chunk_shape[1];

    /* if (ploidy != 2) { */
    /*     std::cerr << "Ploidy must be 2"; */
    /*     return */

    /* } */
    std::cout <<  typeid(chunk_shape).name() << std::endl;
    std::cout << "input shape = " << shape << " chunk shape = " << chunk_shape << std::endl;

    for (auto v_chunk_offset = 0; v_chunk_offset < num_variants;
            v_chunk_offset += v_chunk_size) {

        for (auto s_chunk_offset = 0; s_chunk_offset < num_samples;
                s_chunk_offset += s_chunk_size) {
            std::cout << "Get chunk at " << v_chunk_offset << ", " << s_chunk_offset
                << std::endl;
            auto slice = tensorstore::Dims(0, 1).HalfOpenInterval(
                        {v_chunk_offset, v_chunk_offset + v_chunk_size},
                        {s_chunk_offset, s_chunk_offset + s_chunk_size});
            /* auto x = tensorstore::Read<tensorstore::zero_origin>(input | slice).result(); */

            std::vector<int8_t> vec(v_chunk_size * s_chunk_size * ploidy);
            auto arr = tensorstore::Array(
                vec.data(), {v_chunk_size, s_chunk_size, ploidy},
                tensorstore::c_order);
            /* tensorstore::Read(input , tensorstore::UnownedToShared(arr)).value(); */
                    /* std::cout << "x = " << slice << std::endl; */
            /* count_alleles(arr); */
            // FIXME breaking here. Trying to get tensorstore to read this slice
            // into the std vector memory allocated above.
            auto z = tensorstore::Read(input |
                    tensorstore::Dims(0, 1).HalfOpenInterval(
                        {v_chunk_offset, v_chunk_offset + v_chunk_size},
                        {s_chunk_offset, s_chunk_offset + s_chunk_size})
                    ).result();
            /* std::cout << "Got " << z << std::endl; */

        }
    }


    /* return absl::OkStatus(); */

    /* TENSORSTORE_ASSIGN_OR_RETURN( */
    /*     auto data, */
    /*    tensorstore::Read(input).result()); */

    /* /1* const auto max = data.shape()[data.rank() - 1] - 1; *1/ */
    /* auto dtype = data.dtype(); */
    /* std::cout << "dtype = " << dtype << std::endl; */

    /* /1* Read the data into a std::vector of the right size, based on */
    /*  * https://github.com/google/tensorstore/issues/36#issuecomment-1292757434 */
    /*  *1/ */
    /* std::vector<int8_t> vec(num_variants * num_samples * ploidy); */
    /* auto arr = tensorstore::Array(vec.data(), {num_variants, num_samples, ploidy}, */
    /*         tensorstore::c_order); */
    /* tensorstore::Read(input, tensorstore::UnownedToShared(arr)).value(); */


    /* auto offset = 0; */
    /* for (auto v = 0; v < num_variants; v++) { */
    /*     for (auto s = 0; s < num_samples; s++) { */
    /*         for (auto z = 0; z < ploidy; z++) { */
    /*             /1* For some mysterious reason std::cout won't show these values, */
    /*              * but printf works fine?? *1/ */
    /*             /1* std::cout << offset << " value = " << p[offset] << std::endl; *1/ */
    /*             printf("vec[%d] = %d\n", (int) offset, (int) vec[offset]); */
    /*             offset ++; */
    /*         } */
    /*     } */
    /* } */

    return absl::OkStatus();
}

int main(int argc, char* argv[]) {
    auto path = argv[1];
    tensorstore::Spec input_spec =
        tensorstore::Spec::FromJson(
            {
                {"driver", "zarr"}, {"kvstore", {{"driver", "file"}, {"path", path}}},
                /* {"path", "input"}, */
            })
            .value();
    auto status = Run(input_spec, "tmp.txt");
    if (!status.ok()) {
        std::cerr << status << std::endl;
    }
    return 0;
}
