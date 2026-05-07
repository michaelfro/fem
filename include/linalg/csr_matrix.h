#pragma once

#include <vector>
#include <cstddef>

struct csr_matrix {
    std::vector<double> values;
    std::vector<size_t> col_index;
    std::vector<size_t> row_index;
    int size; // do I need this?
};
