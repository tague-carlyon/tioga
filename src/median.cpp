// Copyright TIOGA Developers. See COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD 3-Clause)
//
// Fast median finder implementation
// Replaces median.F90 with modern C++ implementation
// Provides identical function signature for backward compatibility

#include <algorithm>

extern "C" {

// Fast median finding using simple algorithm
// Focus on correctness over performance initially
void median_(int *ix, double *x, int *n, double *xmed)
{
    // Handle edge cases
    if (*n <= 0) {
        *xmed = 0.0;
        return;
    }
    
    if (*n == 1) {
        *xmed = x[0];
        return;
    }
    
    if (*n == 2) {
        *xmed = 0.5 * (x[0] + x[1]);
        if (x[1] < x[0]) {
            double temp = x[0];
            x[0] = x[1];
            x[1] = temp;
            
            int itemp = ix[0];
            ix[0] = ix[1];
            ix[1] = itemp;
        }
        return;
    }
    
    // For n >= 3, use sorting approach initially
    // Sort the values while keeping track of original indices
    struct IndexedValue {
        double value;
        int original_index;
    };
    
    // Create indexed array
    IndexedValue* indexed = new IndexedValue[*n];
    for (int i = 0; i < *n; i++) {
        indexed[i].value = x[i];
        indexed[i].original_index = ix[i];
    }
    
    // Sort by value
    std::sort(indexed, indexed + *n, [](const IndexedValue& a, const IndexedValue& b) {
        return a.value < b.value;
    });
    
    // Find median position
    int median_pos = *n / 2;
    *xmed = indexed[median_pos].value;
    
    // Update ix array to reflect new order (mimic Fortran behavior)
    for (int i = 0; i < *n; i++) {
        ix[i] = indexed[i].original_index;
    }
    
    delete[] indexed;
}

} // extern "C"