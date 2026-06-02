#include <stdio.h>
#include "matrix.h"

int main() {
    printf("=== Library Operation Test ===\n\n");

    // 1. Create a 2x2 matrix and initialize values
    Matrix* A = mat_create(2, 2);
    if (A != NULL) {
        A->data[0][0] = 1.0; A->data[0][1] = 2.0;
        A->data[1][0] = 3.0; A->data[1][1] = 4.0;
    }

    // 2. Create a vector of size 2 and initialize values
    Vector* v = vec_create(2);
    if (v != NULL) {
        v->data[0] = 1.0;
        v->data[1] = 2.0;
    }

    // 3. Print input values
    printf("[Matrix A]\n");
    mat_print(A);

    printf("[Vector v]\n");
    vec_print(v);

    // 4. Matrix-Vector multiplication (A * v) test
    Vector* result = mat_vec_multiply(A, v);
    printf("[Result: A * v]\n");
    vec_print(result);

    // 5. Memory free (Crucial!)
    mat_free(A);
    vec_free(v);
    vec_free(result);

    printf("Test completed successfully!\n");
    return 0;
}