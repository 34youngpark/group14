#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"

/* lu.c에 구현된 함수 원형 선언 */
int lu_decomposition_partial_pivoting(
    const Matrix* A,
    Matrix** L,
    Matrix** U,
    Matrix** P
);

int main(void)
{
    int n;

    printf("Enter matrix size (n x n): ");

    if (scanf("%d", &n) != 1 || n <= 0) {
        fprintf(stderr, "Error: Invalid matrix size.\n");
        return 1;
    }

    Matrix* A = mat_create(n, n);

    if (A == NULL) {
        fprintf(stderr, "Failed to create matrix A.\n");
        return 1;
    }

    printf("\nEnter elements of matrix A (%d x %d):\n", n, n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {

            printf("A[%d][%d] = ", i, j);

            if (scanf("%lf", &A->data[i][j]) != 1) {
                fprintf(stderr, "Error: Invalid input.\n");
                mat_free(A);
                return 1;
            }
        }
    }

    Matrix* L = NULL;
    Matrix* U = NULL;
    Matrix* P = NULL;

    if (lu_decomposition_partial_pivoting(A, &L, &U, &P) != 0) {
        fprintf(stderr, "LU decomposition failed.\n");
        mat_free(A);
        return 1;
    }

    printf("\n===== Matrix A =====\n");
    mat_print(A);

    printf("\n===== Matrix P =====\n");
    mat_print(P);

    printf("\n===== Matrix L =====\n");
    mat_print(L);

    printf("\n===== Matrix U =====\n");
    mat_print(U);

    Matrix* PA = mat_multiply(P, A);
    Matrix* LU = mat_multiply(L, U);

    if (PA != NULL && LU != NULL) {

        printf("\n===== P * A =====\n");
        mat_print(PA);

        printf("\n===== L * U =====\n");
        mat_print(LU);

    } else {
        fprintf(stderr, "Matrix multiplication failed.\n");
    }

    mat_free(PA);
    mat_free(LU);

    mat_free(A);
    mat_free(L);
    mat_free(U);
    mat_free(P);

    return 0;
}
