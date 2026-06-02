#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"

#define EPSILON 1e-12

/*
 * LU Decomposition with Partial Pivoting
 *
 * PA = LU
 *
 * 입력:
 *   A : 정사각 행렬
 *
 * 출력:
 *   L : Lower triangular matrix
 *   U : Upper triangular matrix
 *   P : Permutation matrix
 *
 * 반환값:
 *   성공 시 0
 *   실패 시 -1
 *
 * 주의:
 *   - 입력 행렬 A는 변경되지 않음
 *   - L, U, P는 함수 내부에서 새로 생성됨
 *   - 호출 측에서 반드시 mat_free() 필요
 */

static int is_square_matrix(const Matrix* m) {
    return (m != NULL && m->rows == m->cols);
}

static int find_pivot_row(const Matrix* m, int col, int start_row) {
    int pivot_row = start_row;
    double max_value = fabs(m->data[start_row][col]);

    for (int i = start_row + 1; i < m->rows; i++) {
        double current = fabs(m->data[i][col]);

        if (current > max_value) {
            max_value = current;
            pivot_row = i;
        }
    }

    return pivot_row;
}

int lu_decomposition_partial_pivoting(
    const Matrix* A,
    Matrix** L,
    Matrix** U,
    Matrix** P
) {
    if (A == NULL || L == NULL || U == NULL || P == NULL) {
        fprintf(stderr, "Error: NULL pointer input.\n");
        return -1;
    }

    if (!is_square_matrix(A)) {
        fprintf(stderr, "Error: LU decomposition requires square matrix.\n");
        return -1;
    }

    int n = A->rows;

    *L = mat_identity(n);
    *U = mat_copy(A);
    *P = mat_identity(n);

    if (*L == NULL || *U == NULL || *P == NULL) {
        fprintf(stderr, "Error: Memory allocation failed.\n");

        if (*L) mat_free(*L);
        if (*U) mat_free(*U);
        if (*P) mat_free(*P);

        return -1;
    }

    for (int k = 0; k < n - 1; k++) {

        // Partial Pivoting
        int pivot_row = find_pivot_row(*U, k, k);

        if (fabs((*U)->data[pivot_row][k]) < EPSILON) {
            fprintf(stderr, "Error: Singular matrix detected.\n");

            mat_free(*L);
            mat_free(*U);
            mat_free(*P);

            *L = NULL;
            *U = NULL;
            *P = NULL;

            return -1;
        }

        // Swap rows in U and P
        if (pivot_row != k) {
            mat_swap_rows(*U, k, pivot_row);
            mat_swap_rows(*P, k, pivot_row);

            // L의 이전 column 값들도 같이 swap
            for (int j = 0; j < k; j++) {
                double temp = (*L)->data[k][j];
                (*L)->data[k][j] = (*L)->data[pivot_row][j];
                (*L)->data[pivot_row][j] = temp;
            }
        }

        // Gaussian Elimination
        for (int i = k + 1; i < n; i++) {

            double factor = (*U)->data[i][k] / (*U)->data[k][k];

            (*L)->data[i][k] = factor;

            for (int j = k; j < n; j++) {
                (*U)->data[i][j] -= factor * (*U)->data[k][j];
            }
        }
    }

    // 마지막 pivot 검사
    if (fabs((*U)->data[n - 1][n - 1]) < EPSILON) {
        fprintf(stderr, "Error: Singular matrix detected.\n");

        mat_free(*L);
        mat_free(*U);
        mat_free(*P);

        *L = NULL;
        *U = NULL;
        *P = NULL;

        return -1;
    }

    return 0;
}

/* test main*/

int main(void)
{
    Matrix* A = mat_create(3, 3);

    if (A == NULL) {
        fprintf(stderr, "Failed to create matrix A.\n");
        return 1;
    }

    /*
        테스트 행렬

        A =
        [ 2   1   1 ]
        [ 4  -6   0 ]
        [ -2  7   2 ]
    */

    A->data[0][0] = 2.0;
    A->data[0][1] = 1.0;
    A->data[0][2] = 1.0;

    A->data[1][0] = 4.0;
    A->data[1][1] = -6.0;
    A->data[1][2] = 0.0;

    A->data[2][0] = -2.0;
    A->data[2][1] = 7.0;
    A->data[2][2] = 2.0;

    Matrix* L = NULL;
    Matrix* U = NULL;
    Matrix* P = NULL;

    if (lu_decomposition_partial_pivoting(A, &L, &U, &P) != 0) {
        fprintf(stderr, "LU decomposition failed.\n");
        mat_free(A);
        return 1;
    }

    printf("===== Matrix A =====\n");
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
