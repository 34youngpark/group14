#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* =========================================
 * 1. 자료구조 정의
 * ========================================= */
typedef struct {
    int rows;
    int cols;
    double *block; // 연속 메모리 블록
    double **data; // 2차원 접근용 포인터 배열
} Matrix;

typedef struct {
    int size;
    double *data;
} Vector;

/* =========================================
 * 2. 메모리 할당 및 해제
 * ========================================= */
Matrix* mat_create(int rows, int cols);
Vector* vec_create(int size);
void mat_free(Matrix* m);
void vec_free(Vector* v);

/* =========================================
 * 3. 행렬 기본 연산
 * ========================================= */
Matrix* mat_add(const Matrix* a, const Matrix* b);
Matrix* mat_sub(const Matrix* a, const Matrix* b);
Matrix* mat_multiply(const Matrix* a, const Matrix* b);
Matrix* mat_transpose(const Matrix* m);

/* =========================================
 * 4. 벡터 및 혼합 연산
 * ========================================= */
Vector* mat_vec_multiply(const Matrix* m, const Vector* v);
double vec_dot_product(const Vector* v1, const Vector* v2);
Matrix* vec_outer_product(const Vector* v1, const Vector* v2);
double vec_norm(const Vector* v);

/* =========================================
 * 5. 유틸리티 기능
 * ========================================= */
Matrix* mat_identity(int size);
Matrix* mat_copy(const Matrix* m);
void mat_swap_rows(Matrix* m, int row1, int row2);
void mat_print(const Matrix* m);
void vec_print(const Vector* v);

#endif
