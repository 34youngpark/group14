#include "matrix.h"

/* =========================================
 * 메모리 할당 및 해제
 * ========================================= */
Matrix* mat_create(int rows, int cols) {
    if (rows <= 0 || cols <= 0) {
        fprintf(stderr, "[Error] mat_create: Invalid dimensions (%d x %d).\n", rows, cols);
        return NULL;
    }
    Matrix* m = (Matrix*)malloc(sizeof(Matrix));
    if (!m) return NULL;
    m->rows = rows;
    m->cols = cols;
    m->block = (double*)calloc(rows * cols, sizeof(double));
    if (!m->block) {
        free(m);
        return NULL;
    }
    m->data = (double**)malloc(rows * sizeof(double*));
    if (!m->data) {
        free(m->block);
        free(m);
        return NULL;
    }
    for (int i = 0; i < rows; i++) {
        m->data[i] = m->block + (i * cols);
    }
    return m;
}

Vector* vec_create(int size) {
    if (size <= 0) {
        fprintf(stderr, "[Error] vec_create: Invalid size (%d).\n", size);
        return NULL;
    }
    Vector* v = (Vector*)malloc(sizeof(Vector));
    if (!v) return NULL;
    v->size = size;
    v->data = (double*)calloc(size, sizeof(double));
    if (!v->data) {
        free(v);
        return NULL;
    }
    return v;
}

void mat_free(Matrix* m) {
    if (m) {
        if (m->data) free(m->data);
        if (m->block) free(m->block);
        free(m);
    }
}

void vec_free(Vector* v) {
    if (v) {
        if (v->data) free(v->data);
        free(v);
    }
}

/* =========================================
 * 행렬 기본 연산
 * ========================================= */
Matrix* mat_add(const Matrix* a, const Matrix* b) {
    if (!a || !b || a->rows != b->rows || a->cols != b->cols) {
        fprintf(stderr, "[Error] mat_add: Dimension mismatch.\n");
        return NULL;
    }
    Matrix* result = mat_create(a->rows, a->cols);
    for (int i = 0; i < a->rows * a->cols; i++) result->block[i] = a->block[i] + b->block[i];
    return result;
}

Matrix* mat_sub(const Matrix* a, const Matrix* b) {
    if (!a || !b || a->rows != b->rows || a->cols != b->cols) {
        fprintf(stderr, "[Error] mat_sub: Dimension mismatch.\n");
        return NULL;
    }
    Matrix* result = mat_create(a->rows, a->cols);
    for (int i = 0; i < a->rows * a->cols; i++) result->block[i] = a->block[i] - b->block[i];
    return result;
}

Matrix* mat_multiply(const Matrix* a, const Matrix* b) {
    if (!a || !b || a->cols != b->rows) {
        fprintf(stderr, "[Error] mat_multiply: Dimension mismatch.\n");
        return NULL;
    }
    Matrix* result = mat_create(a->rows, b->cols);
    for (int i = 0; i < a->rows; i++) {
        for (int k = 0; k < a->cols; k++) {
            double temp = a->data[i][k];
            for (int j = 0; j < b->cols; j++) result->data[i][j] += temp * b->data[k][j];
        }
    }
    return result;
}

Matrix* mat_transpose(const Matrix* m) {
    if (!m) return NULL;
    Matrix* result = mat_create(m->cols, m->rows);
    for (int i = 0; i < m->rows; i++) {
        for (int j = 0; j < m->cols; j++) result->data[j][i] = m->data[i][j];
    }
    return result;
}

/* =========================================
 * 벡터 및 혼합 연산
 * ========================================= */
Vector* mat_vec_multiply(const Matrix* m, const Vector* v) {
    if (!m || !v || m->cols != v->size) {
        fprintf(stderr, "[Error] mat_vec_multiply: Dimension mismatch.\n");
        return NULL;
    }
    Vector* result = vec_create(m->rows);
    for (int i = 0; i < m->rows; i++) {
        double sum = 0.0;
        for (int j = 0; j < m->cols; j++) sum += m->data[i][j] * v->data[j];
        result->data[i] = sum;
    }
    return result;
}

double vec_dot_product(const Vector* v1, const Vector* v2) {
    if (!v1 || !v2 || v1->size != v2->size) {
        fprintf(stderr, "[Error] vec_dot_product: Dimension mismatch.\n");
        return 0.0;
    }
    double sum = 0.0;
    for (int i = 0; i < v1->size; i++) sum += v1->data[i] * v2->data[i];
    return sum;
}

Matrix* vec_outer_product(const Vector* v1, const Vector* v2) {
    if (!v1 || !v2) return NULL;
    Matrix* result = mat_create(v1->size, v2->size);
    for (int i = 0; i < v1->size; i++) {
        for (int j = 0; j < v2->size; j++) result->data[i][j] = v1->data[i] * v2->data[j];
    }
    return result;
}

double vec_norm(const Vector* v) {
    if (!v) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < v->size; i++) sum += v->data[i] * v->data[i];
    return sqrt(sum);
}

/* =========================================
 * 유틸리티 기능
 * ========================================= */
Matrix* mat_identity(int size) {
    Matrix* result = mat_create(size, size);
    if (!result) return NULL;
    for (int i = 0; i < size; i++) result->data[i][i] = 1.0;
    return result;
}

Matrix* mat_copy(const Matrix* m) {
    if (!m) return NULL;
    Matrix* result = mat_create(m->rows, m->cols);
    for (int i = 0; i < m->rows * m->cols; i++) result->block[i] = m->block[i];
    return result;
}

void mat_swap_rows(Matrix* m, int row1, int row2) {
    if (!m || row1 < 0 || row1 >= m->rows || row2 < 0 || row2 >= m->rows || row1 == row2) return;
    double* temp = m->data[row1];
    m->data[row1] = m->data[row2];
    m->data[row2] = temp;
}

void mat_print(const Matrix* m) {
    if (!m) return;
    for (int i = 0; i < m->rows; i++) {
        for (int j = 0; j < m->cols; j++) printf("%9.4f ", m->data[i][j]);
        printf("\n");
    }
    printf("\n");
}

void vec_print(const Vector* v) {
    if (!v) return;
    for (int i = 0; i < v->size; i++) printf("%9.4f\n", v->data[i]);
    printf("\n");
}
