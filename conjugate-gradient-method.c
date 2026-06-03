#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// 행렬 구조체 정의
typedef struct {
    int rows;
    int cols;
    double** data;
} Matrix;

// 벡터 구조체 정의
typedef struct {
    int len;
    double* data;
} Vector;

// 행렬 생성 함수
Matrix create_matrix(int rows, int cols) {
    Matrix mat;
    mat.rows = rows;
    mat.cols = cols;
    mat.data = (double**)malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; i++) {
        mat.data[i] = (double*)calloc(cols, sizeof(double));
    }
    return mat;
}

// 벡터 생성 함수
Vector create_vector(int len) {
    Vector vec;
    vec.len = len;
    vec.data = (double*)calloc(len, sizeof(double));
    return vec;
}

// 벡터를 0으로 초기화
void vector_set_zero(Vector vec) {
    for (int i = 0; i < vec.len; i++) vec.data[i] = 0;
}

// 벡터 복사
void vector_copy(Vector src, Vector dest) {
    for (int i = 0; i < src.len; i++) dest.data[i] = src.data[i];
}

// 벡터 뺄셈: c = a - b
void vector_sub(Vector a, Vector b, Vector c) {
    for (int i = 0; i < a.len; i++) c.data[i] = a.data[i] - b.data[i];
}

// 벡터 덧셈: c = a + b
void vector_add(Vector a, Vector b, Vector c) {
    for (int i = 0; i < a.len; i++) c.data[i] = a.data[i] + b.data[i];
}

// 벡터 스칼라 곱: res = vec * s
void vector_scale(Vector vec, double s, Vector res) {
    for (int i = 0; i < vec.len; i++) res.data[i] = vec.data[i] * s;
}

// 스칼라 곱을 포함한 벡터 덧셈: res = a + s*b
void vector_add_scaled(Vector a, double s, Vector b, Vector res) {
    for (int i = 0; i < a.len; i++) res.data[i] = a.data[i] + s * b.data[i];
}

// 벡터 내적 계산
double vector_dot(Vector a, Vector b) {
    double sum = 0;
    for (int i = 0; i < a.len; i++) sum += a.data[i] * b.data[i];
    return sum;
}

// 행렬 × 벡터 곱셈: res = A * v
void matrix_vector_mult(Matrix A, Vector v, Vector res) {
    for (int i = 0; i < A.rows; i++) {
        res.data[i] = 0;
        for (int j = 0; j < A.cols; j++) {
            res.data[i] += A.data[i][j] * v.data[j];
        }
    }
}

// 메모리 해제 함수들
void free_matrix(Matrix m) {
    for (int i = 0; i < m.rows; i++) free(m.data[i]);
    free(m.data);
}

void free_vector(Vector v) {
    free(v.data);
}

// ===================== 공액 구배법(CG) 알고리즘 =====================
void conjugate_gradient(Matrix A, Vector b, Vector x, int n, int max_iter, double tol) {
    Vector r = create_vector(n);  // 잔차 벡터
    Vector p = create_vector(n);  // 탐색 방향 벡터
    Vector Ap = create_vector(n); // 임시 벡터 A*p

    vector_set_zero(x);
    matrix_vector_mult(A, x, Ap);
    vector_sub(b, Ap, r);
    vector_copy(r, p);

    double r_old = vector_dot(r, r);

    for (int k = 0; k < max_iter; k++) {
        matrix_vector_mult(A, p, Ap);
        double alpha = r_old / vector_dot(p, Ap);

        vector_add_scaled(x, alpha, p, x);
        vector_add_scaled(r, -alpha, Ap, r);

        double r_new = vector_dot(r, r);

        // 수렴 조건 판단
        if (r_new < tol * tol) {
            printf("수렴 성공! 반복 횟수: %d\n", k + 1);
            break;
        }

        double beta = r_new / r_old;
        vector_scale(p, beta, p);
        vector_add(r, p, p);
        r_old = r_new;
    }

    // 사용한 메모리 해제
    free_vector(r);
    free_vector(p);
    free_vector(Ap);
}

// ===================== 메인 함수 =====================
int main() {
    int n = 3;

    // 대칭 양정치 행렬 A 생성
    Matrix A = create_matrix(n, n);
    A.data[0][0] = 4; A.data[0][1] = 1; A.data[0][2] = 1;
    A.data[1][0] = 1; A.data[1][1] = 5; A.data[1][2] = 2;
    A.data[2][0] = 1; A.data[2][1] = 2; A.data[2][2] = 6;

    // 우변 벡터 b 생성
    Vector b = create_vector(n);
    b.data[0] = 6;
    b.data[1] = 8;
    b.data[2] = 9;

    // 해 벡터 x
    Vector x = create_vector(n);

    // CG 알고리즘 실행
    conjugate_gradient(A, b, x, n, 1000, 1e-6);

    // 결과 출력
    printf("\n계산 결과:\n");
    for (int i = 0; i < n; i++) {
        printf("x[%d] = %.6f\n", i, x.data[i]);
    }

    // 최종 메모리 해제
    free_matrix(A);
    free_vector(b);
    free_vector(x);
    return 0;
}
