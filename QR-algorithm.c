#include "matrix.h"


/* Givens rotation: (a, b) -> (r, 0) 을 만드는 (c, s) 계산 */
static void givens_params(double a, double b, double *c, double *s) {
    if (b == 0.0) {
        *c = 1.0;
        *s = 0.0;
    } else {
        double r = sqrt(a * a + b * b);
        *c = a / r;
        *s = -b / r;
    }
}

/* H = G * H: 행 i, i+1에 왼쪽 Givens 적용 */
static void apply_givens_left(Matrix *H, int i, double c, double s, int col_start, int col_end) {
    for (int j = col_start; j <= col_end; j++) {
        double tmp   =  c * H->data[i][j] - s * H->data[i+1][j];
        H->data[i+1][j] =  s * H->data[i][j] + c * H->data[i+1][j];
        H->data[i][j]   = tmp;
    }
}

/* H = H * G^T: 열 i, i+1에 오른쪽 Givens 적용 */
static void apply_givens_right(Matrix *H, int i, double c, double s, int row_start, int row_end) {
    for (int k = row_start; k <= row_end; k++) {
        double tmp      =  c * H->data[k][i] - s * H->data[k][i+1];
        H->data[k][i+1] =  s * H->data[k][i] + c * H->data[k][i+1];
        H->data[k][i]   = tmp;
    }
}

/* Householder 반사로 A를 upper Hessenberg 형태로 축소 (in-place) */
static void reduce_to_hessenberg(Matrix *A) {
    int n = A->rows;

    for (int k = 0; k < n - 2; k++) {
        int sub_size = n - k - 1;
        Vector *x = vec_create(sub_size);
        for (int i = 0; i < sub_size; i++) {
            x->data[i] = A->data[k + 1 + i][k];
        }

        double x_norm = vec_norm(x);
        if (x_norm < 1e-15) {
            vec_free(x);
            continue;
        }

        Vector *v = vec_create(sub_size);
        for (int i = 0; i < sub_size; i++) {
            v->data[i] = x->data[i];
        }
        /* 부호 선택으로 수치 안정성 확보 */
        v->data[0] += (x->data[0] >= 0.0 ? 1.0 : -1.0) * x_norm;

        double v_norm = vec_norm(v);
        if (v_norm < 1e-15) {
            vec_free(x);
            vec_free(v);
            continue;
        }
        for (int i = 0; i < sub_size; i++) {
            v->data[i] /= v_norm;
        }

        /* 왼쪽 적용: A[k+1:, k:] -= 2v(v^T A[k+1:, k:]) */
        for (int j = k; j < n; j++) {
            double dot = 0.0;
            for (int i = 0; i < sub_size; i++) {
                dot += v->data[i] * A->data[k + 1 + i][j];
            }
            for (int i = 0; i < sub_size; i++) {
                A->data[k + 1 + i][j] -= 2.0 * v->data[i] * dot;
            }
        }

        /* 오른쪽 적용: A[:, k+1:] -= 2(A[:, k+1:]v)v^T */
        for (int i = 0; i < n; i++) {
            double dot = 0.0;
            for (int j = 0; j < sub_size; j++) {
                dot += A->data[i][k + 1 + j] * v->data[j];
            }
            for (int j = 0; j < sub_size; j++) {
                A->data[i][k + 1 + j] -= 2.0 * dot * v->data[j];
            }
        }

        vec_free(x);
        vec_free(v);
    }
}

/* 오른쪽 하단 2x2 블록의 eigenvalue 중 H[m][m]에 가장 가까운 값을 shift로 반환 */
static double wilkinson_shift(const Matrix *H, int m) {
    double a = H->data[m-1][m-1];
    double b = H->data[m-1][m];
    double c = H->data[m][m-1];
    double d = H->data[m][m];

    double trace  = a + d;
    double det    = a * d - b * c;
    double disc   = (trace * trace / 4.0) - det;

    double sqrt_disc;
    if (disc < 0.0) {
        /* 복소 eigenvalue 케이스 — Rayleigh quotient shift로 fallback */
        return d;
    }
    sqrt_disc = sqrt(disc);

    double lambda1 = trace / 2.0 + sqrt_disc;
    double lambda2 = trace / 2.0 - sqrt_disc;

    return (fabs(lambda1 - d) <= fabs(lambda2 - d)) ? lambda1 : lambda2;
}


Vector* qr_algorithm(const Matrix *A, int max_iter, double tol) {
    if (A == NULL || A->rows != A->cols) {
        fprintf(stderr, "[qr_algorithm] 오류: NULL이거나 정방 행렬이 아닙니다.\n");
        return NULL;
    }

    int n = A->rows;
    Vector *eigenvalues = vec_create(n);

    Matrix *H = mat_copy(A);
    reduce_to_hessenberg(H);

    int m = n - 1;

    for (int iter = 0; iter < max_iter && m > 0; iter++) {
        /* deflation: 수렴한 서브대각 원소를 0으로 확정 */
        while (m > 0 &&
               fabs(H->data[m][m-1]) < tol * (fabs(H->data[m-1][m-1]) + fabs(H->data[m][m]))) {
            H->data[m][m-1] = 0.0;
            eigenvalues->data[m] = H->data[m][m];
            m--;
        }
        if (m == 0) {
            eigenvalues->data[0] = H->data[0][0];
            break;
        }

        double shift = wilkinson_shift(H, m);

        /* implicit Q: 첫 번째 Givens로 bulge 생성 후 아래로 chase */
        double c, s;
        givens_params(H->data[0][0] - shift, H->data[1][0], &c, &s);
        apply_givens_left(H, 0, c, s, 0, m);
        apply_givens_right(H, 0, c, s, 0, (2 <= m ? 2 : m));

        for (int k = 1; k < m; k++) {
            givens_params(H->data[k][k-1], H->data[k+1][k-1], &c, &s);
            H->data[k][k-1]   = c * H->data[k][k-1] - s * H->data[k+1][k-1];
            H->data[k+1][k-1] = 0.0;

            int col_end = (k + 2 <= m) ? k + 2 : m;
            apply_givens_left(H, k, c, s, k, m);
            apply_givens_right(H, k, c, s, 0, col_end);
        }
    }

    if (m == 0) {
        eigenvalues->data[0] = H->data[0][0];
    } else if (m > 0) {
        fprintf(stderr, "[qr_algorithm] 경고: %d번 반복 후 완전 수렴하지 않았습니다."
                        " 대각 원소를 근사 eigenvalue로 반환합니다.\n", max_iter);
        for (int i = 0; i <= m; i++) {
            eigenvalues->data[i] = H->data[i][i];
        }
    }

    mat_free(H);
    return eigenvalues;
}
