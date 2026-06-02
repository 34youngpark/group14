#include "matrix.h"

/* =========================================
   내부 상수
 * ========================================= */
#define MAX_ITER  10000
#define TOLERANCE 1e-10

/* =========================================
   static 내부 헬퍼: partial pivoting 가우스 소거법

   * (A - sigma*I) * z = q 형태의 선형계를 풀기 위해
     inverse_power_iteration 내부에서만 사용.
  ========================================= */
static Vector* solve_linear_system(const Matrix* A, const Vector* b) {
    int n = A->rows;
    if (A->cols != n || b->size != n) {
        fprintf(stderr, "[solve_linear_system] 오류: 차원 불일치\n");
        return NULL;
    }

    /* 확대 행렬 [A | b] 구성 */
    Matrix* aug = mat_create(n, n + 1);
    if (!aug) return NULL;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            aug->data[i][j] = A->data[i][j];
        aug->data[i][n] = b->data[i];
    }

    /* 전진 소거 (partial pivoting) */
    for (int col = 0; col < n; col++) {
        int pivot = col;
        double maxval = fabs(aug->data[col][col]);
        for (int row = col + 1; row < n; row++) {
            if (fabs(aug->data[row][col]) > maxval) {
                maxval = fabs(aug->data[row][col]);
                pivot = row;
            }
        }
        if (maxval < 1e-14) {
            fprintf(stderr, "[solve_linear_system] 오류: 특이 행렬(singular matrix)\n");
            mat_free(aug);
            return NULL;
        }

        if (pivot != col) {
            for (int j = 0; j <= n; j++) {
                double tmp          = aug->data[col][j];
                aug->data[col][j]   = aug->data[pivot][j];
                aug->data[pivot][j] = tmp;
            }
        }

        for (int row = col + 1; row < n; row++) {
            double factor = aug->data[row][col] / aug->data[col][col];
            for (int j = col; j <= n; j++)
                aug->data[row][j] -= factor * aug->data[col][j];
        }
    }

    /* 후진 대입 */
    Vector* x = vec_create(n);
    if (!x) { mat_free(aug); return NULL; }

    for (int i = n - 1; i >= 0; i--) {
        double sum = aug->data[i][n];
        for (int j = i + 1; j < n; j++)
            sum -= aug->data[i][j] * x->data[j];
        x->data[i] = sum / aug->data[i][i];
    }

    mat_free(aug);
    return x;
}

/* =========================================
 * Power Iteration (거듭제곱 반복법)
 
    행렬 A의 지배 고유값(절댓값 최대)과
    대응하는 고유벡터를 계산.
 
   Parameter:
      A          - 입력 정방 행렬 (n x n)
      max_iter   - 최대 반복 횟수 (0이면 기본값 사용)
      tol        - 수렴 허용 오차  (0.0이면 기본값 사용)
      eigenvalue - 출력: 지배 고유값
      iters_done - 출력: 수행된 반복 횟수 (NULL 가능)

    ∴ output: 고유벡터 (호출자가 vec_free로 해제), 오류 시 NULL
 ========================================= */
Vector* power_iteration(const Matrix* A, int max_iter, double tol,
                        double* eigenvalue, int* iters_done) {
    if (!A || A->rows != A->cols) {
        fprintf(stderr, "[power_iteration] 오류: NULL이거나 정방 행렬이 아닙니다.\n");
        return NULL;
    }
    if (!eigenvalue) {
        fprintf(stderr, "[power_iteration] 오류: eigenvalue 포인터가 NULL입니다.\n");
        return NULL;
    }

    int n = A->rows;
    if (max_iter <= 0) max_iter = MAX_ITER;
    if (tol <= 0.0)    tol      = TOLERANCE;

    /* 초기 벡터: [1, 0, ..., 0] */
    Vector* q = vec_create(n);
    if (!q) return NULL;
    q->data[0] = 1.0;

    double lambda = 0.0;
    int iter = 0;

    for (iter = 0; iter < max_iter; iter++) {
        /* z = A * q */
        Vector* z = mat_vec_multiply(A, q);
        if (!z) { vec_free(q); return NULL; }

        /* Rayleigh Quotient: lambda = q^T * A * q */
        double lambda_new = vec_dot_product(q, z);

        /* 정규화: q = z / ||z|| */
        double norm = vec_norm(z);
        if (norm < 1e-14) {
            fprintf(stderr, "[power_iteration] 경고: 영벡터 발생 (불변 부분공간?)\n");
            vec_free(z);
            vec_free(q);
            return NULL;
        }
        for (int i = 0; i < n; i++)
            z->data[i] /= norm;

        vec_free(q);
        q = z;

        /* 수렴 판정 */
        if (fabs(lambda_new - lambda) < tol) {
            lambda = lambda_new;
            iter++;
            break;
        }
        lambda = lambda_new;
    }

    *eigenvalue = lambda;
    if (iters_done) *iters_done = iter;
    return q;
}

/* =========================================
 * Inverse Power Iteration (역반복법)
 
    shift sigma에 가장 가까운 고유값과 대응하는 고유벡터를 계산.
    sigma = 0.0 으로 설정하면 절댓값 최소 고유값을 구합니다.
 
   Parameters:
     A          - 입력 정방 행렬 (n x n)
     sigma      - shift 값 (최솟값 탐색 시 0.0)
     max_iter   - 최대 반복 횟수 (0이면 기본값 사용)
     tol        - 수렴 허용 오차  (0.0이면 기본값 사용)
     eigenvalue - 출력: sigma에 가장 가까운 고유값
     iters_done - 출력: 수행된 반복 횟수 (NULL 가능)
 
  ∴ output : 고유벡터 (호출자가 vec_free로 해제), 오류 시 NULL
 * ========================================= */
Vector* inverse_power_iteration(const Matrix* A, double sigma,
                                int max_iter, double tol,
                                double* eigenvalue, int* iters_done) {
    if (!A || A->rows != A->cols) {
        fprintf(stderr, "[inverse_power_iteration] 오류: NULL이거나 정방 행렬이 아닙니다.\n");
        return NULL;
    }
    if (!eigenvalue) {
        fprintf(stderr, "[inverse_power_iteration] 오류: eigenvalue 포인터가 NULL입니다.\n");
        return NULL;
    }

    int n = A->rows;
    if (max_iter <= 0) max_iter = MAX_ITER;
    if (tol <= 0.0)    tol      = TOLERANCE;

    /* B = A - sigma * I (shift 적용) */
    Matrix* I = mat_identity(n);
    if (!I) return NULL;

    Matrix* B = mat_create(n, n);
    if (!B) { mat_free(I); return NULL; }

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            B->data[i][j] = A->data[i][j] - sigma * I->data[i][j];
    mat_free(I);

    /* 초기 벡터: [1, 0, ..., 0] */
    Vector* q = vec_create(n);
    if (!q) { mat_free(B); return NULL; }
    q->data[0] = 1.0;

    double lambda = 0.0;
    int iter = 0;

    for (iter = 0; iter < max_iter; iter++) {
        /* (A - sigma*I) * z = q 를 가우스 소거법으로 풀기 */
        Vector* z = solve_linear_system(B, q);
        if (!z) { vec_free(q); mat_free(B); return NULL; }

        /* Rayleigh Quotient: mu = q^T * z */
        double mu = vec_dot_product(q, z);

        /* 정규화: q = z / ||z|| */
        double norm = vec_norm(z);
        if (norm < 1e-14) {
            fprintf(stderr, "[inverse_power_iteration] 경고: 영벡터 발생\n");
            vec_free(z);
            vec_free(q);
            mat_free(B);
            return NULL;
        }
        for (int i = 0; i < n; i++)
            z->data[i] /= norm;

        vec_free(q);
        q = z;

        /* A의 고유값 복원: lambda = sigma + 1/mu */
        double lambda_new = sigma + 1.0 / mu;

        /* 수렴 판정 */
        if (fabs(lambda_new - lambda) < tol) {
            lambda = lambda_new;
            iter++;
            break;
        }
        lambda = lambda_new;
    }

    mat_free(B);

    *eigenvalue = lambda;
    if (iters_done) *iters_done = iter;
    return q;
}
