#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

/* 엔진 파일(eigen.c)에 있는 함수를 불러옵니다. */
extern Vector* power_iteration(const Matrix* A, int max_iter, double tol,
                               double* eigenvalue, int* iters_done);
extern Vector* inverse_power_iteration(const Matrix* A, double sigma,
                                       int max_iter, double tol,
                                       double* eigenvalue, int* iters_done);

int main(void) {
    int n;

    printf("==================================================\n");
    printf("   사용자 입력형 Power / Inverse Power Iteration\n");
    printf("==================================================\n\n");

    /* 0. 사용 안내 (입력 전, 사용자 친화 가이드) */
    printf("[ 사용 안내 ]\n");
    printf("이 프로그램은 정방 행렬의 고유값(eigenvalue)과 고유벡터를 구합니다.\n\n");
    printf("진행 순서:\n");
    printf("  1) 행렬 크기 N 을 입력합니다.        (예: 2)\n");
    printf("  2) 행렬 성분을 좌→우, 위→아래 순으로 입력합니다.\n");
    printf("  3) 사용할 알고리즘을 선택합니다.\n");
    printf("       - Power Iteration         : 절댓값이 가장 큰 고유값\n");
    printf("       - Inverse Power Iteration : sigma에 가장 가까운 고유값\n");
    printf("         (sigma=0 이면 절댓값이 가장 작은 고유값)\n\n");
    printf("예시) A = [ 4  1 ]   →  지배 고유값 5, 최소 고유값 2\n");
    printf("          [ 2  3 ]\n");
    printf("  N=2, 성분: 4 1 2 3, 알고리즘: 1 (Power) 선택 시  λ = 5\n");
    printf("--------------------------------------------------\n\n");

    /* 1. 행렬 크기 입력 */
    printf("행렬의 크기(N x N)에서 N을 입력하세요 (예: 3): ");
    if (scanf("%d", &n) != 1 || n <= 0) {
        fprintf(stderr, "올바른 양의 정수를 입력해주세요.\n");
        return 1;
    }

    /* 2. 행렬 생성 */
    Matrix* A = mat_create(n, n);
    if (!A) {
        fprintf(stderr, "행렬 생성에 실패했습니다.\n");
        return 1;
    }

    /* 3. 행렬 성분 입력 */
    printf("\n--- 행렬의 각 성분을 차례대로 입력하세요 ---\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("A[%d][%d] 성분 입력: ", i, j);
            if (scanf("%lf", &A->data[i][j]) != 1) {
                fprintf(stderr, "숫자 입력에 오류가 발생했습니다.\n");
                mat_free(A);
                return 1;
            }
        }
    }

    /* 4. 입력 행렬 확인 출력 */
    printf("\n[확인] 입력된 행렬 A:\n");
    for (int i = 0; i < n; i++) {
        printf("    [ ");
        for (int j = 0; j < n; j++)
            printf("%8.4f ", A->data[i][j]);
        printf("]\n");
    }

    /* 5. 알고리즘 선택 */
    int choice;
    double sigma = 0.0;
    printf("\n--------------------------------------------------\n");
    printf("실행할 알고리즘을 선택하세요:\n");
    printf("  1. Power Iteration         (지배 고유값)\n");
    printf("  2. Inverse Power Iteration (sigma에 가장 가까운 고유값)\n");
    printf("선택 (1 or 2): ");
    if (scanf("%d", &choice) != 1 || (choice != 1 && choice != 2)) {
        fprintf(stderr, "1 또는 2를 입력해주세요.\n");
        mat_free(A);
        return 1;
    }

    if (choice == 2) {
        printf("sigma 값을 입력하세요 (절댓값 최소 고유값 탐색 시 0): ");
        if (scanf("%lf", &sigma) != 1) {
            fprintf(stderr, "숫자 입력에 오류가 발생했습니다.\n");
            mat_free(A);
            return 1;
        }
    }

    /* 6. 계산 수행 */
    double eigenvalue = 0.0;
    int    iters      = 0;
    Vector* evec      = NULL;

    printf("\n[진행] 고유값 계산 중 (max_iter: %d, tol: %.0e)...\n",
           10000, 1e-10);

    if (choice == 1) {
        evec = power_iteration(A, 0, 0.0, &eigenvalue, &iters);
    } else {
        evec = inverse_power_iteration(A, sigma, 0, 0.0, &eigenvalue, &iters);
    }

    /* 7. 결과 출력 */
    printf("==================================================\n");
    if (evec) {
        if (choice == 1)
            printf("  [Power Iteration] 지배 고유값\n");
        else
            printf("  [Inverse Power Iteration] sigma = %.4f 에 가장 가까운 고유값\n", sigma);

        printf("  λ = %11.8f  (반복 횟수: %d)\n\n", eigenvalue, iters);
        printf("  대응 고유벡터:\n");
        for (int i = 0; i < n; i++)
            printf("    v[%d] = %11.8f\n", i, evec->data[i]);
        vec_free(evec);
    } else {
        fprintf(stderr, "  오류: 고유값 계산에 실패했습니다.\n");
    }
    printf("==================================================\n");

    mat_free(A);
    return 0;
}
