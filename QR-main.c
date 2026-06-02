#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

// 엔진 파일(QR-algorithm.c)에 있는 함수를 불러옵니다.
extern Vector* qr_algorithm(const Matrix *A, int max_iter, double tol);

int main() {
    int n;

    printf("==================================================\n");
    printf("   사용자 입력형 QR Algorithm 고유값 계산기\n");
    printf("==================================================\n\n");

    // 1. 행렬의 크기 입력 받기
    printf("행렬의 크기(N x N)에서 N을 입력하세요 (예: 3): ");
    if (scanf("%d", &n) != 1 || n <= 0) {
        fprintf(stderr, "올바른 양의 정수를 입력해주세요.\n");
        return 1;
    }

    // 2. 행렬 생성
    Matrix *A = mat_create(n, n);
    if (A == NULL) {
        fprintf(stderr, "행렬 생성에 실패했습니다.\n");
        return 1;
    }

    // 3. 사용자로부터 행렬 성분 입력 받기
    printf("\n--- 행렬의 각 성분을 차례대로 입력하세요 ---\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("A[%d][%d] 성분 입력: ", i, j);
            // %lf는 double 자료형을 입력받을 때 사용합니다.
            if (scanf("%lf", &A->data[i][j]) != 1) {
                fprintf(stderr, "숫자 입력에 오류가 발생했습니다.\n");
                mat_free(A);
                return 1;
            }
        }
    }

    // 4. 입력된 행렬 확인차 출력
    printf("\n[확인] 입력된 행렬 A:\n");
    for (int i = 0; i < n; i++) {
        printf("    [ ");
        for (int j = 0; j < n; j++) {
            printf("%8.4f ", A->data[i][j]);
        }
        printf("]\n");
    }
    printf("\n");

    // 5. QR 알고리즘 실행
    printf("[진행] 고유값 계산 중 (max_iter: 1000, tol: 1e-9)...\n");
    Vector *eigenvalues = qr_algorithm(A, 1000, 1e-9);

    // 6. 결과 출력
    if (eigenvalues != NULL) {
        printf("\n 계산 완료! 고유값(Eigenvalues) 결과:\n");
        for (int i = 0; i < n; i++) {
            printf("    λ[%d] = %11.8f\n", i, eigenvalues->data[i]);
        }
        vec_free(eigenvalues);
    } else {
        printf("\n❌ 오류: 고유값 계산에 실패했습니다.\n");
    }

    // 메모리 정리
    mat_free(A);
    printf("==================================================\n");
    return 0;
}
