#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"

// qr-algorithm.c 파일의 함수를 가져옵니다.
extern Vector* qr_algorithm(const Matrix *A, int max_iter, double tol);

int main() {
    int n = 3;
    printf("==================================================\n");
    printf("   QR Algorithm 고유값(Eigenvalues) 계산 테스트\n");
    printf("==================================================\n\n");

    // 1. 3x3 테스트용 행렬 생성
    Matrix *A = mat_create(n, n);
    if (A == NULL) {
        fprintf(stderr, "행렬 생성에 실패했습니다.\n");
        return 1;
    }

    /* * 테스트용 대칭 행렬 설정 
     * A = [  2.0, -1.0,  0.0 ]
     * [ -1.0,  2.0, -1.0 ]
     * [  0.0, -1.0,  2.0 ]
     * * 이론적 실제 고유값:
     * λ1 = 2 + sqrt(2) ≒ 3.41421356
     * λ2 = 2.00000000
     * λ3 = 2 - sqrt(2) ≒ 0.58578643
     */
    A->data[0][0] = 2.0;  A->data[0][1] = -1.0; A->data[0][2] = 0.0;
    A->data[1][0] = -1.0; A->data[1][1] = 2.0;  A->data[1][2] = -1.0;
    A->data[2][0] = 0.0;  A->data[2][1] = -1.0; A->data[2][2] = 2.0;

    // 2. 입력 행렬 출력
    printf("[1] 입력 행렬 A:\n");
    for (int i = 0; i < n; i++) {
        printf("    [ ");
        for (int j = 0; j < n; j++) {
            printf("%8.4f ", A->data[i][j]);
        }
        printf("]\n");
    }
    printf("\n");

    // 3. QR 알고리즘 엔진 작동 (최대 1000번 반복, 허용 오차 1e-9)
    printf("[2] QR 알고리즘 연산 시작 (max_iter: 1000, tol: 1e-9)...\n");
    Vector *eigenvalues = qr_algorithm(A, 1000, 1e-9);

    // 4. 결과 출력 및 검증
    if (eigenvalues != NULL) {
        printf("\n[3] 연산 결과 (계산된 고유값):\n");
        for (int i = 0; i < n; i++) {
            printf("    λ[%d] = %11.8f\n", i, eigenvalues->data[i]);
        }
        
        printf("\n[4] 이론적 실제값 (참고용 정답지):\n");
        printf("    λ_1 ≒  3.41421356\n");
        printf("    λ_2  =  2.00000000\n");
        printf("    λ_3 ≒  0.58578643\n");
        printf("\n💡 계산된 값과 실제값이 거의 일치하면 알고리즘이 성공적으로 구현된 것입니다.\n");
        
        // 메모리 해제
        vec_free(eigenvalues);
    } else {
        printf("\n❌ 오류: 고유값 계산에 실패했거나 NULL이 반환되었습니다.\n");
    }

    mat_free(A);
    printf("==================================================\n");
    return 0;
}
