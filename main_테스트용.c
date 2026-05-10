#include <stdio.h>
#include "matrix.h"

int main() {
    printf("=== 라이브러리 정상 작동 테스트 ===\n\n");

    // 1. 2x2 행렬 생성 및 값 넣기
    Matrix* A = mat_create(2, 2);
    if (A != NULL) {
        A->data[0][0] = 1.0; A->data[0][1] = 2.0;
        A->data[1][0] = 3.0; A->data[1][1] = 4.0;
    }

    // 2. 크기가 2인 벡터 생성 및 값 넣기
    Vector* v = vec_create(2);
    if (v != NULL) {
        v->data[0] = 1.0;
        v->data[1] = 2.0;
    }

    // 3. 입력값 출력 확인
    printf("[Matrix A]\n");
    mat_print(A);

    printf("[Vector v]\n");
    vec_print(v);

    // 4. 행렬-벡터 곱셈 (A * v) 연산 테스트
    Vector* result = mat_vec_multiply(A, v);
    printf("[Result: A * v]\n");
    vec_print(result);

    // 5. 메모리 해제 (매우 중요!)
    mat_free(A);
    vec_free(v);
    vec_free(result);

    printf("테스트가 성공적으로 종료되었습니다!\n");
    return 0;
}
