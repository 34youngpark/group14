**고급프로그래밍 14조
C언어 기반 선형대수 솔버 프로젝트**

개요
  본 프로젝트는 C언어 기반의 선형대수 솔버 구현 프로젝트입니다. 공용으로 사용하는 matrix.c 파일을 기반으로 하여 조원 각각이 하나의 수치해석 솔버를 구현합니다.

규칙
  각자의 코드(linear.c, nonlinear.c 등) 간에 어떠한 함수 호출이나 의존성도 없어야 합니다. 오직 본인에게 할당된 알고리즘 로직에만 집중합니다.
  **공용 파일 수정 금지: matrix.c 및 matrix.h 파일은 절대 건드리지 않습니다.**



  

**팀원 및 담당 주제**

최지훈: LU Decomposition with Partial Pivoting 

전혜빈: Conjugate Gradient Method

이항석: Multiderivative(Multivariate) Newton-Raphson

김민교: Broyden's Method 

안의빈: Power Iteration, Inverse Power Iteration

박세영: QR Algorithm with Shifts 





**matrix.c API 사용설명서**

자신의 .c 파일 최상단에 #include "matrix.h"를 선언한 뒤, 아래의 함수들을 사용할 수 있습니다.

**추가로 필요한 함수는 다른 함수와 충돌하지 않도록 matrix.h와 겹치지 않는 이름으로, static 함수로 만들어 사용합니다.**

  1. 핵심 자료구조
    본 라이브러리는 0-based 인덱스를 사용하며, 빠른 연산을 위해 내부적으로 연속된 메모리 공간을 사용합니다. 값에 접근할 때는 직관적인 2차원/1차원 배열 방식을 사용하세요.

    행렬 (Matrix): m->data[row][col] 형태로 접근
    
    벡터 (Vector): v->data[index] 형태로 접근

  2. 메모리 관리 주의사항
    mat_create, vec_create를 포함해 연산의 결과값으로 반환된 모든 새로운 행렬과 벡터는 동적 할당된 메모리입니다. 함수 종료 전에 반드시 mat_free(), vec_free()를 호출하여 메모리를 해제해야 합니다.

  3. 사용 가능 함수 목록
     
     (1)생성 및 파괴
        ```Matrix* mat_create(int rows, int cols); // 0으로 초기화된 행렬 생성
     
        Vector* vec_create(int size);           // 0으로 초기화된 벡터 생성
     
        void mat_free(Matrix* m);               // 행렬 메모리 해제
     
        void vec_free(Vector* v);               // 벡터 메모리 해제
        ```
     
     (2) 행렬 기본 연산 (Matrix Operations)
     
        차원(Dimension)이 맞지 않으면 에러 메시지를 출력하고 NULL을 반환합니다.
       ``` Matrix* mat_add(const Matrix* a, const Matrix* b);      // A + B
     
        Matrix* mat_sub(const Matrix* a, const Matrix* b);      // A - B
     
        Matrix* mat_multiply(const Matrix* a, const Matrix* b); // A * B
     
        Matrix* mat_transpose(const Matrix* m);                 // 전치 행렬 (A^T)
     ```
       
     (3)벡터 및 혼합 연산 (Vector & Mixed Operations)
        ```Vector* mat_vec_multiply(const Matrix* m, const Vector* v);    // 행렬-벡터 곱셈 (Ax)
     
        double vec_dot_product(const Vector* v1, const Vector* v2);    // 내적 (Dot Product)
     
        Matrix* vec_outer_product(const Vector* v1, const Vector* v2); // 외적 (Outer Product - Broyden용)
     
        double vec_norm(const Vector* v);                              // 벡터의 L2 노름 (크기, 수렴 판정용)
     ```
        
     (4) 유틸리티 함수
        ```Matrix* mat_identity(int size);                         // 단위 행렬(I) 생성
     
        Matrix* mat_copy(const Matrix* m);                      // 원본 훼손 방지용 깊은 복사(Deep Copy)
     
        void mat_swap_rows(Matrix* m, int row1, int row2);      // 두 행 교환 (LU 부분 피벗용, O(1) 속도)
     
        void mat_print(const Matrix* m);                        // 디버깅용 행렬 콘솔 출력
     
        void vec_print(const Vector* v);                        // 디버깅용 벡터 콘솔 출력
     ```


