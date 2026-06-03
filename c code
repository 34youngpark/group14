#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h" 

#define MAXITER 100 
#define TOLERANCE 1e-6 
#define H 1e-6 

/* ========================= API 확장용 static 헬퍼 함수 정의 ========================= */

// 벡터 뺄셈: r = a - b (새 벡터 동적 할당)
static Vector* vec_subtract(const Vector* a, const Vector* b) 
{
    Vector* r = vec_create(a->size);
    if (r == NULL) return NULL;
    for (int i = 0; i < a->size; i++) {
        r->data[i] = a->data[i] - b->data[i];
    }
    return r; 
} 

// 벡터 덧셈: r = a + b (새 벡터 동적 할당)
static Vector* vec_addition(const Vector* a, const Vector* b) 
{
    Vector* r = vec_create(a->size);
    if (r == NULL) return NULL;
    for (int i = 0; i < a->size; i++) {
        r->data[i] = a->data[i] + b->data[i];
    }
    return r; 
} 

// 벡터 스칼라 곱: r = c * v (새 벡터 동적 할당)
static Vector* vec_scalar_mul(const Vector* v, double c) 
{
    Vector* r = vec_create(v->size);
    if (r == NULL) return NULL;
    for (int i = 0; i < v->size; i++) {
        r->data[i] = c * v->data[i];
    }
    return r; 
}

/* ========================= 가우시안 소거법 (Linear Solver) ========================= */
// 매개변수로 넘어온 A와 b를 내부 알고리즘 진행 과정에서 직접 수정(In-place modification)하므로,
// 원본 보호를 원할 경우 반드시 호출부에서 mat_copy나 복사본을 만들어 전달해야 합니다.
static Vector* solve_linear_system(Matrix* A, Vector* b) 
{
    int n = A->rows;
    Vector* x = vec_create(n);
    if (x == NULL) return NULL;
    
    // 1. 전진 소거 (Forward Elimination) 및 부분 피보팅
    for (int k = 0; k < n - 1; k++)
    { 
        int pivot = k;
        for (int i = k + 1; i < n; i++) {
            if (fabs(A->data[i][k]) > fabs(A->data[pivot][k])) {
                pivot = i;
            }
        }
        
        if (pivot != k) {
            mat_swap_rows(A, k, pivot); // API 제공 함수 활용 (O(1))
            double temp = b->data[k];
            b->data[k] = b->data[pivot];
            b->data[pivot] = temp;
        }
        
        for (int i = k + 1; i < n; i++) {
            if (fabs(A->data[k][k]) < 1e-12) continue; // 제로 피벗 방지
            double factor = A->data[i][k] / A->data[k][k]; 
            for (int j = k; j < n; j++) {
                A->data[i][j] -= factor * A->data[k][j];
            }
            b->data[i] -= factor * b->data[k];
        }
    }
    
    // 2. 후진 대입 (Backward Substitution)
    for (int i = n - 1; i >= 0; i--) {
        if (fabs(A->data[i][i]) < 1e-12) {
            // 특이 행렬(Singular Matrix) 대응 안전장치
            x->data[i] = 0.0;
            continue;
        }
        double sum = b->data[i]; 
        for (int j = i + 1; j < n; j++) {
            sum -= A->data[i][j] * x->data[j];
        }
        x->data[i] = sum / A->data[i][i];
    }
    return x; 
} 

/* ========================= 대상 비선형 함수 F(x) ========================= */ 
static Vector* function_F(const Vector* x) 
{ 
    if (x == NULL || x->size != 2) return NULL; 
    Vector* f = vec_create(2); 
    if (f == NULL) return NULL; 
    
    double x0 = x->data[0]; 
    double x1 = x->data[1];
    
    // 수렴 테스트용 시스템식 설정
    f->data[0] = x0 * x0 + x1 * x1 - 4.0; 
    f->data[1] = x0 - x1; 
    return f; 
} 

/* ========================= 수치 자코비안 계산 (초기 B0용) ========================= */
static Matrix* compute_jacobian(const Vector* x) 
{ 
    int n = x->size; 
    Matrix* J = mat_create(n, n); 
    if (J == NULL) return NULL;
    
    Vector* fx = function_F(x); 
    if (fx == NULL) { 
        mat_free(J); 
        return NULL; 
    } 

    for (int j = 0; j < n; j++) 
    { 
        Vector* xh = vec_create(n);
        if (xh == NULL) {
            vec_free(fx); mat_free(J); return NULL;
        }
        for (int i = 0; i < n; i++) {
            xh->data[i] = x->data[i];
        }
        xh->data[j] += H; 
        
        Vector* fxh = function_F(xh);
        if (fxh == NULL) {
            vec_free(xh); vec_free(fx); mat_free(J); return NULL;
        }
        
        for (int i = 0; i < n; i++) {
            J->data[i][j] = (fxh->data[i] - fx->data[i]) / H;
        }
        vec_free(xh);
        vec_free(fxh);
    }
    vec_free(fx); 
    return J;
}

/* ========================= 브로이덴 메서드 본체 ========================= */ 
void broyden_method(Vector* x)
{ 
    int n = x->size; 
    Matrix* B = compute_jacobian(x);
    if (B == NULL) return;

    for (int iter = 0; iter < MAXITER; iter++) 
    {
        Vector* Fx = function_F(x);
        if (Fx == NULL) { mat_free(B); return; }

        printf("Iteration %d: ", iter);
        vec_print(x); // API 디버깅 출력 함수 활용
        
        // 1. 오차 판단 (L2 노름 기준)
        if (vec_norm(Fx) < TOLERANCE) 
        {
            printf("Converged!\n");
            vec_free(Fx);
            mat_free(B); 
            return;
        }

        // 2. 우변 선형시스템 풀이용 벡터 생성 (-F(x))
        Vector* minusF = vec_scalar_mul(Fx, -1.0);
        Matrix* Bcopy = mat_copy(B); // 원본 보존을 위한 깊은 복사
        
        // 중요 오류 수정: 가우시안 소거법 내부에서 minusF가 직접 변형되어 연산되므로, 
        // 껍데기만 생성하여 가짜 루프를 돌던 이전의 잘못된 `rhs` 변수를 완전히 들어냈습니다.
        Vector* s = solve_linear_system(Bcopy, minusF);
        
        // 3. 변수 업데이트 계산
        Vector* x_new = vec_addition(x, s);
        Vector* Fx_new = function_F(x_new);
        
        Vector* y = vec_subtract(Fx_new, Fx);
        Vector* Bs = mat_vec_multiply(B, s);
        Vector* y_minus_Bs = vec_subtract(y, Bs); 
        
        // 4. Broyden 공식에 의한 근사 자코비안 분모 확인 및 업데이트
        double denom = vec_dot_product(s, s); 
        if (fabs(denom) > 1e-12)
        {
            Matrix* update = vec_outer_product(y_minus_Bs, s);
            if (update != NULL) {
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        B->data[i][j] += update->data[i][j] / denom;
                    }
                }
                mat_free(update); 
            }
        }

        // 해 전이
        for (int i = 0; i < n; i++) {
            x->data[i] = x_new->data[i]; 
        }

        // 5. 메모리 관리 지침 엄수: 루프 내 동적 할당 객체 완전 소멸
        vec_free(Fx); 
        vec_free(minusF);
        vec_free(s);
        vec_free(x_new);
        vec_free(Fx_new); 
        vec_free(y); 
        vec_free(Bs); 
        vec_free(y_minus_Bs);
        mat_free(Bcopy);
    }
    
    mat_free(B);
    printf("Did not converge.\n"); 
}

/* ========================= 메인 엔트리 포인트 ========================= */
int main(void)
{
    // 1차원 벡터 할당 규칙 적용
    Vector* x = vec_create(2);
    if (x == NULL) {
        printf("Memory allocation failed.\n");
        return 1;
    }

    printf("초기 추정값 x0, x1을 공백으로 구분하여 입력해 주세요 (예: 1.0 1.0):\n");
    if (scanf("%lf %lf", &x->data[0], &x->data[1]) != 2) {
        printf("올바르지 않은 입력 양식입니다.\n");
        vec_free(x);
        return 1;
    }

    printf("\n[연산 프로세스 기동]\n");
    broyden_method(x);

    printf("\n[최종 추정 오차 범위 내 결과 수렴]\n");
    printf("최종 해 벡터 x: ");
    vec_print(x);

    // 프로그램 종료 직전 최종 릴리스 확인
    vec_free(x);
    return 0;
}
