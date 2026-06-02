# Power Iteration & Inverse Power Iteration

고급프로그래밍 14조 · 담당: 안의빈 (2024199084)

실수 정방 행렬의 **고유값(eigenvalue)** 과 **고유벡터(eigenvector)** 를 구하는
거듭제곱 반복법(Power Iteration) 및 역반복법(Inverse Power Iteration)을 C언어로 구현한 브랜치입니다.

---

## 파일 구성

| 파일 | 역할 |
|---|---|
| `eigen.c` | 알고리즘 핵심 구현 (`power_iteration`, `inverse_power_iteration` 및 내부 헬퍼) |
| `eigen-main.c` | 사용자 입력형 CLI 인터페이스. 행렬과 알고리즘을 직접 선택하여 결과를 출력 |
| `matrix.h` | 팀 공용 자료구조 및 API 선언 |
| `matrix.c` | 팀 공용 라이브러리 구현 |

---

## 알고리즘 개요

### 1. Power Iteration (거듭제곱 반복법)

행렬 A의 **지배 고유값** (절댓값이 가장 큰 고유값)을 구합니다.

1. 초기 벡터 `q = [1, 0, ..., 0]`
2. `z = A * q` 계산
3. Rayleigh Quotient `λ = qᵀz` 로 고유값 추정
4. `q = z / ‖z‖` 로 정규화
5. `|λ_new − λ_old| < tol` 이면 수렴, 아니면 2로 돌아감

### 2. Inverse Power Iteration (역반복법)

행렬 A의 **shift σ에 가장 가까운 고유값**을 구합니다.  
`σ = 0.0` 으로 설정하면 절댓값이 가장 작은 고유값을 구합니다.

1. `B = A − σI` 구성
2. 초기 벡터 `q = [1, 0, ..., 0]`
3. `(A − σI) * z = q` 를 가우스 소거법(부분 피벗)으로 풀기
4. Rayleigh Quotient `μ = qᵀz`
5. 실제 고유값 복원: `λ = σ + 1/μ`
6. `|λ_new − λ_old| < tol` 이면 수렴, 아니면 3으로 돌아감

---

## 함수 인터페이스

```c
#include "matrix.h"

Vector* power_iteration(const Matrix* A, int max_iter, double tol,
                        double* eigenvalue, int* iters_done);

Vector* inverse_power_iteration(const Matrix* A, double sigma,
                                int max_iter, double tol,
                                double* eigenvalue, int* iters_done);
```

| 매개변수 | 설명 |
|---|---|
| `A` | 입력 정방 행렬 (원본 훼손 없음) |
| `sigma` | shift 값 (Inverse Power Iteration 전용) |
| `max_iter` | 최대 반복 횟수 (0이면 기본값 10000 사용) |
| `tol` | 수렴 판정 허용 오차 (0.0이면 기본값 1e-10 사용) |
| `eigenvalue` | 출력용 고유값 포인터 (필수) |
| `iters_done` | 출력용 반복 횟수 포인터 (NULL 가능) |

반환값: 고유벡터를 담은 `Vector`. 오류 시 `NULL` 반환.  
반환된 `Vector`의 메모리 해제는 호출자가 `vec_free()`로 수행합니다.

---

## 빌드 및 실행

```bash
gcc -o eigen_test eigen-main.c eigen.c matrix.c -lm
./eigen_test
```

실행 예시:

```
행렬의 크기(N x N)에서 N을 입력하세요 (예: 3): 2
A[0][0] 성분 입력: 4
A[0][1] 성분 입력: 1
A[1][0] 성분 입력: 2
A[1][1] 성분 입력: 3

실행할 알고리즘을 선택하세요:
  1. Power Iteration         (지배 고유값)
  2. Inverse Power Iteration (sigma에 가장 가까운 고유값)
선택 (1 or 2): 1

  λ =  5.00000000  (반복 횟수: 27)
  v[0] =  0.70710678
  v[1] =  0.70710678
```

---

## 한계

- **복소 고유값**을 갖는 행렬(비대칭)은 수렴이 보장되지 않습니다.
- 지배 고유값이 **2개 이상 동일한 절댓값**을 가지면 Power Iteration이 수렴하지 않을 수 있습니다.
- Inverse Power Iteration에서 `σ`가 실제 고유값과 정확히 일치하면 `(A − σI)`가 특이 행렬이 되어 실패합니다.
