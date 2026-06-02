# QR Algorithm with Shifts

고급프로그래밍 14조 · 담당: 박세영 (2024199021)

실수 대칭 정방 행렬의 **고유값(eigenvalue)** 을 구하는 QR 알고리즘을 C언어로 구현한 브랜치입니다.

---

## 파일 구성

| 파일 | 역할 |
|---|---|
| `QR-algorithm.c` | QR 알고리즘 핵심 구현 (`qr_algorithm` 함수 및 내부 유틸리티) |
| `QR-main.c` | 사용자 입력형 CLI 인터페이스. 행렬을 직접 입력받아 고유값을 출력 |
| `matrix.h` | 팀 공용 자료구조 및 API 선언 |
| `matrix.c` | 팀 공용 라이브러리 구현 |
| `qr_test` | 컴파일된 실행 파일 |

---

## 알고리즘 개요

총 3단계로 동작합니다.

**1단계 — Hessenberg 축소**  
Householder 반사를 이용해 입력 행렬을 upper Hessenberg 형태로 변환합니다. 실수 대칭 행렬의 경우 결과는 tridiagonal 행렬이 됩니다. 이 전처리를 통해 이후 QR iteration의 연산량을 크게 줄입니다.

**2단계 — QR iteration (Wilkinson shift + Givens rotation)**  
Explicit Q 행렬 없이 Givens rotation만으로 `H ← Q^T H Q` 변환을 수행합니다(implicit Q). 수렴을 가속하기 위해 매 반복마다 활성 부분 행렬의 오른쪽 하단 2×2 블록 고유값 중 `H[m][m]`에 가장 가까운 값을 Wilkinson shift로 선택합니다.

**3단계 — Deflation**  
서브대각 원소 `|H[m][m-1]|`이 `tol × (|H[m-1][m-1]| + |H[m][m]|)` 이하로 수렴하면 해당 고유값을 확정하고 활성 부분 행렬 크기를 줄입니다. 모든 고유값이 확정되면 종료합니다.

---

## 함수 인터페이스

```c
#include "matrix.h"

Vector* qr_algorithm(const Matrix *A, int max_iter, double tol);
```

| 매개변수 | 설명 |
|---|---|
| `A` | 입력 정방 행렬 (내부에서 `mat_copy`로 복사하므로 원본 훼손 없음) |
| `max_iter` | 최대 반복 횟수 (권장: `1000` 이상) |
| `tol` | 수렴 판정 허용 오차 (권장: `1e-10`) |

반환값: 고유값을 담은 `Vector` (크기 n). 오류 시 `NULL` 반환.  
반환된 `Vector`의 메모리 해제는 호출자가 `vec_free()`로 수행합니다.

---

## 빌드 및 실행

```bash
gcc -o qr_test QR-main.c QR-algorithm.c matrix.c -lm
./qr_test
```

실행하면 다음과 같이 동작합니다.

```
행렬의 크기(N x N)에서 N을 입력하세요 (예: 3): 3
A[0][0] 성분 입력: 4
...
λ[0] =  5.00000000
λ[1] =  3.00000000
λ[2] = -1.00000000
```

---

## 한계

- **실수 대칭 행렬**에 대해서만 올바르게 동작합니다. 비대칭 행렬은 복소 고유값을 가질 수 있으며, 이 경우 수렴이 보장되지 않습니다.
- `max_iter` 내에 수렴하지 못하면 경고 메시지를 출력하고 현재 대각 원소를 근사값으로 반환합니다.
