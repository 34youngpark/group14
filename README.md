## QR Algorithm with Shifts (박세영, 2024199021)

### 알고리즘 개요

1. **Hessenberg 축소**: Householder 반사를 이용해 입력 행렬을 upper Hessenberg 형태로 변환
2. **QR iteration**: Wilkinson shift + Givens rotation으로 고유값 수렴 가속
3. **Deflation**: 하단 서브대각 원소가 `tol` 이하로 수렴하면 고유값 확정 후 활성 부분 행렬 축소
4. **반환**: 실수 고유값(eigenvalue)들을 `Vector`로 반환

> ** 실수 **대칭 행렬**에 대해서만 올바르게 작동합니다. 비대칭 행렬은 복소 고유값을 가질 수 있어 지원하지 않습니다.

### 사용법

```c
#include "matrix.h"

Vector* eigenvalues = qr_algorithm(A, 1000, 1e-10);
// ...
vec_free(eigenvalues);
```

### 함수 시그니처

```c
Vector* qr_algorithm(const Matrix *A, int max_iter, double tol);
```

| 매개변수 | 설명 |
|---|---|
| `A` | 정방 실수 행렬 (원본 훼손 없음 — 내부적으로 복사) |
| `max_iter` | 최대 반복 횟수 (권장: 1000 이상) |
| `tol` | 수렴 판정 허용 오차 (권장: 1e-10) |

반환: 고유값을 담은 `Vector` (크기 n). 오류 시 `NULL` 반환.  
메모리 해제는 호출자가 `vec_free()`로 수행.
