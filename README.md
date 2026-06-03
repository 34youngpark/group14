이항석
# Multivariable Newton-Raphson

고급프로그래밍 14조 · 담당: LI HANGXI (2026193019)

비선형 연립방정식(non-linear system)의 해를 구하기 위한 Multivariable Newton-Raphson Method를 C언어로 구현한 브랜치입니다.

## 파일 구성

| 파일               | 역할                        |
| ---------------- | ------------------------- |
| Newton-Raphson.c | Newton-Raphson 알고리즘 핵심 구현 |
| Newton-main.c    | 사용자 입력형 CLI 인터페이스         |
| Newton-Raphson.h | 함수 선언 헤더 파일               |
| matrix.h         | 팀 공용 자료구조 및 API 선언        |
| matrix.c         | 팀 공용 라이브러리 구현             |
| newton_test      | 컴파일된 실행 파일                |

## 알고리즘 개요

본 프로그램은 다음 비선형 연립방정식을 풉니다.

f₁(x,y) = x² + y² − 4

f₂(x,y) = x − y

해는 다음과 같습니다.

(x,y) = (√2, √2)

Newton-Raphson Method는 반복적으로 Jacobian 행렬을 계산하여 해를 근사합니다.

### 1단계 — 함수값 계산

현재 추정값 x에 대해 F(x)를 계산합니다.

### 2단계 — Jacobian 계산

수치 미분(Finite Difference)을 이용하여 Jacobian 행렬 J(x)를 계산합니다.

J(i,j) = dFi/dxj

### 3단계 — 선형 시스템 풀이

다음 식을 풉니다.

J(x)Δx = −F(x)

가우스 소거법(Gaussian Elimination with Partial Pivoting)을 사용하여 Δx를 계산합니다.

### 4단계 — 해 갱신

x ← x + Δx

를 수행하고 수렴할 때까지 반복합니다.

## 함수 인터페이스

```c
void multivariable_newton(Vector* x);
```

매개변수

| 매개변수 | 설명             |
| ---- | -------------- |
| x    | 초기 추정값을 저장한 벡터 |

반환값

없음 (입력 벡터 x가 최종 해로 갱신됨)

## 빌드 및 실행

```bash
gcc -o newton_test Newton-main.c Newton-Raphson.c matrix.c -lm
./newton_test
```

실행 예시

```text
Enter initial guess.

x = 1
y = 1

Iteration 0
...

Converged!

Final Solution:
[1.414214 1.414214]
```

## 한계

현재 구현은 2변수 비선형 연립방정식에 대해 작성되었습니다.

Jacobian은 수치 미분(Finite Difference) 방식으로 계산하므로 H 값 선택에 따라 정확도에 영향을 받을 수 있습니다.

최대 반복 횟수(MAX_ITER) 내에 수렴하지 못하면 "Did not converge." 메시지를 출력하고 종료합니다.
