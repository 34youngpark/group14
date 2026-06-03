최지훈
# LU Decomposition with Partial Pivoting 구현 설명

본 코드는 부분 피벗팅(Partial Pivoting)을 적용한 LU 분해를 수행하는 프로그램입니다. 입력된 정사각 행렬 A를 순열 행렬(P), 하삼각행렬(L), 상삼각행렬(U)로 분해하여 다음 식을 만족하도록 합니다.

PA = LU

먼저 입력 행렬이 정사각 행렬인지 확인하고.

```c
if (!is_square_matrix(A)) {
    fprintf(stderr, "Error: LU decomposition requires square matrix.\n");
    return -1;
}
```

이후 L과 P는 단위행렬로 초기화하고, U는 입력 행렬 A를 복사하여 생성합니다.

```c
*L = mat_identity(n);
*U = mat_copy(A);
*P = mat_identity(n);
```

각 열마다 부분 피벗팅을 수행하여 절댓값이 가장 큰 원소를 피벗으로 선택합니다.

```c
int pivot_row = find_pivot_row(*U, k, k);
```

피벗 행이 현재 행과 다를 경우 U와 P의 행을 교환하며, 이미 계산된 L의 값도 함께 교환합니다.

```c
mat_swap_rows(*U, k, pivot_row);
mat_swap_rows(*P, k, pivot_row);
```

이후 가우스 소거법을 이용하여 피벗 아래 원소를 제거합니다.

double factor = (*U)->data[i][k] / (*U)->data[k][k];

(*L)->data[i][k] = factor;

계산된 factor는 L 행렬에 저장되고, 소거 결과는 U 행렬에 반영됩니다.

(*U)->data[i][j] -= factor * (*U)->data[k][j];

모든 반복이 끝나면 최종적으로 P, L, U 행렬이 생성되며, 이를 통해 입력 행렬 A를 다음과 같이 표현할 수 있습니다.

PA = LU

또한 피벗 값이 0에 가까운 경우 특이행렬(Singular Matrix)로 판단하여 오류를 발생시키도록 구현하였습니.
