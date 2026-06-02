# Group 14 — Linear Algebra Solver: Plain-English Guide

> For: 안의빈 · Power Iteration & Inverse Power Iteration (`eigen.c`)

---

## 1. What is this project, in one sentence?

This project finds **eigenvalues and eigenvectors** of a matrix using two iterative algorithms, written in C with a shared matrix/vector library.

---

## 2. What is an eigenvalue / eigenvector? (Plain English)

Think of a matrix as a **transformation** — it rotates and stretches vectors in space.

Most vectors change direction when you multiply them by a matrix. But a few special vectors only **scale** — they don't rotate. Those are eigenvectors.

```
A · v = λ · v

A  = the matrix
v  = the eigenvector  (special direction — doesn't rotate)
λ  = the eigenvalue   (how much it stretches/shrinks along that direction)
```

**Example:** Matrix `A = [[4,1],[2,3]]` has:
- Eigenvalue λ₁ = **5**, eigenvector **[1, 1] / √2** → stretched 5×
- Eigenvalue λ₂ = **2**, eigenvector **[1,-2] / √5** → stretched 2×

Eigenvalues tell you the "natural frequencies" of the matrix — useful in physics, machine learning, graph theory, etc.

---

## 3. Power Iteration (`power_iteration` in eigen.c)

**Goal:** Find the eigenvalue with the **largest absolute value** (the "dominant" eigenvalue).

### The Idea

Any vector can be written as a combination of eigenvectors:

```
q₀ = c₁v₁ + c₂v₂ + ... + cₙvₙ
```

When you multiply by A repeatedly, each component grows by its eigenvalue:

```
Aᵏq₀ = c₁λ₁ᵏv₁ + c₂λ₂ᵏv₂ + ...
```

After many steps, the term with the **biggest λ dominates** — the others become negligible by comparison. The vector converges to the dominant eigenvector.

### Step-by-Step Algorithm

```
1. Start:   q = [1, 0, ..., 0]   ← initial guess (any nonzero vector)

2. Multiply:  z = A · q           ← stretch q by A

3. Estimate: λ = qᵀ · z           ← Rayleigh Quotient: how much did q grow in its own direction?

4. Normalize: q = z / ‖z‖         ← keep q as a unit vector (prevent overflow)

5. Check:  if |λ_new − λ_old| < tol → DONE
           else → go to step 2
```

### Code Location: `eigen.c:93`

```c
Vector* power_iteration(const Matrix* A, int max_iter, double tol,
                        double* eigenvalue, int* iters_done);
```

| Parameter    | What it is                                      |
|-------------|--------------------------------------------------|
| `A`          | The input matrix                                 |
| `max_iter`   | Max iterations before giving up (default 10000) |
| `tol`        | Convergence threshold (default 1e-10)           |
| `eigenvalue` | **Output**: the dominant eigenvalue              |
| `iters_done` | **Output**: how many iterations it took          |
| return value | The eigenvector (caller must `vec_free` it)      |

---

## 4. Inverse Power Iteration (`inverse_power_iteration` in eigen.c)

**Goal:** Find the eigenvalue **closest to a given number σ** (sigma).

If σ = 0, this finds the eigenvalue with the **smallest absolute value**.

### The Key Trick: Shift + Invert

If the eigenvalues of A are λ₁, λ₂, ..., then the eigenvalues of `(A − σI)⁻¹` are:

```
1/(λ₁−σ),  1/(λ₂−σ),  ...
```

The eigenvalue of A closest to σ produces the **largest** value of `1/(λ−σ)`. So if we do Power Iteration on `(A−σI)⁻¹`, it converges to the closest eigenvalue. 

### Why Not Compute the Inverse Directly?

Computing a matrix inverse is expensive and numerically unstable. Instead, at each step we just **solve a linear system**:

```
(A − σI) · z = q     ←  same as  z = (A−σI)⁻¹ · q
```

The `solve_linear_system` helper in `eigen.c` does this using **Gaussian elimination with partial pivoting**.

### Step-by-Step Algorithm

```
1. Build: B = A − σ·I              ← shift the matrix

2. Start: q = [1, 0, ..., 0]

3. Solve: B · z = q                ← Gaussian elimination (eigen.c:15)

4. Rayleigh Quotient: μ = qᵀ · z   ← μ ≈ 1/(λ − σ)

5. Recover: λ = σ + 1/μ            ← the actual eigenvalue of A

6. Normalize: q = z / ‖z‖

7. Check: if |λ_new − λ_old| < tol → DONE
          else → go to step 3
```

### Code Location: `eigen.c:169`

```c
Vector* inverse_power_iteration(const Matrix* A, double sigma,
                                int max_iter, double tol,
                                double* eigenvalue, int* iters_done);
```

---

## 5. What is the Rayleigh Quotient?

It's the formula for estimating eigenvalue from a current vector:

```
λ ≈ (qᵀAq) / (qᵀq)
```

Since q is always normalized (‖q‖ = 1), the denominator is 1, so:

```
λ = qᵀ · z    where z = A·q
```

This is a scalar — it measures how much the vector grew "in its own direction" (the dot product captures the projection).

---

## 6. The Gaussian Elimination Helper (`solve_linear_system`)

This is inside `eigen.c:15` and is only used by Inverse Power Iteration.

It solves `A·x = b` using:

1. **Build augmented matrix** `[A | b]` — attach b as an extra column
2. **Forward elimination** — zero out entries below the diagonal, with **partial pivoting** (swap rows to put the largest value on the diagonal first, for numerical stability)
3. **Back substitution** — solve from bottom row upward

```
[4  1 | 5]        [4  1 | 5]         [4  1 | 5]
[2  3 | 7]  →  0  row ops  →   [0  2.5| 4.5]
                                 ↓ back sub → x = [0.8, 1.8]
```

---

## 7. File Structure

```
group14-main/
├── matrix.h        ← shared API: Matrix + Vector types, function declarations
├── matrix.c        ← shared library: all basic matrix/vector operations
│
├── eigen.c         ← CORE: Power Iteration + Inverse Power Iteration
├── eigen-main.c    ← CLI program: user types in a matrix, picks algorithm, sees result
│
├── testmain.c      ← simple test of the matrix library (A*v multiplication)
│
├── visualizer.html ← interactive visual explanation (open in browser)
├── GUIDE.md        ← this file
└── README-eigen.md ← original Korean README for this component
```

---

## 8. Key Data Structures

### `Matrix` (matrix.h:11)

```c
typedef struct {
    int rows, cols;
    double *block;    // flat array of all values (rows*cols doubles)
    double **data;    // pointers into block → lets you write A->data[i][j]
} Matrix;
```

The `block` is one contiguous chunk of memory. The `data` pointers point into it so you can use `A->data[row][col]` syntax conveniently. This is a common C trick for 2D arrays.

### `Vector` (matrix.h:18)

```c
typedef struct {
    int size;
    double *data;    // flat array: v->data[i]
} Vector;
```

### Memory Rule

**Every function that returns a `Matrix*` or `Vector*` allocates new memory.** The caller is responsible for freeing it:

```c
Vector* result = power_iteration(A, 0, 0.0, &lambda, &iters);
// ... use result ...
vec_free(result);   // ← MUST do this, or memory leak
mat_free(A);
```

---

## 9. How to Build & Run

```bash
# Build the CLI tool
gcc -o eigen_test eigen-main.c eigen.c matrix.c -lm

# Run it
./eigen_test
```

**Example session:**
```
행렬의 크기: 2
A[0][0] = 4   A[0][1] = 1
A[1][0] = 2   A[1][1] = 3

알고리즘 선택: 1   (Power Iteration)

→ λ = 5.00000000  (27 iterations)
→ v = [0.7071, 0.7071]   (= [1,1]/√2)
```

**Test the matrix library:**
```bash
gcc -o testmain testmain.c matrix.c -lm
./testmain
```

**Open the visual:**
```
open visualizer.html    # macOS
# or just double-click the file in Finder
```

---

## 10. Limitations

| Situation | What happens |
|---|---|
| Complex eigenvalues (non-symmetric matrix) | Power Iteration won't converge |
| Two eigenvalues with the same magnitude | Power Iteration oscillates, no convergence |
| σ exactly equals an eigenvalue | Inverse Power: `(A−σI)` becomes singular → crash |
| Very large/small values in A | May need more iterations |

---

## 11. Quick Reference: What calls what

```
eigen-main.c   (user input + output)
    ↓  calls
eigen.c        power_iteration()
               inverse_power_iteration()
                   ↓  calls internally
               solve_linear_system()   (static — not exposed)
    ↓  both use
matrix.c       mat_vec_multiply()
               vec_dot_product()
               vec_norm()
               mat_create() / mat_free() / vec_create() / vec_free()
               mat_identity() / mat_copy()
```
