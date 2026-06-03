/*
    Multivariate Newton-Raphson Solver
    C Linear Algebra Project
*/

#include <stdio.h>
#include <math.h>
#include "matrix.h"
#include "Newton-Raphson.h"

#define MAX_ITER 100
#define TOLERANCE 1e-6
#define H 1e-6

/*
    Nonlinear Function F(x)

    f1(x,y) = x^2 + y^2 - 4
    f2(x,y) = x - y

    Solution:
        (sqrt(2), sqrt(2))
*/

static Vector* function_F(const Vector* x)
{
    if (x == NULL || x->size != 2)
    {
        return NULL;
    }

    Vector* f = vec_create(2);

    if (f == NULL)
    {
        return NULL;
    }

    double x0 = x->data[0];
    double x1 = x->data[1];

    f->data[0] = x0 * x0 + x1 * x1 - 4.0;
    f->data[1] = x0 - x1;

    return f;
}

/*
    Numerical Jacobian Matrix

    J(i,j) = dFi/dxj
*/

static Matrix* compute_jacobian(const Vector* x)
{
    int n = x->size;

    Matrix* J = mat_create(n, n);

    if (J == NULL)
    {
        return NULL;
    }

    Vector* fx = function_F(x);

    if (fx == NULL)
    {
        mat_free(J);
        return NULL;
    }

    for (int j = 0; j < n; j++)
    {
        Vector* xh = vec_create(n);

        if (xh == NULL)
        {
            vec_free(fx);
            mat_free(J);
            return NULL;
        }

        /* copy x */
        for (int i = 0; i < n; i++)
        {
            xh->data[i] = x->data[i];
        }

        xh->data[j] += H;

        Vector* fxh = function_F(xh);

        if (fxh == NULL)
        {
            vec_free(xh);
            vec_free(fx);
            mat_free(J);
            return NULL;
        }

        for (int i = 0; i < n; i++)
        {
            J->data[i][j] =
                (fxh->data[i] - fx->data[i]) / H;
        }

        vec_free(xh);
        vec_free(fxh);
    }

    vec_free(fx);

    return J;
}

/*
    Gaussian Elimination with Partial Pivoting

    Solve:
        A*x = b
*/

static Vector* solve_linear_system(Matrix* A, Vector* b)
{
    int n = A->rows;

    Vector* x = vec_create(n);

    if (x == NULL)
    {
        return NULL;
    }

    /* Forward Elimination */
    for (int k = 0; k < n - 1; k++)
    {
        int max_row = k;

        /* Partial Pivoting */
        for (int i = k + 1; i < n; i++)
        {
            if (fabs(A->data[i][k]) >
                fabs(A->data[max_row][k]))
            {
                max_row = i;
            }
        }

        /* Row Swap */
        if (max_row != k)
        {
            mat_swap_rows(A, k, max_row);

            double temp = b->data[k];
            b->data[k] = b->data[max_row];
            b->data[max_row] = temp;
        }

        /* Singular Matrix Check */
        if (fabs(A->data[k][k]) < 1e-12)
        {
            printf("Singular matrix detected.\n");

            vec_free(x);
            return NULL;
        }

        /* Elimination */
        for (int i = k + 1; i < n; i++)
        {
            double factor =
                A->data[i][k] / A->data[k][k];

            for (int j = k; j < n; j++)
            {
                A->data[i][j] -=
                    factor * A->data[k][j];
            }

            b->data[i] -= factor * b->data[k];
        }
    }

    /* Back Substitution */
    for (int i = n - 1; i >= 0; i--)
    {
        if (fabs(A->data[i][i]) < 1e-12)
        {
            printf("Singular matrix detected.\n");

            vec_free(x);
            return NULL;
        }

        double sum = b->data[i];

        for (int j = i + 1; j < n; j++)
        {
            sum -= A->data[i][j] * x->data[j];
        }

        x->data[i] =
            sum / A->data[i][i];
    }

    return x;
}

/*
    Multivariable Newton-Raphson Method

    J(x) * dx = -F(x)

    x_new = x + dx
*/

void multivariable_newton(Vector* x)
{
    if (x == NULL)
    {
        return;
    }

    int n = x->size;

    for (int iter = 0; iter < MAX_ITER; iter++)
    {
        Vector* fx = function_F(x);

        if (fx == NULL)
        {
            printf("Function evaluation failed.\n");
            return;
        }

        printf("Iteration %d\n", iter);

        vec_print(x);

        printf("Residual Norm = %f\n",
               vec_norm(fx));

        /* Convergence Check */
        if (vec_norm(fx) < TOLERANCE)
        {
            printf("Converged!\n");

            vec_free(fx);
            return;
        }

        /* Jacobian */
        Matrix* J = compute_jacobian(x);

        if (J == NULL)
        {
            printf("Jacobian computation failed.\n");

            vec_free(fx);
            return;
        }

        /* -F(x) */
        Vector* minus_fx = vec_create(n);

        if (minus_fx == NULL)
        {
            vec_free(fx);
            mat_free(J);
            return;
        }

        for (int i = 0; i < n; i++)
        {
            minus_fx->data[i] = -fx->data[i];
        }

        /*
            Copy matrix and vector
            because Gaussian elimination
            modifies them
        */

        Matrix* J_copy = mat_copy(J);

        Vector* b_copy = vec_create(n);

        if (J_copy == NULL || b_copy == NULL)
        {
            vec_free(fx);
            vec_free(minus_fx);

            mat_free(J);

            return;
        }

        for (int i = 0; i < n; i++)
        {
            b_copy->data[i] =
                minus_fx->data[i];
        }

        /* Solve J*dx = -F */
        Vector* dx =
            solve_linear_system(J_copy, b_copy);

        if (dx == NULL)
        {
            printf("Linear solve failed.\n");

            vec_free(fx);
            vec_free(minus_fx);
            vec_free(b_copy);

            mat_free(J);
            mat_free(J_copy);

            return;
        }

        /* x = x + dx */
        for (int i = 0; i < n; i++)
        {
            x->data[i] += dx->data[i];
        }

        /* Additional convergence check */
        if (vec_norm(dx) < TOLERANCE)
        {
            printf("Step size converged.\n");

            vec_free(fx);
            vec_free(minus_fx);
            vec_free(b_copy);
            vec_free(dx);

            mat_free(J);
            mat_free(J_copy);

            return;
        }

        /* Memory Free */
        vec_free(fx);
        vec_free(minus_fx);
        vec_free(b_copy);
        vec_free(dx);

        mat_free(J);
        mat_free(J_copy);

        printf("\n");
    }

    printf("Did not converge.\n");
}
