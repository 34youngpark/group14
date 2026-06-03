#include <stdio.h>
#include "matrix.h"
#include "Newton-Raphson.h"

int main()
{
    Vector* x = vec_create(2);

    if (x == NULL)
    {
        printf("Memory allocation failed.\n");
        return 1;
    }

    printf("========================================\n");
    printf(" Multivariable Newton-Raphson Solver\n");
    printf("========================================\n");

    printf("This program solves:\n");
    printf("f1(x,y) = x^2 + y^2 - 4\n");
    printf("f2(x,y) = x - y\n\n");

    printf("Enter initial guess.\n");

    printf("x = ");
    scanf("%lf", &x->data[0]);

    printf("y = ");
    scanf("%lf", &x->data[1]);

    printf("\n");

    multivariable_newton(x);

    printf("\nFinal Solution:\n");
    vec_print(x);

    vec_free(x);

    return 0;
}
