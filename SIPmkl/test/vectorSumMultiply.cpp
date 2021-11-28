#include <mkl.h>
#include <iostream>
using namespace std;

int scond1(
    float h[], int inch,
    float x[], int incx,
    float y[], int incy,
    int nh, int nx, int iy0, int ny)
{
    int status;
    VSLConvTaskPtr task;
    vslsConvNewTask1D(&task, VSL_CONV_MODE_DIRECT, nh, nx, ny);
    vslConvSetStart(task, &iy0);
    status = vslsConvExec1D(task, h, inch, x, incx, y, incy);
    vslConvDeleteTask(&task);
    return status;
}

int main(int argc, char const *argv[])
{
    float *a, *b, *c, *d;
    int n = 3;
    a = (float *)mkl_malloc(n * 1 * sizeof(float), 64);
    b = (float *)mkl_malloc(n * 1 * sizeof(float), 64);
    c = (float *)mkl_malloc(n * 1 * sizeof(float), 64);
    d = (float *)mkl_malloc((n + n - 1) * sizeof(float), 64);
    printf("a vector:\n");
    for (int i = 0; i < n; i++)
    {
        a[i] = i;
        printf("%2.0f", a[i]);
    }
    printf("\n");
    printf("b vector:\n");
    for (int i = 0; i < n; i++)
    {
        b[i] = i;
        printf("%2.0f", b[i]);
    }
    printf("\n");
    vsMul(n, a, b, c);
    printf("vector a*b:\n");
    for (int i = 0; i < n; i++)
    {
        printf("%3.0f", c[i]);
    }
    printf("\n");
    vsAdd(n, a, b, c);
    printf("vector a+b:\n");
    for (int i = 0; i < n; i++)
    {
        printf("%3.0f", c[i]);
    }
    printf("\n");
    printf("convoltion a b:\n");
    scond1(a, 1, b, 1, d, 1, n, n, 0, (n + n - 1));
    for (int i = 0; i < n+n-1; i++)
    {
        printf("%5.0f", d[i]);
    }
    printf("\n");
    mkl_free(a);
    mkl_free(b);
    mkl_free(c);
    mkl_free(d);
    return 0;
}
