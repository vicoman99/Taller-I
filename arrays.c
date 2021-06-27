/** @file arrays.c
 *  @brief Implementation of all the arrays and matrixes
 * 
 *  @author Victor Coman
 *  @date May 2021
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include "arrays.h"

// vectores

//calculo de la norma como longitud de un vector v de tamaño n
double norma(double *v, int n)
{
    double sum;
    int i;

    sum = 0.0;
    for (i = 0; i < n; i++)
        sum += v[i] * v[i];

    return (sqrt(sum));
}

//calculo del producto escalar de dos vectores
double dot(double *v1, int n1, double *v2, int n2)
{
    double sum;
    int i;

    if (n1 != n2)
    {
        printf("dot: n1 != n2\n");
        exit(EXIT_FAILURE);
    }

    sum = 0.0;
    for (i = 0; i < n1; i++)
        sum += v1[i] * v2[i];

    return (sum);
}

//creacion de un vector con el tamaño n
double *vector(int n)
{
    double *v;

    v = (double *)calloc(n, sizeof(double));
    if (v == NULL)
    {
        printf("vector: memory not allocated\n");
        exit(EXIT_FAILURE);
    }

    return v;
}

//calculo del producto vectorial
double *crossProd(double *vector_a, int n1, double *vector_b, int n2)
{
    if (n1 != n2 || n1 != 3)
    {
        printf("vector: ERROR, both must be of size 3 \n");
        exit(EXIT_FAILURE);
    }
    double *res;
    res = vector(3);
    res[0] = vector_a[1] * vector_b[2] - vector_a[2] * vector_b[1];
    res[1] = vector_a[2] * vector_b[0] - vector_a[0] * vector_b[2];
    res[2] = vector_a[0] * vector_b[1] - vector_a[1] * vector_b[0];

    return res;
}

//liberacion de un vector
void freeVector(double *v, int n)
{
    free(v);
}

//suma de vectores
double *sumV(double *v1, int n1, double *v2, int n2)
{
    double *result;
    int i, j;

    if (n1 != n2)
    {
        printf("sum: different number of rows or columns\n");
        exit(EXIT_FAILURE);
    }
    result = vector(n1);
    for (i = 0; i < n1; ++i)
    {
        result[i] = v1[i] + v2[i];
    }

    return result;
}

//mostrar el vector
void printVector(double *v, int n)
{
    int i;

    for (i = 0; i < n; i++)
        printf("%5.19lf ", v[i]);
    printf("\n");
}

//multiplicacion de un vector por un escalar
double *vec_x_esc(double *v, int n, double k)
{
    int i;
    double *L = (double *)calloc(n, sizeof(double));
    if (L == NULL)
        exit(EXIT_FAILURE);

    for (i = 0; i < n; i++)
        L[i] = v[i] * k;

    return L;
}

//comparacion de dos vectores
int compareV(double *v1, int n1, double *v2, int n2)
{
    int i, j;

    if (n1 != n2)
    {
        printf("vec: different dimension\n");
        return 0;
    }

    for (i = 0; i < n1; ++i)
        if (fabs(v1[i] - v2[i]) > pow(10, -3))
        {
            return 0;
        }

    return 1;
}

//crea un vector donde todas sus componentes son ceros
double *zerosV(int n)
{
    int i;
    double *L = (double *)calloc(n, sizeof(double));
    if (L == NULL)
        exit(EXIT_FAILURE);

    for (i = 0; i < n; i++)
    {
        L[i] = 0;
    }

    return L;
}

//matrices

//crea una matriz donde todas sus componentes son ceros
double **zeros(int n, int m)
{
    int i, dim = n * m;
    double **L = (double **)calloc(dim, sizeof(double));
    if (L == NULL)
        exit(EXIT_FAILURE);

    for (i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            L[i][j] = 0;
        }
    }

    return L;
}

//calculo de la matriz traspuesta
double **trasp(double **m, int n)
{
    double **result;
    int i, j;

    result = array(n, n);

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            result[i][j] = m[j][i];

    return result;
}

//asigna un vector de tamaño n como columna de una matriz en la posicion k
void asignar(double *m, int nf, int nc, int k, double *v, int n)
{
    int i;

    if (nf != n)
        exit(EXIT_FAILURE);

    for (i = 0; i <= n; i++)
        m[i * nc + k] = v[i];
}

// Calculo de la inversa usando el método de Gauss Jordan
double **inv(double **m, int n)
{
    double **mat, ratio, **result;
    int i, j, k;

    mat = array(n, 2 * n);
    result = array(n, n);

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            mat[i][j] = m[i][j];

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i == j)
            {
                mat[i][j + n] = 1;
            }
            else
            {
                mat[i][j + n] = 0;
            }
        }
    }

    /* Applying Gauss Jordan Elimination */

    for (i = 0; i < n; i++)
    {
        if (mat[i][i] == 0.0)
        {
            printf("inv: error\n");
            exit(EXIT_FAILURE);
        }
        for (j = 0; j < n; j++)
        {
            if (i != j)
            {
                ratio = mat[j][i] / mat[i][i];
                for (k = 0; k < 2 * n; k++)
                {
                    mat[j][k] = mat[j][k] - ratio * mat[i][k];
                }
            }
        }
    }

    /* Row Operation to Make Principal Diagonal to 1 */

    for (i = 0; i < n; i++)
    {
        for (j = n; j < 2 * n; j++)
        {
            mat[i][j] = mat[i][j] / mat[i][i];
        }
    }

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            result[i][j] = mat[i][j + n];

    freeArray(mat, n, 2 * n);

    return result;
}

//producto de matrices
double **prod(double **mat1, int nf1, int nc1, double **mat2, int nf2, int nc2)
{
    double **result;
    int i, j, k;
    //se comprueba el numero de filas de la primera matriz tenga el mismo num de columnas que la segunda
    if (nf1 != nc2)
    {
        printf("prod: nf1 != nc2\n");
        exit(EXIT_FAILURE);
    }
    result = array(nf1, nc2);
    for (i = 0; i < nf1; i++)
    {
        for (j = 0; j < nc2; j++)
        {
            result[i][j] = 0; //inicializar
            for (k = 0; k < nc1; k++)
            {
                result[i][j] = result[i][j] + mat1[i][k] * mat2[k][j];
            }
        }
    }

    return result;
}

//producto de una matriz con un vector

double *matXvec(double **mat, int nf, int nc, double *vec, int n)
{
    double *res;
    res = vector(n);

    int i;
    //se comprueba el numero de filas de la primera matriz tenga el mismo num de columnas que la segunda
    if (nf != n)
    {
        printf("problema de dimensiones\n");
        exit(EXIT_FAILURE);
    }

    double **mataux = array(nf,nc);

    for(i = 0; i < nf; i++)
    {
        mataux[i][0] = vec[i];
    }

    double **resmat = array(nf, nc);
    resmat = prod(mat, nf, nc, mataux, nf, nc);

    for(i = 0; i < nf; i++)
    {
        res[i] = resmat[i][0];
    }

    freeArray(mataux, nf, nc);
    freeArray(resmat, nf, nc);

    return res;
}

double *vecXmat(double *vec, int n, double **mat, int nf, int nc)
{
    double *res;
    res = vector(n);

    int i;
    //se comprueba el numero de filas de la primera matriz tenga el mismo num de columnas que la segunda
    if (nf != n)
    {
        printf("problema de dimensiones\n");
        exit(EXIT_FAILURE);
    }

    double **mataux = array(nf,nc);

    for(i = 0; i < nf; i++)
    {
        mataux[0][i] = vec[i];
    }

    double **resmat = array(nf, nc);
    resmat = prod(mataux, nf, nc, mat, nf, nc);

    for(i = 0; i < nf; i++)
    {
        res[i] = resmat[0][i];
    }

    freeArray(mataux, nf, nc);
    freeArray(resmat, nf, nc);

    return res;
}


//producto de un vector con otro vector traspuesto, el vector traspuesto tendrá solo 1 fila para que sea
//posible la multiplicacion, y tener mismo numero de columnas que filas del vector, se genera una matriz
double **vecXvecTrans(double *vec, int n, double *vecT, int nc)
{
    int i;
    //se comprueba el numero de filas de la primera matriz tenga el mismo num de columnas que la segunda
    if (nc != n)
    {
        printf("problema de dimensiones\n");
        exit(EXIT_FAILURE);
    }

    double **mataux = array(n,n);
    double **mataux2 = array(n,n);

    for(i = 0; i < n; i++)
    {
        mataux[i][0] = vec[i];
        mataux2[0][i] = vecT[i];
    }

   // printArray(mataux, 3, 3);
   // printf("\r\n");
   // printArray(mataux2, 3, 3);
   // printf("\r\n");

    double **resmat = array(nc, nc);
    resmat = prod(mataux, n, n, mataux2, n, n);

    freeArray(mataux, n, n);
    freeArray(mataux2, n, n);

    return resmat;
}

//multiplicacion de una matriz por un escalar
double **mat_x_esc(double **mat, int n, double k)
{
    int i, j;
    int dim = n * n;
    double **L = (double **)calloc(dim, sizeof(double));
    if (L == NULL)
        exit(EXIT_FAILURE);

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            L[i][j] = mat[i][j] * k;

    return L;
}

//matriz identidad
double **eye(int n)
{
    double **result;
    int i;

    result = array(n, n);
    for (i = 0; i < n; ++i)
        result[i][i] = 1.0;

    return result;
}

//suma de matrices
double **sum(double **mat1, int nf1, int nc1, double **mat2, int nf2, int nc2)
{
    double **result;
    int i, j;

    if ((nf1 != nf2) || (nc1 != nc2))
    {
        printf("sum: different number of rows or columns\n");
        exit(EXIT_FAILURE);
    }
    result = array(nf1, nc2);
    for (i = 0; i < nf1; ++i)
        for (j = 0; j < nc1; ++j)
        {
            result[i][j] = mat1[i][j] + mat2[i][j];
        }

    return result;
}

//comparacion de dos matrices
int compare(double **mat1, int nf1, int nc1, double **mat2, int nf2, int nc2)
{
    int i, j;

    if ((nf1 != nf2) || (nc1 != nc2))
    {
        printf("sum: different number of rows or columns\n");
        return 0;
    }

    for (i = 0; i < nf1; ++i)
        for (j = 0; j < nc1; ++j)
        {
            if (fabs(mat1[i][j] - mat2[i][j]) > pow(10, -3))
                return 0;
        }

    return 1;
}

//creacion de una matriz
double **array(int nf, int nc)
{
    double **m;
    int i;

    m = (double **)calloc(nf, sizeof(double *));
    if (m == NULL)
    {
        printf("array: memory not allocated\n");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < nf; i++)
    {
        m[i] = (double *)calloc(nc, sizeof(double));
        if (m[i] == NULL)
        {
            printf("array: memory not allocated\n");
            exit(EXIT_FAILURE);
        }
    }

    return m;
}

//liberacion de una matriz
void freeArray(double **mat, int nf, int nc)
{
    int i;

    for (i = 0; i < nf; i++)
        free(mat[i]);
    free(mat);
}

//se muestra la matriz
void printArray(double **mat, int nf, int nc)
{
    int i, j;

    for (i = 0; i < nf; i++)
    {
        for (j = 0; j < nc; j++)
        {
            printf("%5.19lf ", mat[i][j]);
        }
        printf("\n");
    }
}

//finds

int find1(double *v1, int n1, double *v2, int n2, double JD)
{
    int i1 = 0, i2 = 0, ind1 = 0, ind2 = 0;
    if (n1 != n2)
    {
        printf("find1: n1 != n2\n");
        exit(EXIT_FAILURE);
    }

    while (i1 < n1)
    {
        if (v1[i1] > JD)
        {
            ind1 = i1;
            break;
        }
        ++i1;
    }
    while (i2 < n2)
    {
        if (v2[i2] > JD)
        {
            ind2 = i2;
            break;
        }
        ++i2;
    }
    return (ind1 < ind2 ? ind1 : ind2);
}

//i = find(PC(:,1)<=JD & JD<=PC(:,2),1,'first');
int find2(double *v, int n, double JD)
{
    int i = 0;
    while (i < n)
    {
        if (fabs(v[i] - JD) < pow(10, -10))
        {
            break;
        }
        ++i;
    }
    return i;
}
