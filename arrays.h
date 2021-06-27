/** @file arrays.h
 *  @brief Specification of all the arrays and matrixes
 *  @author Victor Coman
 *  @date May 2021
 */

#ifndef ARRAY_h_
#define ARRAY_h_

// vectores

/** norm of a vector
 *  @param [in] v vector
 *  @param [in] n size of the vector
 *  @return a double
 */
double norma(double *v, int n);

/** Scalar product of 2 vectors
 *  @param [in] v1 vector 1
 *  @param [in] n1 size of the vector 1
 *  @param [in] v2 vector 2
 *  @param [in] n2 size of the vector 2
 *  @return a double
 */
double dot(double *v1, int n1, double *v2, int n2);

/** Creation of a vector
 *  @param [in] v vector
 *  @param [in] n size of the vector
 *  @return a double
 */
double *vector(int n);

/** Cross product of 2 vectors
 *  @param [in] vector_a vector 1
 *  @param [in] n1 size of the vector 1
 *  @param [in] vector_b vector 2
 *  @param [in] n2 size of the vector 2
 *  @return a vector
 */
double *crossProd(double *vector_a, int n1, double *vector_b, int n2);

/** Frees a vector
 *  @param [in] v vector
 *  @param [in] n size of the vector
 *  @return 0
 */
void freeVector(double *v, int n);

/** sum of 2 vectors
 *  @param [in] v1 vector 1
 *  @param [in] n1 size of the vector 1
 *  @param [in] v2 vector 2
 *  @param [in] n2 size of the vector 2
 *  @return a vector
 */
double *sumV(double *v1, int n1, double *v2, int n2);

/** product of a vector with a scalar
 *  @param [in] v vector 1
 *  @param [in] n size of the vector 1
 *  @param [in] k scalar number
 *  @return a vector
 */
double *vec_x_esc(double *v, int n, double k);

/** prints a vector
 *  @param [in] v vector
 *  @param [in] n size of the vector
 *  @return 0
 */
void printVector(double *v, int n);

/** creates a vector full of zeros
 *  @param [in] n size of the vector
 *  @return a vector
 */
double *zerosV(int n);


// matrices

/** creates a matrix full of zeros
 *  @param [in] n number of rows
 *  @param [in] m number of columns
 *  @return a matrix
 */
double **zeros(int n, int m);

/** repalce a column of the matrix with a vector
 *  @param [in] m matrix
 *  @param [in] nf number of rows
 *  @param [in] nc number of columns
 *  @param [in] k number of the column that will be changed
 *  @param [in] v vector
 *  @param [in] n size of the vector
 *  @return 0
 */
void asignar(double *m, int nf, int nc, int k, double *v, int n);

/** makes the transpose of a matrix
 *  @param [in] m matrix
 *  @param [in] n number of rows/columns
 *  @return a matrix
 */
double **trasp(double **m, int n);

/** makes the inverse of a matrix
 *  @param [in] m matrix
 *  @param [in] n number of rows/columns
 *  @return a matrix
 */
double **inv(double **m, int n);

/** product of 2 matrix
 *  @param [in] mat1 matrix 1
 *  @param [in] nf1 number of rows
 *  @param [in] nc1 number of columns
 *  @param [in] mat2 matrix 2
 *  @param [in] nf2 number of rows
 *  @param [in] nc2 number of columns
 *  @return a matrix
 */
double **prod(double **mat1, int nf1, int nc1, double **mat2, int nf2, int nc2);

/** product of a matrix and a vector
 *  @param [in] mat matrix
 *  @param [in] nf number of rows
 *  @param [in] nc number of columns
 *  @param [in] vec vector
 *  @param [in] n size of the vector
 *  @return a vector
 */
double *matXvec(double **mat, int nf, int nc, double *vec, int n);

/** product of a row vector and a matrix
 *  @param [in] mat matrix
 *  @param [in] nf number of rows
 *  @param [in] nc number of columns
 *  @param [in] vec vector
 *  @param [in] n size of the vector
 *  @return a vector
 */
double *vecXmat(double *vec, int n, double **mat, int nf, int nc);

/** product of a vector and a matrix
 *  @param [in] vec vector
 *  @param [in] n number of rows
 *  @param [in] vecT row vector
 *  @param [in] nc size of the vector
 *  @return a matrix
 */
double **vecXvecTrans(double *vec, int n, double *vecT, int nc);

/** product of a matrix and a scalar
 *  @param [in] mat matrix
 *  @param [in] n number of rows/columns
 *  @param [in] k scalar
 *  @return a matrix
 */
double **mat_x_esc(double **mat, int n, double k);

/** creates an identity matrix
 *  @param [in] n number of rows/columns
 *  @return a matrix
 */
double **eye(int n);

/** sum of 2 matrix
 *  @param [in] mat1 matrix 1
 *  @param [in] nf1 number of rows
 *  @param [in] nc1 number of columns
 *  @param [in] mat2 matrix 2
 *  @param [in] nf2 number of rows
 *  @param [in] nc2 number of columns
 *  @return a matrix
 */
double **sum(double **mat1, int nf1, int nc1, double **mat2, int nf2, int nc2);

/** compares 2 vectors to see if they are the same
 *  @param [in] v1 matrix 1
 *  @param [in] n1 number of rows
 *  @param [in] v2 matrix 2
 *  @param [in] n2 number of rows
 *  @return an integer, 0 if they are not and 1 if they are
 */
int compareV(double *v1, int n1, double *v2, int n2);

/** compares 2 matrix to see if they are the same
 *  @param [in] mat1 matrix 1
 *  @param [in] nf1 number of rows
 *  @param [in] nc1 number of columns
 *  @param [in] mat2 matrix 2
 *  @param [in] nf2 number of rows
 *  @param [in] nc2 number of columns
 *  @return an integer, 0 if they are not and 1 if they are
 */
int compare(double **mat1, int nf1, int nc1, double **mat2, int nf2, int nc2);

/** creates a matrix
 *  @param [in] nf number of rows
 *  @param [in] nc number of columns
 *  @return a matrix
 */
double **array(int nf, int nc);

/** frees a matrix
 *  @param [in] mat matrix
 *  @param [in] nf number of rows
 *  @param [in] nc number of columns
 *  @return 0
 */
void freeArray(double **mat, int nf, int nc);

/** prints a matrix
 *  @param [in] mat matrix
 *  @param [in] nf number of rows
 *  @param [in] nc number of columns
 *  @return 0
 */
void printArray(double **mat, int nf, int nc);


//finds

/** finds JD
 *  @param [in] v1 vector 1
 *  @param [in] n1 size
 *  @param [in] v2 vector 2
 *  @param [in] n2 size
 *  @param [in] JD what we want to find
 */
int find1(double *v1, int n1, double *v2, int n2, double JD);


/** finds JD
 *  @param [in] v1 vector 1
 *  @param [in] n1 size
 *  @param [in] JD what we want to find
 */
int find2(double *v, int n, double JD);

#endif /* ARRAY_h_ */
