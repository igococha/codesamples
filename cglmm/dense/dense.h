#ifndef DENSE_H
#define DENSE_H

#include <stdio.h>
#include <stddef.h>

#include "pm_rng.h"

/* **** Structures and datatypes **** */

typedef struct DenseBlock DenseBlock;

struct DenseBlock {
  double* data;
  size_t size;
};

typedef enum { DENSE_ROW_MAJOR, DENSE_COLUMN_MAJOR } DenseOrder;
typedef enum { DENSE_TRANS, DENSE_TRANS_NO } DenseTrans;
typedef enum { DENSE_UP, DENSE_LO } DenseUPLO;

typedef struct {
  size_t rows;
  size_t columns;
  size_t tda;
  double* data;
  DenseBlock* block;
  DenseOrder order;
  int owner;
} Dense;


typedef struct {
  size_t size;
  size_t stride;
  double* data;
  DenseBlock* block;
  int owner;
} DenseVector;

/* ****** end structures *** */


/* ***** Dense Block interface ******* */

DenseBlock* dense_block_new(size_t size);
void dense_block_free(DenseBlock** b);


/* ******** Dense matrix Interface ********** */


#define DENSE_GET_DATA(M) M->data
#define DENSE_GET_ROWS(M) M->rows
#define DENSE_GET_COLUMNS(M) M->columns
#define DENSE_GET_TDA(M) M->tda
#define DENSE_SET(M,I,J,V) *((M->order==DENSE_ROW_MAJOR)?\
  &(M->data[I*M->tda+J]) : &(M->data[J*M->tda+I])) = V
#define DENSE_GET(M,I,J) (M->order==DENSE_ROW_MAJOR)?\
  (M->data[I*M->tda+J]) : (M->data[J*M->tda+I])
#define DENSE_GET_RMAJOR(M,I,J) (M->data[I*M->tda+J])
#define DENSE_GET_CMAJOR(M,I,J) (M->data[J*M->tda+I]])

Dense* dense_new(size_t rows, size_t cols);
Dense* dense_new_order(size_t rows, size_t cols, DenseOrder order);
Dense* dense_new_sub(Dense* m, size_t row, size_t col, size_t rows, size_t columns);
void dense_free(Dense** m);
void dense_set_order(Dense *m, DenseOrder order);
void dense_print(Dense* m, FILE* fp);

int dense_copy(Dense *t, Dense *s);
int dense_copy_trans(Dense *t, Dense *s);
void dense_set_all(Dense* m, double v);
void dense_set_diag(Dense* m, double v);
int dense_set_column(Dense* m, int col_idx, double v);
int dense_set_column_intarray(Dense* m, int col_idx, int* data);
int dense_set_column_doublearray(Dense* m, int col_idx, double* data);
int dense_set_column_range(Dense* m, int col_idx, double v,int from,int to);
int dense_set_column_intarray_range(Dense* m, int col_idx, int* data, int from, int to);
int dense_set_column_doublearray_range(Dense* m, int col_idx, double* data, int from, int to);

/* B = alpha*A + beta*B  */
int dense_add(Dense *B, double alpha, DenseTrans transA, Dense *A, double beta);

/* BLAS 3 */

/* DGEMM. C = alpha*op(A)*op(B) + beta*C */
void dense_mbym(Dense* C, double alpha, DenseTrans transA, Dense *A, DenseTrans transB, Dense *B, double beta);
int dense_mtbym(Dense *A, Dense *B);

void dense_mbyc(Dense* m, double alpha);
void dense_mbyc_diag(Dense *A, double c);
void dense_mbyc_nondiag(Dense *A, double c);

/* LAPACK */
int dense_cholesky(Dense *m, DenseUPLO uplo);
int dense_inverse_cholesky(Dense *m, DenseUPLO uplo);
int dense_solve_triangular(Dense* A,DenseUPLO uplo, DenseTrans trans, Dense* b);
int dense_solve_positive(Dense *A, DenseUPLO uplo, Dense* b);
int dense_solve_cholesky(Dense *A, DenseUPLO uplo, Dense* b);


/* Function i.e. f(A)  */
int dense_diag_max(Dense *A, double *max);
int dense_triangular_logdet(Dense *A, double *logdet);

void dense_upd_mult_all(Dense* m, double alpha);
void dense_mult_diag(Dense *A, double c);
void dense_mult_nondiag(Dense *A, double c);

/* ************** Dense Vector ************ */

#define DENSE_VECTOR_SET(V,I,C) (V->data[V->stride * I]=C)
#define DENSE_VECTOR_GET(V,I) (V->data[V->stride * I])

DenseVector* dense_vector_new(size_t size);
DenseVector* dense_vector_new_sub(DenseVector *v,size_t i, size_t size);
DenseVector* dense_vector_matrix_row(Dense* matrix,int row);
DenseVector* dense_vector_matrix_column(Dense* matrix,int col);
void dense_vector_free(DenseVector** vector);
void dense_vector_print(DenseVector* vector,FILE* fp);

void dense_vector_set_all(DenseVector* vector, double val);
void dense_vector_set_range(DenseVector* vector, double val,int from,int to);
void dense_vector_set_intarray(DenseVector* vector, int* data);
void dense_vector_set_doublearray(DenseVector* vector, double* data);
void dense_vector_set_intarray_range(DenseVector* vector, int* data, int from, int to);
void dense_vector_set_doublearray_range(DenseVector* vector, double* data, int from, int to);
int dense_vector_copy_idx(DenseVector *vdest, size_t fd, DenseVector *vsource, size_t fs, size_t n);


/* Random */
void dense_vector_random(DenseVector* vector, pm_rng* rng);
void dense_vector_ran_gaussian(DenseVector* vector, pm_rng* rng,double mean, double sigma);



/* BLAS 1 */
void dense_vector_scale(DenseVector* vector, double alpha);
void dense_vector_copy(DenseVector *vdest, DenseVector *vsource);
void dense_vector_add(DenseVector *vdest, double alpha, DenseVector* vsource);
double dense_vector_dot(DenseVector* v1, DenseVector *v2);
double dense_vector_asum(DenseVector* vector);

/* BLAS 2 */
/* DGEMV - y <- alpha*A*x + beta*y */
void dense_vector_mbyv_inc(DenseVector* y, double alpha, DenseTrans transA, Dense* A, DenseVector* x, double beta);
double dense_vector_vtbyv(DenseVector *v);

/* LAPACK */
int dense_vector_solve_triangular(Dense *A, DenseUPLO uplo, DenseTrans trans, DenseVector* b);
int dense_vector_solve_positive(Dense *A, DenseUPLO uplo, DenseVector* b);
int dense_vector_solve_cholesky(Dense *A, DenseUPLO uplo, DenseVector* b);





#endif
