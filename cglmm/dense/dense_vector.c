#include "dense.h"
#include "cblas.h"
#include "lapacke.h"
#include "pm_error.h"
#include "pm_randist.h"

#include <stdlib.h>
#include <assert.h>

DenseVector* dense_vector_new(size_t size) {
  DenseVector* vector;
  DenseBlock* block;
  block = dense_block_new(size);
  vector = (DenseVector*)malloc(sizeof(DenseVector));
  vector->size = size;
  vector->stride = 1;
  vector->data = block->data;
  vector->block = block;
  vector->owner = 1;

  return vector;
}

/* These version reuse the matrix's block */
/* Deallocation of the matrix will result in invalid pointers */
DenseVector* dense_vector_matrix_row(Dense* matrix,int row) {
  DenseVector* vector;

  vector = (DenseVector*)malloc(sizeof(DenseVector));
  vector->block = matrix->block;
  vector->size = matrix->columns;
  vector->owner = 0;
  if (matrix->order == DENSE_ROW_MAJOR) {
    vector->stride = 1;
    vector->data = matrix->data + row*(matrix->tda);
  } else { /* Column major */
    vector->stride = matrix->tda;
    vector->data = matrix->data + row;
  }
  return vector;
}

/* Creates sub-vector  - re-uses block */
DenseVector* dense_vector_new_sub(DenseVector *v,size_t i, size_t size) {
  DenseVector *sub;

  assert((v!=NULL)&&(i<v->size)&&(i+size <= v->size));
  
  sub = (DenseVector*)malloc(sizeof(DenseVector));
  if (sub==NULL) return NULL;
  sub->block = v->block;
  sub->owner = 0;
  sub->size = size;
  sub->stride = v->stride;
  sub->data = v->data + i*(v->stride);
  return sub;
}

DenseVector* dense_vector_matrix_column(Dense* matrix,int col) {
  DenseVector* vector;

  vector = (DenseVector*)malloc(sizeof(DenseVector));
  vector->block = matrix->block;
  vector->size = matrix->rows;
  vector->owner = 0;
  if (matrix->order == DENSE_COLUMN_MAJOR) {
    vector->stride = 1;
    vector->data = matrix->data + col*(matrix->tda);
  } else { /* Column major */
    vector->stride = matrix->tda;
    vector->data = matrix->data + col;
  }
  return vector;
}

void dense_vector_free(DenseVector** vector) {
  if (*vector==NULL) return;
  if ((*vector)->owner) {
    dense_block_free(&(*vector)->block);
  }
  free(*vector);
  *vector = NULL;
  return;
}


void dense_vector_print(DenseVector* vector,FILE* fp) {
  int i;
  double* p;
  p = vector->data;
  fprintf(fp, "[ ");
  for(i=0; i < vector->size; i++) {
    fprintf(fp,"%f ",*p);
    p += vector->stride;
  }
  fprintf(fp,"]");
  return;
}

void dense_vector_set_all(DenseVector* vector, double val) {
  int i;
  double *p = vector->data;
  for(i=0; i < vector->size; i++) {
    *p = val;
    p += vector->stride;
  }
}
void dense_vector_set_range(DenseVector* vector, double val, int from, int to) {
  int i;
  double *p = vector->data + vector->stride*from;
  for(i=from; i <= to; i++) {
    *p = val;
    p += vector->stride;
  }
}

void dense_vector_set_intarray(DenseVector* vector, int* data) {
  int i;
  double *p = vector->data;
  for(i=0; i < vector->size; i++,data++) {
    *p = (double)*data;
    p += vector->stride;
  }
}

void dense_vector_set_doublearray(DenseVector* vector, double* data) {
  int i;
  double *p = vector->data;
  for(i=0; i < vector->size; i++,data++) {
    *p = *data;
    p += vector->stride;
  }
}

void dense_vector_set_intarray_range(DenseVector* vector, int* data, int from, int to) {
  int i;
  double *p = vector->data + vector->stride*from;
  for(i=from; i <= to; i++,data++) {
    *p = (double)*data;
    p += vector->stride;
  }
}

void dense_vector_set_doublearray_range(DenseVector* vector, double* data, int from, int to) {
  int i;
  double *p = vector->data + vector->stride*from;
  for(i=from; i <= to; i++,data++) {
    *p = *data;
    p += vector->stride;
  }
}

int dense_vector_copy_idx(DenseVector *vdest, size_t fd, DenseVector *vsource, size_t fs, size_t n) {
  double *p, *q;
  int i;
  size_t p_stride, q_stride;

  if ((fd+n > vdest->size)&&(fs+n > vsource->size)){
    PM_ERROR(PM_EOUTOFBOUNDS,"Vector indices out of bounds");
  }
  
  p = &vdest->data[fd];
  q = &vsource->data[fs];
  p_stride = vdest->stride;
  q_stride = vsource->stride;
  for(i=0; i < n; i++) {
    *p = *q;
    p += p_stride;
    q += q_stride;
  }
  return 0;
}


void dense_vector_random(DenseVector* vector, pm_rng* rng) {
  pm_rng_uniform_seq(rng,vector->data,vector->size,vector->stride);
}

void dense_vector_ran_gaussian(DenseVector* vector, pm_rng* rng,double mean, double sigma) {
  int i;
  double *data;

  data = vector->data;
  /* wasteful since box-muller generates two gaussians at a time 
     and the rejection method might prevent optimizations (vectorization) */
  if (mean==0.0) {
    for(i=0; i < vector->size; i++) {
      *data = pm_ran_gaussian(rng,sigma);
      data += vector->stride;
    }
  } else {
    for(i=0; i < vector->size; i++) {
      *data = mean + pm_ran_gaussian(rng,sigma);
      data += vector->stride;
    }
  }


}


/* BLAS 1 */

void dense_vector_scale(DenseVector* vector, double alpha) {
  cblas_dscal((int)vector->size,alpha,vector->data,(int)vector->stride);
  return;
}

void dense_vector_copy(DenseVector *vdest, DenseVector *vsource) {
  int size;
  size = (vdest->size<vsource->size?vdest->size:vsource->size);
  cblas_dcopy((int)size,vsource->data,(int)vsource->stride,vdest->data,(int)vdest->stride);
}


void dense_vector_add(DenseVector *vdest, double alpha, DenseVector* vsource) {
  int size;
  size = (vdest->size<vsource->size?vdest->size:vsource->size);
  cblas_daxpy((int)size,alpha,vsource->data,(int)vsource->stride,vdest->data,(int)vdest->stride);
}

double dense_vector_dot(DenseVector* v1, DenseVector *v2) {
  int size;
  double d;
  size = (v1->size < v2->size ? v1->size:v2->size);
  d =  cblas_ddot(size,v1->data,(int)v1->stride,v2->data,(int)v2->stride);
  return d;
}


double dense_vector_asum(DenseVector *vector) {
  return cblas_dasum((int)vector->size,vector->data,(int)vector->stride);
}


/* Blas 2 */


void dense_vector_mbyv_inc(DenseVector* y, double alpha, DenseTrans transA, Dense* A, DenseVector* x, double beta) {
  enum CBLAS_ORDER order;
  enum CBLAS_TRANSPOSE tA;

  order = (A->order == DENSE_COLUMN_MAJOR)? CblasColMajor:CblasRowMajor;
  tA = (transA==DENSE_TRANS) ? CblasTrans : CblasNoTrans;

  /* a row major matrix will require extra memory for workspace */ 
  cblas_dgemv(order,tA,A->rows,A->columns,alpha,A->data,A->tda,x->data,x->stride,beta,y->data,y->stride);

  return;

}

double dense_vector_vtbyv(DenseVector *v) {
  double c,prod,*p;
  int i;
  p = v->data;
  prod = 0;
  for(i=0; i < v->size;i++) {
    c = *p;
    prod += c*c;
    p += v->stride;
  }
  return prod;
}

int dense_vector_solve_triangular(Dense* A,DenseUPLO uplo, DenseTrans trans, DenseVector* b) {
  lapack_int status;
  enum CBLAS_ORDER order;
  char ul,tr;

  /* this may only work with stride = 1 */
  /* temporary memory should be used for this */

  order = (A->order == DENSE_COLUMN_MAJOR)? CblasColMajor:CblasRowMajor;  
  ul = (uplo==DENSE_LO)?'L':'U';
  tr = (trans==DENSE_TRANS)?'T':'N';
  
  status = LAPACKE_dtrtrs(order,ul,tr,'N',A->rows,1,A->data,A->tda,b->data,b->size);


  return PM_SUCCESS;

}

int dense_vector_solve_positive(Dense *A, DenseUPLO uplo, DenseVector* b) {
  lapack_int status;
  enum CBLAS_ORDER order;
  char ul;

  /* todo: check A and b have same order */
  
  order = (A->order == DENSE_COLUMN_MAJOR)? CblasColMajor:CblasRowMajor;  
  ul = (uplo==DENSE_LO)?'L':'U';
 
  status = LAPACKE_dposv(order,ul,A->rows,1,A->data,A->tda,b->data,b->size);
  
  return PM_SUCCESS;
}

int dense_vector_solve_cholesky(Dense *A, DenseUPLO uplo, DenseVector* b) {
  lapack_int status;
  enum CBLAS_ORDER order;
  char ul;

  /* todo: check A and b have same order */
  
  order = (A->order == DENSE_COLUMN_MAJOR)? CblasColMajor:CblasRowMajor;  
  ul = (uplo==DENSE_LO)?'L':'U';
 
  status = LAPACKE_dpotrs(order,ul,A->rows,1,A->data,A->tda,b->data,b->size);
  
  return PM_SUCCESS;
}

