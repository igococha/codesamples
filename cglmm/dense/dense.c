#include "dense.h"
#include "cblas.h"
#include "lapacke.h"
#include "pm_error.h"

#include <stdlib.h>
#include <math.h>
#include <assert.h>

/* DenseBlock implementation */


DenseBlock* dense_block_new(size_t size) {
  DenseBlock* block = (DenseBlock*)malloc(sizeof(DenseBlock));
  double* data = (double*)malloc(size * sizeof(double));
  block->data = data;
  block->size = size;
  return block;
}

void dense_block_free(DenseBlock** b) {
  free((*b)->data);
  free(*b);
  return;
}

/* Dense implementation */
Dense* dense_new(size_t rows, size_t columns) {
  return dense_new_order(rows,columns,DENSE_COLUMN_MAJOR);
}

Dense* dense_new_order(size_t rows, size_t columns, DenseOrder order) {
  Dense* m;
  DenseBlock* block;
  m = (Dense*)malloc(sizeof(Dense));
  m->rows = rows;
  m->columns = columns;
  m->tda = (order==DENSE_ROW_MAJOR)?columns:rows;
  m->block = dense_block_new(rows*columns);
  m->data = m->block->data;
  m->order = order;
  m->owner = 1;
  return m;
}

/* Creates Submatrix: a matrix that shares the same block, order and tda  */
Dense* dense_new_sub(Dense* m, size_t row, size_t col, size_t rows, size_t columns) {
  Dense* sub;

  assert((m!=NULL)&&(row<m->rows)&&(col<m->columns));
  assert((row+rows <= m->rows)&&(col+columns <= m->columns));
  
  sub = (Dense*)malloc(sizeof(Dense));
  if (sub==NULL) return NULL;
  sub->rows = rows;
  sub->columns = columns;  
  sub->order = m->order;
  sub->tda = m->tda;
  sub->block = m->block;
  sub->owner = 0;
  if (m->order == DENSE_COLUMN_MAJOR) {
    sub->data = m->data + col*m->tda + row;
  } else {
    sub->data = m->data + row*m->tda + col;
  }
  return sub;
}

void dense_free(Dense** m) {
  if (*m == NULL) return;
  if ((*m)->owner) {
    dense_block_free(&(*m)->block);
  }
  free(*m);
  *m=  NULL;
  return;
}

/* 
Changes the matrix's order. It does not reorder the data .
Valid only for square matrices (so the tda remains valid)
The end result for square matrices is the matrix's transpose with the order swapped.
*/
void dense_set_order(Dense* m, DenseOrder order) {
  if (m->order != order) {
    if (m->rows==m->columns) {
      m->order = order;
    }
  }
}

void dense_print(Dense* m, FILE* fp) {
  int i,j,ni, nj, tda, step, next;
  double* pbuffer, *pstart;
  tda = m->tda;
  pstart = m->data;
  if (m->order == DENSE_ROW_MAJOR) {
    step=1; next=tda;
  } else {
    step=tda; next=1;
  }

  for(i=0; i < m->rows; i++) {
    pbuffer = pstart;
    for(j=0; j < m->columns; j++) {
      fprintf(fp,"%.5f ",*pbuffer);
      pbuffer+=step;
    }
    fprintf(fp,"\n");
    pstart += next;
  }
  return;
}

/* different orderings allowed */
int dense_copy(Dense *t, Dense *s) {
  int ttda,sstep,stda,i,j,ni,nj;
  double *tdata,*tstart,*sdata,*sstart;

  pm_error_clear(NULL);
  if ((t->rows != s->rows)||(t->columns != s->columns)) {
    PM_ERROR_NOMSG(PM_EMATRIXSIZE);
  }
  if (t->order==DENSE_COLUMN_MAJOR) {
    ni = t->columns;
    nj = t->rows;
  } else {
    ni = t->rows;
    nj = t->columns;
  }
  ttda =t->tda;
  
  if (t->order == s->order) {
    sstep = 1; stda = s->tda;
  } else {
    sstep = s->tda; stda = 1;
  }

  tstart = t->data;
  sstart = s->data;
  for(i=0;i < ni; i++) {
    tdata = tstart;  sdata = sstart;
    for(j=0;j < nj; j++) {
      *(tdata++) = *sdata;
      sdata += sstep;
    }
    tstart += ttda;
    sstart += stda;
  }
  return PM_SUCCESS;
}

/* Copy transpose */
/* different orderings allowed */
int dense_copy_trans(Dense *t, Dense *s) {
  int ttda,sstep,stda,i,j,ni,nj;
  double *tdata,*tstart,*sdata,*sstart;

  pm_error_clear(NULL);
  if ((t->rows != s->columns)||(t->columns != s->rows)) {
    PM_ERROR_NOMSG(PM_EMATRIXSIZE);
  }
  if (t->order==DENSE_COLUMN_MAJOR) {
    ni = t->columns;
    nj = t->rows;
  } else {
    ni = t->rows;
    nj = t->columns;
  }
  ttda =t->tda;
  
  if (t->order != s->order) {
    sstep = 1; stda = s->tda;
  } else {
    sstep = s->tda; stda = 1;
  }

  tstart = t->data;
  sstart = s->data;
  for(i=0;i < ni; i++) {
    tdata = tstart;  sdata = sstart;
    for(j=0;j < nj; j++) {
      *(tdata++) = *sdata;
      sdata += sstep;
    }
    tstart += ttda;
    sstart += stda;
  }
  return PM_SUCCESS;
}

/* copies lower (upper) triangular section to opposite corner */
int dense_copy_triangular(Dense *m, DenseUPLO uplo) {
  size_t i,j;
  double *f,*fstart,*fdata,*t,*tstart,*tdata;
  size_t fstep,fjump,tstep,tjump;

  if (m==NULL) {
    PM_ERROR(PM_ENP,"Argument to function is a null matrix");
  }
  if (m->rows != m->columns) {
    PM_ERROR(PM_EARGUMENTS,"Operation valind only on square matrices");
  }

  if (((m->order==DENSE_COLUMN_MAJOR)&&(uplo==DENSE_LO)) ||
      ((m->order==DENSE_ROW_MAJOR)&&(uplo==DENSE_UP))) {
    fstep = 1; fjump = m->tda;
  } else {
    fjump = 1; fstep = m->tda;
  }
  tstep = fjump;
  tjump = fstep;
  
  fstart = m->data+fstep;  /* skip diagonal */
  tstart = m->data+tstep;  /* skip diagonal */  
  for(i=1; i < m->rows; i++) {
    fdata = fstart;
    tdata = tstart;
    for(j=0; j < m->rows-i; j++) {
      *tdata = *fdata;
      fdata += fstep;
      tdata += tstep;
    }
    fstart += (fjump+fstep);
    tstart += (tjump+tstep);
  }
  
  return PM_SUCCESS;
  
}

void dense_set_all(Dense* m, double v) {
  int i,j,ni, nj, tda,step,next;
  double* pbuffer, *pstart;
  tda = m->tda;
  pstart = m->data;
  if (m->order == DENSE_ROW_MAJOR) {
    ni=m->rows; nj=m->columns; step=1; next=tda;
  } else {
    ni=m->columns; nj=m->rows; step=tda; next=1;
  }

  for(i=0; i < ni; i++) {
    pbuffer = pstart;
    for(j=0; j < nj; j++) {
      *pbuffer=v;
      pbuffer++;
    }    
    pstart += tda;
  }
  return;
}

void dense_set_diag(Dense* A, double v) {
  int i,min;
  double *pdiag;
  
  min = (A->rows < A->columns)?A->rows:A->columns;

  pdiag = A->data;
  for(i=0; i<min; i++) {
    *pdiag = v;
    pdiag += A->tda+1;
  }
  return; 
}


int dense_set_column(Dense* m, int col_idx, double v) {
  double *p;
  int i,step,jump;

  /* error testing for col_idx? */
  
  if (m->order == DENSE_ROW_MAJOR) {
    step = m->tda;
    jump = 1;
  } else { /* Column Major */
    step = 1;
    jump = m->tda;
  }  
  p = m->data+col_idx*jump;
  for(i= 0; i < m->rows; i++) {
    *p = v;
    p += step;
  }  
  return 0;
}

int dense_set_column_intarray(Dense* m, int col_idx, int* data) {
  double *p;
  int i,step,jump;
  int *v;

  /* error testing for col_idx? */
  
  if (m->order == DENSE_ROW_MAJOR) {
    step = m->tda;
    jump = 1;
  } else { /* Column Major */
    step = 1;
    jump = m->tda;
  }
  v = data;
  p = m->data+col_idx*jump;
  for(i= 0; i < m->rows; i++,v++) {
    *p = (double)*v;
    p += step;
  }  
  return 0;
}


int dense_set_column_doublearray(Dense* m, int col_idx, double* data) {
  double *p;
  int i,step,jump;
  double *v;

  /* error testing for col_idx? */
  
  if (m->order == DENSE_ROW_MAJOR) {
    step = m->tda;
    jump = 1;
  } else { /* Column Major */
    step = 1;
    jump = m->tda;
  }
  v = data;
  p = m->data+col_idx*jump;
  for(i= 0; i < m->rows; i++,v++) {
    *p = *v;
    p += step;
  }  
  return 0;
}

int dense_set_column_intarray_range(Dense* m, int col_idx, int* data,int from, int to) {
  double *p;
  int i,step,jump;
  int *v;

  /* error testing for col_idx? */
  
  if (m->order == DENSE_ROW_MAJOR) {
    step = m->tda;
    jump = 1;
  } else { /* Column Major */
    step = 1;
    jump = m->tda;
  }
  v = data;
  p = m->data+col_idx*jump + step*from;
  for(i=from; i <= to; i++,v++) {
    *p = (double)*v;
    p += step;
  }  
  return 0;
}


int dense_set_column_range(Dense* m, int col_idx, double v,int from,int to) {
  double *p;
  int i,step,jump;

  /* error testing for col_idx? */
  
  if (m->order == DENSE_ROW_MAJOR) {
    step = m->tda;
    jump = 1;
  } else { /* Column Major */
    step = 1;
    jump = m->tda;
  }  
  p = m->data + col_idx*jump + from*step;
  for(i= from; i <= to; i++) {
    *p = v;
    p += step;
  }  
  return 0;
}

int dense_set_column_doublearray_range(Dense* m, int col_idx, double* data,int from, int to) {
  double *p;
  int i,step,jump;
  double *v;

  /* error testing for col_idx? */
  
  if (m->order == DENSE_ROW_MAJOR) {
    step = m->tda;
    jump = 1;
  } else { /* Column Major */
    step = 1;
    jump = m->tda;
  }
  v = data;
  p = m->data+col_idx*jump + step*from;
  for(i=from; i <= to; i++,v++) {
    *p = *v;
    p += step;
  }  
  return 0;
}

/* B = alpha*A + beta*B  */
int dense_add(Dense *B, double alpha, DenseTrans transA, Dense *A, double beta) {
  int Btda,Astep,Atda,i,j,ni,nj;
  double *Bdata, *Bstart, *Adata, *Astart,c;

  pm_error_clear(NULL);
  if (transA == DENSE_TRANS_NO) {
    if ((B->rows != A->rows)||(B->columns != A->columns)) {
      PM_ERROR_NOMSG(PM_EMATRIXSIZE);
    }
  } else {
    if ((B->rows != A->columns)||(B->columns != A->rows)) {
      PM_ERROR_NOMSG(PM_EMATRIXSIZE);
    }
  }

  if (B->order==DENSE_COLUMN_MAJOR) {
    ni = B->columns;
    nj = B->rows;
  } else {
    ni = B->rows;
    nj = B->columns;
  }
  Btda =B->tda;
  
  if (( B->order == A->order) && (transA == DENSE_TRANS_NO)) {
    Astep = 1; Atda = A->tda;
  } else {
    Astep = A->tda; Atda = 1;
  }

  Bstart = B->data;
  Astart = A->data;
 
  if (beta != 0) {
    for(i=0;i < ni; i++) {
      Bdata = Bstart;  Adata = Astart;
      for(j=0;j < nj; j++) {
	c = beta*(*Bdata);
	*(Bdata++) = alpha*(*Adata) + c;
	Adata += Astep;
      }
      Bstart += Btda;
      Astart += Atda;
    }
  } else { /* beta == 0 */
    for(i=0;i < ni; i++) {
      Bdata = Bstart;  Adata = Astart;
      for(j=0;j < nj; j++) {
	*(Bdata++) = alpha*(*Adata);
	Adata += Astep;
      }
      Bstart += Btda;
      Astart += Atda;
    }
  }
  return PM_SUCCESS;


  return 0;
    }


void dense_mbyc(Dense* m, double alpha) {
  int i,j,ni, nj, tda,step,next;
  double* pbuffer, *pstart;
  tda = m->tda;
  pstart = m->data;
  if (m->order == DENSE_ROW_MAJOR) {
    ni=m->rows; nj=m->columns; step=1; next=tda;
  } else {
    ni=m->columns; nj=m->rows; step=tda; next=1;
  }

  for(i=0; i < ni; i++) {
    pbuffer = pstart;
    for(j=0; j < nj; j++) {
      *pbuffer *= alpha;
      pbuffer++;
    }    
    pstart += tda;
  }
  return;

}

void dense_mbyc_diag(Dense *A, double alpha) {
  int i,min;
  double *pdiag;
  
  min = (A->rows < A->columns)?A->rows:A->columns;

  pdiag = A->data;
  for(i=0; i<min; i++) {
    *pdiag *= alpha;
    pdiag += A->tda+1;
  }
  return;
}

void dense_mbyc_nondiag(Dense *A, double alpha) {
  int i,j,ni, nj, tda,step,next,min;
  double* pbuffer, *pstart;
  
  tda = A->tda;
  pstart = A->data;
  if (A->order == DENSE_ROW_MAJOR) {
    ni=A->rows; nj=A->columns; step=1; next=tda;
  } else {
    ni=A->columns; nj=A->rows; step=tda; next=1;
  }
  min = (ni<nj)?ni:nj;
  for(i=0; i < ni; i++) {
    pbuffer = pstart;
    /* before diagonal */
    for(j=0; j < i; j++) {
      *pbuffer *= alpha;
      pbuffer++;
    }
    /* skip j=i */
    pbuffer++;
    /* after diagonal */
    for(j=i+1; j < min; j++) {
      *pbuffer *= alpha;
      pbuffer++;
    }
    /* update row/column start */
    pstart += tda;
  }
  return;
}


/* BLAS 3 */


void dense_mbym(Dense* C, double alpha, DenseTrans transA, Dense *A, DenseTrans transB, Dense *B, double beta) {
  enum CBLAS_ORDER order;
  enum CBLAS_TRANSPOSE tA,tB;
  int k;

  if (A->order != B->order) return;
  
  order = (A->order == DENSE_COLUMN_MAJOR)? CblasColMajor:CblasRowMajor;
  
  tA = (transA==DENSE_TRANS) ? CblasTrans : CblasNoTrans;
  tB = (transB==DENSE_TRANS) ? CblasTrans : CblasNoTrans;
  k = (transA==DENSE_TRANS) ? A->rows : A->columns;

  cblas_dgemm(order,tA,tB,C->rows,C->columns,k,alpha,A->data,A->tda,B->data,B->tda,beta,C->data, C->tda );
  

}

/* X = A^T * A */
int dense_mtbym(Dense *X, Dense *A) {
  int m,n,i,j,k,step,tda,stepx,tdax;
  double *p1,*p2,*p1start,*p2start,*px,*pxstart,*pxtr;
  double result;

  m = A->rows;
  n = A->columns;

  /* X must be n x n */
  if ((X->rows != X->columns)||(X->rows!=n)) {
    PM_ERROR_NOMSG(PM_EMATRIXSIZE);
  }
  step = 1; tda = A->tda;
  stepx = 1; tdax = X->tda;
  
  p1start = A->data;
  pxstart = X->data; /* diagonal is written twice */
  for(i=0; i < n; i++) {
    p2start = p1start;
    px = pxstart;
    pxtr = pxstart;
    for(j=i; j < n; j++) {
      /* dot product */
      result = 0.0;
      p1 = p1start;
      p2 = p2start;
      for(k=0; k < m; k++) {
	result += (*p1)*(*p2);
	p1 += step;
	p2 += step;
      }
      *px = result;
      *pxtr = result;
      px += stepx;
      pxtr += tdax; /* transpose */
      p2start += tda;
    }
    p1start += tda; /* next column */
    pxstart += tdax+stepx; /* next diagonal */
  }
  
  return PM_SUCCESS;
}

/*  LAPACK */

/*lapack_int LAPACKE_dgetrf( int matrix_order, lapack_int m, lapack_int n,
                           double* a, lapack_int lda, lapack_int* ipiv );
*/

int dense_cholesky(Dense *m, DenseUPLO uplo) {
  enum CBLAS_ORDER order;
  lapack_int status;
  char ul;

  order = (m->order == DENSE_COLUMN_MAJOR)? CblasColMajor:CblasRowMajor;  
  ul = (uplo==DENSE_LO)?'L':'U';
  
  status = LAPACKE_dpotrf(order,ul,m->rows,m->data,m->tda);
  
  return PM_SUCCESS;
}

int dense_inverse_cholesky(Dense *m, DenseUPLO uplo) {
  enum CBLAS_ORDER order;
  lapack_int status;
  char ul;

  order = (m->order == DENSE_COLUMN_MAJOR)? CblasColMajor:CblasRowMajor;  
  ul = (uplo==DENSE_LO)?'L':'U';
  
  status = LAPACKE_dpotri(order,ul,m->rows,m->data,m->tda);

  dense_copy_triangular(m,uplo);
  
  return PM_SUCCESS;
}


int dense_solve_triangular(Dense* A,DenseUPLO uplo, DenseTrans trans, Dense* b) {
  lapack_int status;
  enum CBLAS_ORDER order;
  char ul, tr;

  order = (A->order == DENSE_COLUMN_MAJOR)? CblasColMajor:CblasRowMajor;  
  ul = (uplo==DENSE_LO)?'L':'U';
  tr = (trans==DENSE_TRANS)?'T':'N';
  
  status = LAPACKE_dtrtrs(order,ul,tr,'N',A->rows,b->columns,A->data,A->tda,b->data,b->tda );


  return PM_SUCCESS;

}

int dense_solve_positive(Dense *A, DenseUPLO uplo, Dense* b) {
  lapack_int status;
  enum CBLAS_ORDER order;
  char ul;

  /* todo: check A and b have same order */
  
  order = (A->order == DENSE_COLUMN_MAJOR)? CblasColMajor:CblasRowMajor;  
  ul = (uplo==DENSE_LO)?'L':'U';
 
  status = LAPACKE_dposv(order,ul,A->rows,b->columns,A->data,A->tda,b->data,b->tda);
  
  return PM_SUCCESS;
}

int dense_solve_cholesky(Dense *A, DenseUPLO uplo, Dense* b) {
  lapack_int status;
  enum CBLAS_ORDER order;
  char ul;

  /* todo: check A and b have same order */
  
  order = (A->order == DENSE_COLUMN_MAJOR)? CblasColMajor:CblasRowMajor;  
  ul = (uplo==DENSE_LO)?'L':'U';
 
  status = LAPACKE_dpotrs(order,ul,A->rows,b->columns,A->data,A->tda,b->data,b->tda);
  
  return PM_SUCCESS;
}


/* Matrix Functions  */
int dense_diag_max(Dense *A, double *max) {
  double *p,m;
  int i,n;
  
  if (A->rows != A->columns) {
    PM_ERROR_NOMSG(PM_EMATRIXSIZE);
  }
  
  n = A->rows;
  p = A->data;
  m = *p;
  p += 1+A->tda;
  for(i=1; i < n; i++) {
    if (*p > m) m = *p;
    p += 1 + A->tda;
    
  }
  *max = m;
  return PM_SUCCESS;
}

int dense_triangular_logdet(Dense *A, double *logdet) {
  double *p,l;
  int i,n;

  if (A->rows != A->columns) {
    PM_ERROR_NOMSG(PM_EMATRIXSIZE);
  }
  
  n = A->rows;
  p = A->data;
  l = 0;
  for(i=0; i < n; i++) {
    l += log(*p);
    p += 1 + A->tda;
    
  }
  *logdet = l;
  return PM_SUCCESS;
}
