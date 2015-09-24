#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>


#include "nlopt.h"

#include "pm_error.h"
#include "pm_mem.h"
#include "phyanimal.h"
#include "phytree_matrix.h"
#include "list_all.h"

struct PhyAnimal {
  ListStr yvar;  /* response variables */
  ListStr fvar; /* fixed variables */
  int n,p,num_tips, all_nodes;
  char* taxa_column;
  Dense* X;
  DenseVector* y;
  Dense* V, *L;
  /* temporaries */
  Dense *A, *A1;
  DenseVector *B, *B1, *E;
  double sigma, logdet, loglh, lambda;
};

static void phyanimal_fit_gls(PhyAnimal model, const double *lambda);
static void phyanimal_fit_lambda(PhyAnimal model);
static double loglh_worker(unsigned n, const double *x, double *grad, void *fdata);

static void str_print(char **ele,void* cl) {
  fprintf((FILE*)cl,"%s ",*ele);
}

static void free_arrays(PhyAnimal m) {
  PM_MEM_FREE(m->taxa_column);
  dense_free(&m->X);
  dense_vector_free(&m->y);
  dense_free(&m->V);
  return;
}

static int allocate_temporaries(PhyAnimal model) {
  int n,p;
  n= model->n;
  p = model->p;
  model->L = dense_new(n,n);
  model->A = dense_new(p+1,p+1);
  model->A1 = dense_new(n,p+1);
  /* vectors */
  model->B = dense_vector_new(p+1);
  model->B1 = dense_vector_new(n);
  model->E = dense_vector_new(n); /* could reuse B1 or use a vector */
  return PM_SUCCESS;
}

static void free_temporaries(PhyAnimal model) {
  dense_free(&model->L);
  dense_free(&model->A);
  dense_free(&model->A1);
  dense_vector_free(&model->B);
  dense_vector_free(&model->B1);
  dense_vector_free(&model->E);
}


/* Returns new Animal model, NULL if there was an error */
PhyAnimal phyanimal_new(char** y, int num_y, char** x, int num_x) {
  PhyAnimal model;
  int i;
  char *name;
  ListStr l;

  pm_error_clear(NULL);
  model = (PhyAnimal) PM_MEM_ALLOC(sizeof(struct PhyAnimal));
  if (model==NULL) return NULL;

  model->yvar = NULL;
  model->fvar = NULL;
  model->taxa_column = NULL;
  model->X = NULL;
  model->y = NULL;
  model->V = NULL;
  /* initialization */
  model->yvar = l = list_str_new();
  for(i=0; i<num_y; i++) {
    name = (char*) PM_MEM_ALLOC(strlen(y[i])+1);
    if (name==NULL) {
      phyanimal_free(&model);
      return NULL;
    }
    strcpy(name,y[i]);
    list_str_push_back(l,name);
  }
  model->fvar = l = list_str_new();
  for(i=0; i<num_x; i++) {
    name = (char*) PM_MEM_ALLOC(strlen(x[i])+1);
    if (name==NULL) {
      phyanimal_free(&model);
      return NULL;
    }
    strcpy(name,x[i]);
    list_str_push_back(l,name);
  }
  model->p = num_x;
  model->n = 0;
  return model;
}

Dense* phyanimal_X(PhyAnimal model) {
  return model->X;
}
DenseVector* phyanimal_y(PhyAnimal model) {
  return model->y;
}

Dense* phyanimal_V(PhyAnimal model) {
  return model->V;
}


void phyanimal_print(PhyAnimal model,FILE* fp) {
  fprintf(fp,"Basic Animal Model:\n");
  list_str_map(model->yvar,str_print ,fp);
  fprintf(fp," ~ ");
  list_str_map(model->fvar,str_print ,fp);
  fprintf(fp,"\n");
}


void phyanimal_free(PhyAnimal *model) {
  ListStr l;
  char *name;
  PhyAnimal m = *model;
  if (*model==NULL) return;
  free_temporaries(*model);
  free_arrays(m);
  /* Variable names */
  if (m->fvar != NULL) {
    l = m->fvar;
    while (!list_str_empty(l)) {
      list_str_pop_front(l,&name);
      PM_MEM_FREE(name);
    }
    list_str_free(&l);
  }
  if (m->yvar != NULL) {
    l = m->yvar;
    while (!list_str_empty(l)) {
      list_str_pop_front(l,&name);
      PM_MEM_FREE(name);
    }
    list_str_free(&l);
  }

  PM_MEM_FREE(*model);
  *model = NULL;
  return;
}


int phyanimal_load_data(PhyAnimal model,DataFrame *df, char *taxa_column, PhyTree *tree, int all_nodes) {
  int status, num_rows, *perm, num_perm;
  DataType datatype;
  void *column;
  char **taxa_names;
  char** names, *name;
  double max;
  
  pm_error_clear(NULL);
  phyanimal_print(model,stdout);
  model->all_nodes = all_nodes;
  /* check that taxacolumn is in dataframe */
  datatype = DF_NOTYPE;
  num_rows = df_num_rows(df);
  column = df_get_column(df,taxa_column,&datatype); /* No new memory */
  if (column != NULL) {
    if (datatype != DF_STRING) {
      PM_ERROR(PM_EDF,
	       "error while loading dataframe\nExpecting column %s to be of type string",
	       taxa_column);
    }
    taxa_names = (char**)column;
  } else {
    PM_ERROR(PM_EDF,"Taxa column %s not in dataframe",taxa_column);
  }
  /* check that taxa_names match tips */
  perm = (int*) PM_MEM_ALLOC(sizeof(int) * num_rows);
  if (perm==NULL) return pm_error_get_type(NULL);
  status = phytree_compare_tips(tree,taxa_names,num_rows,perm,&num_perm);
  if (status != num_rows) {
    if (status >= 0) {
      PM_ERROR(PM_EMODEL,"Column %s does not match tip names.\nFailed at %s",taxa_column, taxa_names[status]);
    } else {
      PM_ERROR(PM_EMODEL,"Column %s does not match tip names. ",taxa_column);
    }
  }
  if (num_perm > 0) {
    printf("Performing %d permutations to dataframe\n",num_perm);
    df_permute(df,perm);
    PM_MEM_FREE(perm);
    taxa_names = (char**)df_get_column(df,taxa_column,&datatype);
    status = phytree_compare_tips(tree,taxa_names,num_rows,NULL,NULL);
    if (status!=num_rows) {
      printf("Horror\n"); exit(0);
    }
  }
  /* Generate V */
  model->V = phytree_compute_matrix(tree,all_nodes);
  if (model->V == NULL) {
    return pm_error_get_type(NULL);
  }

  /* dense_print(model->V,stdout); */
  /* if (model->V != NULL) dense_print(model->V,stdout); */
  /* Normalize V matrix */
  dense_diag_max(model->V,&max);
  dense_mbyc(model->V,1.0/max);
  model->num_tips = num_rows;
  model->n = model->V->rows;      

  /* taxa_column is correct */
  model->taxa_column = (char*) PM_MEM_ALLOC(sizeof(taxa_column)+1);
  if (model->taxa_column != NULL) {
    strcpy(model->taxa_column,taxa_column);
    names = list_str_toarray(model->fvar);    
  }
  if (pm_error_get_type(NULL) != PM_SUCCESS) {
    free_arrays(model);
    return pm_error_get_type(NULL);
  }

  /* generate X */
  printf("Prepending %d rows to X and y\n",model->n-model->num_tips);
  model->X = df_gen_matrix(df,names,list_str_size(model->fvar),1,model->n - model->num_tips);
  PM_MEM_FREE(names);
  if (model->X != NULL) {
    /* printf("X\n"); dense_print(model->X,stdout); */
    /* generate y */
    list_str_front(model->yvar,&name);
    model->y = df_gen_vector(df,name,model->n - model->num_tips);
    /*dense_vector_set_range(model->y,DBL_MAX,0,model->n-model->num_tips-1);*/
    
    /* dense_vector_print(model->y,stdout); */
    /* model->y = df_gen_matrix(df,names,1,0); */
  }

  if (pm_error_get_type(NULL) != PM_SUCCESS) {
    free_arrays(model);
    return pm_error_get_type(NULL);
  }
  status = allocate_temporaries(model);

  return pm_error_get_type(NULL);
}

void phyanimal_fit(PhyAnimal model, int comp_lambda) {
  
  if (comp_lambda) {
    phyanimal_fit_lambda(model);

  } else {
    phyanimal_fit_gls(model,NULL);
    model->loglh = phyanimal_loglh(model);
  }
  printf("Loglh = %f\n",model->loglh);
  printf("Fit finished\n");

}


static void phyanimal_fit_gls(PhyAnimal model, const double *lambda) {
  int status;

  dense_copy(model->L, model->V);
  /* update matrix if lambda is present */
  if (lambda != NULL) {
    printf("--- lambda %f\n",*lambda);
    dense_mbyc_nondiag(model->L,*lambda);
  }
  dense_cholesky(model->L,DENSE_LO); /* L = Chol(V) */
  dense_triangular_logdet(model->L,&model->logdet);
  model->logdet *= 2;

  /*  L*A1 = X   */
  dense_copy(model->A1, model->X);
  dense_solve_triangular(model->L,DENSE_LO,DENSE_TRANS_NO,model->A1);
  /* A = A1^t * A1 */
  dense_mtbym(model->A,model->A1);
  /* L*B1 = y */
  dense_vector_copy(model->B1,model->y);

  dense_vector_solve_triangular(model->L,DENSE_LO,DENSE_TRANS_NO,model->B1);
  
  /* B = A1^t * B1 */
  /*dense_mbym(model->B,1.0,DENSE_TRANS,model->A1,DENSE_TRANS_NO,model->B1,0.0);*/
  dense_vector_mbyv_inc(model->B,1.0,DENSE_TRANS,model->A1,model->B1,0.0);
  
  /* A*beta = B */
  dense_vector_solve_positive(model->A,DENSE_LO,model->B);

  printf("Beta = ");
  dense_vector_print(model->B,stdout);
  
  /* Compute sigma */
  dense_vector_copy(model->E, model->y); /* vectors represented as matrices  */
  dense_vector_mbyv_inc(model->E,-1.0,DENSE_TRANS_NO,model->X,model->B,1.0);
  /* printf("\nError\n");
     dense_vector_print(model->E,stdout); */
  /* sigma: S = E' * V^-1 * E = (Linv E)' (Liinv E) = S1' * S1  */
  dense_vector_solve_triangular(model->L,DENSE_LO,DENSE_TRANS_NO,model->E); /* S1 stored in E */
  /* dense_mtbym(model->S,model->E); */

  model->sigma = dense_vector_vtbyv(model->E); /* v^2 */
  model->sigma /= (model->n - 1);
  printf("\nsigma %f\n",model->sigma);
  return;
}

double phyanimal_loglh(PhyAnimal model) {
  double lh, logdet;
   
  lh = (model->n)*log(2*M_PI*model->sigma) + model->logdet;
  /* (E'*invV*E)/sigma but sigma=(E'*invV*E)/(n-1) */
  lh += (model->n - 1);
  lh /= -2;

  return lh;
}

static void phyanimal_fit_lambda(PhyAnimal model) {
  nlopt_opt opt;
  double lb, ub, lambda, loglh;

  
  opt = nlopt_create(NLOPT_LN_NELDERMEAD, 1); /* one dimension */
  lb = 0.0;
  ub = 1.0;
  nlopt_set_lower_bounds(opt,&lb);
  nlopt_set_upper_bounds(opt,&ub);
  /* No inequality constraints */
  nlopt_set_max_objective(opt, loglh_worker, model);
  nlopt_set_xtol_rel(opt, 1e-4);
  lambda = 0.1;
  if (nlopt_optimize(opt, &lambda, &loglh)<0) {
    printf("NLOPT failed\n");
  } else {
    printf("Maximum lh %f at lambda %f\n",loglh,lambda);
    printf("Sigma = %f\n",model->sigma);
    printf("Sigma_b = %f\nSigma_e = %f\n",model->sigma*lambda,model->sigma*(1-lambda));
  }
  model->loglh = loglh;
  model->lambda = lambda;
  nlopt_destroy(opt);
  return;
}

static double loglh_worker(unsigned n, const double *lambda, double *grad,
			   void *fdata) {
  PhyAnimal model;
  double loglh;
  /* n = 1. One dimensional optimization */
  /* grad should be NULL */
  model = (PhyAnimal)fdata;
  phyanimal_fit_gls(model,lambda);
  loglh = phyanimal_loglh(model);
  printf("loglh %f\n", loglh);
  return loglh;
}
