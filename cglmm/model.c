#include <stdlib.h>
#include <stdio.h>

#include "phylexer.h"
#include "phyparser.h"
#include "phytree.h"
#include "pm_error.h"
#include "pm_mem.h"

#include "list_any.h"
#include "list_all.h"

#include "phytree_matrix.h"
#include "dense.h"
#include "sparse.h"

#include "cblas.h"
#include "lapacke.h"

#include "df.h"

#include "math.h"
#include "nlopt.h"

#include "phyanimal.h"
#include "pm_rng.h"



int main(int argc, char* argv[]) {
  int i,n, status;
  char* datafile, *treefile, *taxacolumn;
  PhyAnimal model;
  DataFrame *df;
  PhyTreeList *treelist;
  double lh;
  pm_rng *rng;
  DenseVector *v,*v2;
  mc2_animal mc2;

  /*  test random number generator */
  /*
  rng = pm_rng_new(pm_rng_default);
  v = dense_vector_new(30);
  dense_vector_print(v,stdout);

  dense_vector_ran_gaussian(v,rng,5,1);
  printf("done\n");
  dense_vector_print(v,stdout);
  exit(0);
  v2 = dense_vector_new(15);
  dense_vector_copy_idx(v2,2,v,0, 100);
  printf("\n");
  dense_vector_print(v2,stdout);
  */

  /*Dense *A, *B;
  A = dense_new(4,3);
  B = dense_new(3,4);
  dense_set_diag(A,2);
  dense_print(A,stdout);
  
  if (dense_copy_trans(B,A)!=0) {
     printf("wrong arguments \n");
    exit(0);
   }
  dense_print(B,stdout);
  */
  
  if (argc < 5) {
    printf("Insufficient number of arguments\n");
    printf("Syntax: <DataFile> <TreeFile> y x1 x2 x3 taxaColumn\n");
    exit(EXIT_FAILURE);
  }
  datafile = argv[1];
  treefile = argv[2];
  taxacolumn = argv[argc-1];
  printf("Data %s Tree %s\n",datafile,treefile);

  if(argc == 5) {
    model = phyanimal_new(&argv[3],1,NULL,0);
  } else {
    model = phyanimal_new(&argv[3],1,&argv[4],argc-5);
  }
  if (model == NULL) {
    pm_error_print(NULL,stdout);
    exit(EXIT_FAILURE);
  }

  phyanimal_print(model,stdout);
  
  /* Read DataFrame */
  df = df_load_new(datafile);
  if (df == NULL) {
    pm_error_print(NULL,stdout);
    exit(EXIT_FAILURE);
  }

  /* read tree */
  treelist = phytreelist_build(treefile,NEWICK);
  if (treelist == NULL) {
    printf("Couldn't build tree\n");
    pm_error_print(NULL,stdout);
    exit(EXIT_FAILURE);
  }
  
  if (treelist->num_trees > 0) {
    status = phyanimal_load_data(model,df,taxacolumn,treelist->trees[0],0);
    if (status != PM_SUCCESS) {
      printf("Couldn't load data\n");
      pm_error_print(NULL,stdout);
      exit(EXIT_FAILURE);
    }
  }
  
  printf("Fitting\n");

  phyanimal_fit(model,0);
  phyanimal_print(model,stdout);

  printf("MCMC\n");
  
  mc2 = mc2_animal_new(model);
  status = mc2_animal_run(mc2,30,5,10);
  if (status != 0) {
    printf("mcmc simulation failed\n");
  }

  mc2_animal_free(&mc2);
  df_free(&df);
  phyanimal_free(&model);
  exit(EXIT_SUCCESS);
}
