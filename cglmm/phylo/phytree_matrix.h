#ifndef PHYTREE_MATRIX_H
#define PHYTREE_MATRIX_H

#include "phytree.h"
#include "dense.h"
#include "sparse.h"

Dense* phytree_compute_matrix(PhyTree *tree, int all_nodes);
Dense* phytree_compute_matrix_tips(PhyTree *tree);
Dense* phytree_compute_matrix_all(PhyTree *tree, int noroot);

Dense* phytree_compute_inverse_dense(PhyTree *tree);
Sparse* phytree_compute_inverse_sparse(PhyTree *tree);

/* Dense* phytree_compute_inverse(PhyTree *tree); */


#endif
