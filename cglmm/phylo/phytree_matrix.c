#include "phytree_matrix.h"
#include "pm_mem.h"
#include <stdlib.h>

/* Local declarations */



/* End local declarations */

Dense* phytree_compute_matrix(PhyTree *tree, int all_nodes) {
  Dense *matrix;
  if (all_nodes) {
    matrix = phytree_compute_matrix_all(tree,1); 
  } else {
    matrix = phytree_compute_matrix_tips(tree);
  }
 
  return matrix;
}

/* Local function definitions */

Dense* phytree_compute_matrix_tips(PhyTree* tree) {
  Dense* matrix;
  double *m, d;
  int ntips, nnodes, ninternal, ndesc, node_id,row,col;
  int i,j,node_id1,n1,idx1,first_idx1,node_id2,n2,idx2,first_idx2;
  NodeId *children;

  ntips = tree->num_tips;
  ninternal = tree->num_internal;
  nnodes = tree->num_nodes;
  
  matrix = dense_new(ntips,ntips);
  if (matrix == NULL) {
    /* phyerror stuff here */
    return NULL;
  }
  dense_set_all(matrix,0);
  m = DENSE_GET_DATA(matrix);
  /* Set diagonal elements */
  row = ntips-1;
  for(node_id = nnodes-1; node_id >= ninternal; node_id--,row--) {
    m[row*ntips + row] = tree->distances[node_id];
  }
  /* traverse internal nodes */
  for(node_id = ninternal-1; node_id >=0; node_id--) {
    d = tree->distances[node_id];
    ndesc = tree->num_children[node_id];
    children = tree->children[node_id];
    /* cartesian product between descendant groups */
    for(i=0; i < ndesc; i++) {
      node_id1 = children[i]; /* changed node->descendants[i]->id; */
      n1 = tree->num_descendants_tips[node_id1]; 
      first_idx1 = tree->pos_first_tip[node_id1];
      /* match against siblings */
      if (i < ndesc-1) {
	for(j=i+1; j < ndesc; j++) {
	  node_id2 = children[j];  /* node->descendants[j]->id; */
	  n2 =  tree->num_descendants_tips[node_id2];
	  first_idx2 = tree->pos_first_tip[node_id2];
	  for(idx1=first_idx1; idx1 < (first_idx1+n1); idx1++) {
	    /* row = preorder[idx1] - num_internal; */
	    row = idx1;
	    for(idx2=first_idx2; idx2<(first_idx2+n2); idx2++) {
	      /* col = preorder[idx2] - num_internal; */
	      col = idx2;
	      m[row*ntips + col] = d;
	      m[col*ntips + row] = d;
	    }
	  }
	}
      }
    }
  }
  return matrix;
}

Dense* phytree_compute_matrix_all(PhyTree* tree, int noroot) {
  Dense* matrix;
  double *m,*pm, d;
  int ntips, nnodes, ninternal, ndesc, node_id,row,col;
  int i,j,node_id1,n1,idx1,first_idx1,node_id2,n2,idx2,first_idx2;
  /* PhyNode* node; */
  NodeId *children;
  int start,dim;

  ntips = tree->num_tips;
  ninternal = tree->num_internal;
  start = (noroot)?1:0;
  nnodes = tree->num_nodes;
  dim = nnodes-start;
  
  matrix = dense_new(dim,dim);
  m = DENSE_GET_DATA(matrix);
  /* set to zero */
  dense_set_all(matrix,0);

  /* Set diagonal elements */
  row=0;
  for(node_id = start; node_id < nnodes; node_id++,row++) {
    m[row*dim + row] = tree->distances[node_id];
  }
  /* Traverse internal nodes */
  for(node_id = 0; node_id < ninternal; node_id++) {
    d = tree->distances[node_id];
    ndesc = tree->num_children[node_id];
    children = tree->children[node_id];
    /* cartesian product between descendant groups */
    for(i=0; i < ndesc; i++) {
      node_id1 = children[i]; /* node->descendants[i]->id; */
      n1 = tree->num_descendants_all[node_id1];
      first_idx1 = tree->pos_preorder[node_id1];
      /* shared path parent - all descendants */
      /* skip for first row/col if root is not included */
      if (node_id >= start) {
	row=node_id-start;
	for(idx1=first_idx1; idx1 < (first_idx1+n1); idx1++) {
	  col = tree->preorder[idx1]-start; /* correct index */
	  m[row*dim + col] = d;
	  m[col*dim + row] = d;
	}
      }
      /* match againts siblings */
      if (i < ndesc-1) {
	for(j=i+1; j < ndesc; j++) {
	  node_id2 = children[j]; /* node->descendants[j]->id; */
	  n2 = tree->num_descendants_all[node_id2];
	  first_idx2 = tree->pos_preorder[node_id2];
	  for(idx1=first_idx1; idx1 < (first_idx1+n1); idx1++) {
	    row = tree->preorder[idx1]-start;
	    for(idx2=first_idx2; idx2<(first_idx2+n2); idx2++) {
	      col = tree->preorder[idx2]-start;
	      m[row*dim + col] = d;
	      m[col*dim + row] = d;
	    }
	  }	
	}
      }
    }
  }
  return matrix;
}


Dense* phytree_compute_inverse_dense(PhyTree *tree) {
  Dense *matrix;
  double *m,v;
  int nnodes;
  int root,node_id, node_idx,parent_id,parent_idx,i,dim;

  nnodes = tree->num_nodes;
  dim = nnodes-1;
  matrix = dense_new(dim,dim);
  dense_set_all(matrix,0);
  m = DENSE_GET_DATA(matrix);

  root = tree->preorder[0]; /* should be zero */
  for(i=1; i < nnodes; i++) {
    node_id = tree->preorder[i];
    parent_id = tree->parent[node_id];
    /* adjust index: we are skipping the root */
    node_idx = node_id-1;
    
    v = 1.0/tree->edge_length[node_id];
    m[node_idx*dim + node_idx] = v;
    if (parent_id != root) {
      parent_idx = parent_id-1;
      m[parent_idx*dim + parent_idx] += v;
      m[parent_idx*dim + node_idx] = -v;
      m[node_idx*dim + parent_idx] = -v;
    }
  }
  
  return matrix;

}

Sparse* phytree_compute_inverse_sparse(PhyTree *tree) {
  Sparse *M,*S;
  double v, *diag;
  int nnodes,dim,nzmax;
  int root,node_id, node_idx,parent_id,parent_idx,i;
  nnodes = tree->num_nodes;
  dim = nnodes-1;
  nzmax = dim + 2*tree->num_internal*tree->max_descendants;
  
  diag = (double*)PM_MEM_ALLOC(dim*sizeof(double));
  if (diag==NULL) return NULL;
  M = sparse_spalloc(dim,dim,nzmax,1,1); /* triplet */

  root = tree->preorder[0]; /* should be zero */
  for(i=1; i < nnodes; i++) {
    node_id = tree->preorder[i];
    parent_id = tree->parent[node_id];
    /* adjust index: we are skipping the root */
    node_idx = node_id-1;    
    v = 1.0/tree->edge_length[node_id];
    /* m[node_idx*dim + node_idx] = v; */
    diag[node_idx] = v;
    if (parent_id != root) {
      parent_idx = parent_id-1;
      /* m[parent_idx*dim + parent_idx] += v; */
      diag[parent_idx] += v;
      /* m[parent_idx*dim + node_idx] = -v; */
      sparse_entry(M,parent_idx,node_idx,-v);
      /* m[node_idx*dim + parent_idx] = -v; */
      sparse_entry(M,node_idx,parent_idx,-v);
    }
  }
  for(i=0; i < dim; i++) {
    sparse_entry(M,i,i,diag[i]);
  }
  PM_MEM_FREE(diag);
  S = sparse_compress(M);
  sparse_free(M);

  return S;
}
