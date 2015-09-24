#include "pm_error.h"
#include "phytree.h"
#include "phyparser.h"
#include "pm_mem.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* ****** PHYNODE interface ******** */
PhyNode* phynode_new(NodeId id, PhyNode* anc, int is_tip) {
  PhyNode* node;
  node = (PhyNode*) PM_MEM_ALLOC(sizeof(PhyNode));
  if (node==NULL) return NULL;
  node->id = id;
  node->ancestor = anc;
  node->tip = is_tip;
  node->num_descendants = 0;
  node->descendants = NULL;
  node->taxa = NULL;
  return node;
}

void phynode_free(PhyNode **pnode) {
  PhyNode *node = *pnode;
  PM_MEM_FREE(node->descendants);
  /* free(node->taxa); hold this until decision is made */
  PM_MEM_FREE(*pnode);
}

void phynode_print(PhyNode* node, PhyTree* tree, FILE* fp) {
  int i;
  if (node->tip) {
    if (node->taxa==NULL) {
      fprintf(fp,"NULL");
    } else {
      fprintf(fp,"%s",node->taxa);
    }
  } else {
    fprintf(fp,"(");
    for(i=0; i < node->num_descendants-1; i++) {
      phynode_print(node->descendants[i],tree,fp);
      fprintf(fp,",");
    }
    phynode_print(node->descendants[i],tree,fp);
    fprintf(fp,")");
  }
  fprintf(fp,":%.2f",node->edge_length);
  //fprintf(fp,":%.2f[%d]",node->edge_length, node->id);
}

int phynode_set_taxa(PhyNode* node, const char* str) {
  int s = strlen(str)+1;
  node->taxa = (char*)PM_MEM_ALLOC(s*sizeof(char));
  if (node->taxa==NULL) {
    return pm_error_get_type(NULL);
  }
  strcpy(node->taxa,str);
  return PM_SUCCESS;
}

void phynode_set_descendants(PhyNode* node, ListAny desc) {
  unsigned short n = list_any_size(desc);;
  node->num_descendants = n;
  if (n < 1) {
    node->descendants = NULL;
  }
  node->descendants = (PhyNode**)list_any_clear_toarray(desc);
  return;
}

/* ****** PHYTREE interface ******* */
static void free_node_arrays(PhyTree* tree);
static int build_structures(PhyTree* tree, PhyNode* root);
static int build_visit_node(PhyNode* node, PhyTree* tree, int *pre_counter,  int* tip_counter);
static void phytree_print_node(NodeId node_id, PhyTree* tree, FILE*fp);

PhyTree* phytree_new(PhyNode* root, ListAny internal_list, ListAny tip_list) {
  PhyTree* tree;
  int i,node_id, children_size, children_idx;
  NodeId *pchildren;
  PhyNode *node,  **current, child;
  PhyNode **internal_nodes, **tips;
  
  tree = (PhyTree*)PM_MEM_ALLOC(sizeof(PhyTree));
  if (tree==NULL) return NULL;
  tree->root = root;
  tree->num_internal = list_any_size(internal_list);
  tree->num_tips = list_any_size(tip_list);

  /* allocate arrays of nodes */
  tree->internal_nodes = (PhyNode**)list_any_clear_toarray(internal_list);
  if (tree->internal_nodes==NULL) {
    PM_MEM_FREE(tree); return NULL;
  }
  internal_nodes = tree->internal_nodes;
  tree->tips = (PhyNode**)list_any_clear_toarray(tip_list);
  if (tree->tips==NULL) {
    free_node_arrays(tree); PM_MEM_FREE(tree); return NULL;
  }
  tips = tree->tips;
  
  tree->num_nodes = tree->num_internal + tree->num_tips;
  /* nodes = (PhyNode**)malloc(tree->num_nodes*sizeof(PhyNode*)); */
  current = internal_nodes;
  tree->binary = 1;
  /* Assign node ids */
  children_size = 0;
  for(node_id=0; node_id < tree->num_internal; node_id++, current++) {
    node = *current;
    node->id = node_id;
    /* nodes[node_id] = node; */
    i = node->num_descendants;
    children_size += i;
    if (i != 2) {
      tree->binary = 0;
    }
    if (i > tree->max_descendants) {
      tree->max_descendants = i;
    }
    
  }
  current = tips;
  for(; node_id < tree->num_nodes; node_id++, current++) {
    node = *current;
    node->id = node_id;
    /* nodes[node_id] = node; */
  }
  
  /* *** Allocate and initialize id-based traversal structures *** */
  tree->children_size = children_size;
  tree->children_all = NULL;
  tree->children = NULL;
  tree->num_children = NULL;
  tree->first_child = NULL;
  tree->parent = NULL;
  tree->edge_length = NULL;
  tree->taxa = NULL;
  tree->marked = NULL;
  tree->distances = NULL;
  
  tree->children_all = (NodeId*)PM_MEM_ALLOC(children_size * sizeof(NodeId));
  if (tree->children_all==NULL) return phytree_free(&tree);
  tree->children = (NodeId**)PM_MEM_ALLOC(tree->num_nodes * sizeof(NodeId**));
  if (tree->children==NULL) return phytree_free(&tree);
  tree->num_children = (int*)PM_MEM_ALLOC(tree->num_nodes * sizeof(int));
  if (tree->num_children==NULL) return phytree_free(&tree);
  tree->first_child = (int*)PM_MEM_ALLOC(tree->num_nodes * sizeof(int));
  if (tree->first_child==NULL) return phytree_free(&tree);
  tree->parent =  (NodeId*)PM_MEM_ALLOC(tree->num_nodes * sizeof(NodeId));
  if (tree->parent==NULL) return phytree_free(&tree);
  tree->edge_length = (double*)PM_MEM_ALLOC(tree->num_nodes * sizeof(double));
  if (tree->edge_length==NULL) return phytree_free(&tree);
  tree->taxa = (char**)PM_MEM_ALLOC(tree->num_nodes * sizeof(char*));
  if (tree->taxa==NULL) return phytree_free(&tree);
  tree->marked = (char*)PM_MEM_ALLOC(tree->num_nodes * sizeof(char));
  if (tree->marked==NULL) return phytree_free(&tree);
  tree->distances = (double*)PM_MEM_ALLOC(tree->num_nodes * sizeof(double));
  if (tree->distances==NULL) return phytree_free(&tree);

  
  children_idx = 0;
  pchildren = tree->children_all;
  current = internal_nodes;
  for(node_id=0; node_id < tree->num_internal; node_id++, current++) {
    node = *current;
    tree->parent[node_id] = (node->ancestor==NULL)?-1:node->ancestor->id;
    tree->edge_length[node_id] = node->edge_length;
    tree->taxa[node_id] = NULL; /* May not be the case */
    tree->first_child[node_id] = children_idx;
    tree->children[node_id] = pchildren;
    tree->num_children[node_id] = node->num_descendants;
    for(i=0; i < node->num_descendants; i++,pchildren++) {
      tree->children_all[children_idx++] = node->descendants[i]->id; 
    }
  }
  current = tips;
  for(; node_id < tree->num_nodes; node_id++, current++) {
    node = *current;
    tree->parent[node_id] = node->ancestor->id;
    tree->edge_length[node_id] = node->edge_length;
    tree->num_children[node_id] = 0;
    tree->taxa[node_id] = node->taxa; /* allocate? */
  }

  /* Allocate and initialize auxiliary structures */
  if (build_structures(tree,root)!=PM_SUCCESS) {
    return phytree_free(&tree);
  }

  /* Use flat structures to compute distances */
  tree->distances[0] = tree->edge_length[0];
  for(i=1; i < tree->num_nodes; i++) {
    node_id = tree->preorder[i];
    tree->distances[node_id] = tree->edge_length[node_id]+tree->distances[tree->parent[node_id]];
  }

  /* Get rid of PhyNodes */
  free_node_arrays(tree);

  return tree;
}

static void free_node_arrays(PhyTree* tree) {
  int i;
  if (tree->internal_nodes != NULL) {
    for(i=0; i < tree->num_internal; i++) {
      phynode_free(&(tree->internal_nodes[i]));
    }
    PM_MEM_FREE(tree->internal_nodes);
  }
  if (tree->tips != NULL) {
    for(i=0; i < tree->num_tips; i++) {
      phynode_free(&(tree->tips[i]));
    }
    PM_MEM_FREE(tree->tips);
  }
  tree->root = NULL;

}


static int build_structures(PhyTree* tree, PhyNode* root) {
  int pre_counter, tip_counter;
  int i,err;
  
  tree->preorder = (NodeId*)PM_MEM_ALLOC(tree->num_nodes * sizeof(NodeId));
  if (tree->preorder==NULL) return pm_error_get_type(NULL);
  tree->pos_preorder = (int*)PM_MEM_ALLOC(tree->num_nodes * sizeof(int));
  if (tree->pos_preorder==NULL) return pm_error_get_type(NULL);
  tree->pos_first_tip = (int*)PM_MEM_ALLOC(tree->num_nodes * sizeof(int));
  if (tree->pos_first_tip==NULL) return pm_error_get_type(NULL);
  /* preorder_tips: m,m+1,m+2,...,n-1 where m=num_internal, n=num_nodes  */
  tree->num_descendants_all = (int*)PM_MEM_ALLOC(tree->num_nodes * sizeof(int));
  if (tree->num_descendants_all==NULL) return pm_error_get_type(NULL);
  tree->num_descendants_tips=(int*)PM_MEM_ALLOC(tree->num_nodes * sizeof(int));
  if (tree->num_descendants_tips==NULL) return pm_error_get_type(NULL);
 
  pre_counter = 0;
  tip_counter = 0;

  build_visit_node(root,tree,&pre_counter,&tip_counter);
  
  return PM_SUCCESS;
}

static int build_visit_node(PhyNode* node, PhyTree* tree, int *pre_counter, int *tip_counter) {
  int i, n, node_id, pre_counter_start, num_tips, nt;

  pre_counter_start = *pre_counter;
  node_id = node->id;
  tree->preorder[*pre_counter] = node_id;
  tree->pos_preorder[node_id] = *pre_counter;
  tree->pos_first_tip[node_id] = *tip_counter;

  (*pre_counter)++;

  if (node->tip) {
    (*tip_counter)++;
    tree->num_descendants_all[node_id] = 1;
    tree->num_descendants_tips[node_id] = 1;
    num_tips=1;
    /* printf("Tip %d \n",node->id); */ 
  } else {
    /* printf("Node %d \n",node->id); */ 
    num_tips = 0;
    for(i=0; i < node->num_descendants; i++) {
      nt = build_visit_node(node->descendants[i],tree,pre_counter,tip_counter);
      num_tips += nt;
    }
    tree->num_descendants_all[node_id] = *pre_counter - pre_counter_start;
    tree->num_descendants_tips[node_id] = num_tips;
  }
  return num_tips;
}


PhyTree* phytree_free(PhyTree** ptree) {
  PhyTree* tree = *ptree;
  PhyNode** nodes;
  int i;
  /* Auxiliary arrays */
  PM_MEM_FREE(tree->num_descendants_tips);
  PM_MEM_FREE(tree->num_descendants_all);
  PM_MEM_FREE(tree->distances);
  PM_MEM_FREE(tree->pos_first_tip);
  PM_MEM_FREE(tree->pos_preorder);
  PM_MEM_FREE(tree->preorder);

  PM_MEM_FREE(tree->marked);
  for(i=0; i < tree->num_nodes; i++) { PM_MEM_FREE(tree->taxa[i]); }
  PM_MEM_FREE(tree->taxa);
  PM_MEM_FREE(tree->edge_length);
  PM_MEM_FREE(tree->parent);
  PM_MEM_FREE(tree->first_child);
  PM_MEM_FREE(tree->num_children);
  PM_MEM_FREE(tree->children);

  free_node_arrays(tree);

  PM_MEM_FREE(*ptree);

  return NULL;
}

void phytree_print(PhyTree* tree, FILE* fp) {
  phytree_print_node(0,tree,fp);
  fprintf(fp,"\n");
}

void phytree_print_tips(PhyTree* tree, FILE* fp) {
  int node_id;
  for(node_id=tree->num_internal; node_id < tree->num_nodes; node_id++) {
    printf("%s ",tree->taxa[node_id]);
  }
  printf("\n");
}

static void phytree_print_node(NodeId node_id, PhyTree* tree, FILE*fp) {
  int i;
  NodeId *children;
  char** taxa = tree->taxa;
  if (node_id >= tree->num_internal) { /* tip */
    if (taxa[node_id]==NULL) {
      fprintf(fp,"NULL");
    } else {
      fprintf(fp,"%s",taxa[node_id]);
    }
  } else {
    fprintf(fp,"(");
    children = tree->children[node_id];
    for(i=0; i < tree->num_children[node_id]-1; i++) {
      phytree_print_node(children[i],tree,fp);
      fprintf(fp,",");
    }
    phytree_print_node(children[i],tree,fp);
    fprintf(fp,")");
  }
  fprintf(fp,":%.2f",tree->edge_length[node_id]);
  //fprintf(fp,":%.2f[%d]",node->edge_length, node->id);
}

int phytree_compare_tips(PhyTree* tree, char** data, int size, int *perm, int *num_perm) {
  int node_id, first_tip,num_nodes,i,j;
  char *marked;
  
  num_nodes = tree->num_nodes;
  if (size != tree->num_tips) return -1; /* size mistmatch */
  first_tip = tree->num_internal;
  if (perm==NULL) { /* simple version , one to one map */
    for(node_id = first_tip,i=0; i < size; node_id++,i++) {
      if (strcmp(tree->taxa[node_id],data[i])!=0)
	break;
    }
  } else { /* perm is not NULL */
    *num_perm = 0;
    for(node_id=first_tip; node_id < num_nodes; node_id++) {
      tree->marked[node_id]=0; /* initialize array of marked nodes */      
    }
    /* for each taxa data[i] in input array */
    for(i=0; i < size; i++) {
      /* lookup in tree */
      for(node_id = first_tip,j=0; node_id < num_nodes; node_id++,j++) {
	if (!tree->marked[node_id]) {
	  if (strcmp(tree->taxa[node_id],data[i])==0) { /* Found */
	    perm[i] = j;
	    if (i != j) (*num_perm)++;
	    tree->marked[node_id] = 1; /* don't check again */
	    break;
	  }
	}
      } /* inner for */
      if (node_id == num_nodes)
	break;
    } 
  }
  return i;
}


/* ****** PHYTREELIST interface ******* */
PhyTreeList* phytreelist_build(char* filename, PhyFormat format) {
  PhyParser phyparser;
  PhyTreeList* tree_list = NULL;
  FILE* fp;
  int i;

  pm_error_clear(NULL);

  fp = fopen(filename,"r");
  if (fp == NULL) {
    PM_ERROR_VAL(PM_EFILE,NULL,"Couldn't open tree file: %s",filename);
  }
  phyparser = phyparser_new(fp,format);
  if (phyparser != NULL) {
    tree_list = phyparser_parse(phyparser,pm_error_default());
  } else {
    fclose(fp);
    return NULL;
  }
  fclose(fp);
  phyparser_free(&phyparser);


  return tree_list;
}

PhyTreeList* phytreelist_new(ListAny tree_list) {
  int num_trees = list_any_size(tree_list);
  PhyTreeList* result;
  result = (PhyTreeList*)PM_MEM_ALLOC(sizeof(PhyTreeList));
  if (result==NULL) return NULL;
  result->num_trees = num_trees;
  result->trees = (PhyTree**)list_any_clear_toarray(tree_list);
  if (result->trees==NULL) {
    PM_MEM_FREE(result); return NULL;
  }
  return result;
}
