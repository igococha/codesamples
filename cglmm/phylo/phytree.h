#ifndef PHYTREE_H
#define PHYTREE_H


#include "list_any.h"


#include <stdio.h>

typedef enum { NEWICK, NEXUS } PhyFormat;

typedef int NodeId;

typedef struct PhyNode PhyNode;

/* Used mainly as temporary structure */
/* It should be discarded after construction of flat structures */
struct PhyNode {
  NodeId id;
  int tip;
  PhyNode* ancestor;
  int num_descendants;
  PhyNode** descendants;
  float edge_length;
  char* taxa;
};

typedef struct PhyTree {
  /* Temporary Node based tree structures */
  PhyNode *root, **internal_nodes, **tips;
  /* PhyNode **nodes, **internal_nodes, **tips; */
  unsigned int num_nodes, num_internal, num_tips;
  /* Id-based traversal structures */
  NodeId *children_all, **children, *parent;
  int *first_child, children_size, *num_children;
  /* traversals */
  int binary, max_descendants;
  NodeId *preorder;
  int *num_descendants_all, *num_descendants_tips;
  int *pos_preorder, *pos_first_tip;
  double *distances, *edge_length;
  char **taxa, *marked; 
  float* dist_root;
} PhyTree;

typedef struct PhyTreeList {
  PhyTree** trees;
  int num_trees;
} PhyTreeList;

/* ****** PHYNODE interface ******** */
PhyNode* phynode_new(NodeId id, PhyNode* anc, int is_tip);
void phynode_free(PhyNode** pnode);
void phynode_print(PhyNode* node, PhyTree* tree, FILE* fp);
int phynode_set_taxa(PhyNode* node, const char* str);
void phynode_set_descendants(PhyNode* node, ListAny desc);

/* ****** PHYTREE interface ******* */
PhyTree* phytree_new(PhyNode* root, ListAny internal_list,  ListAny tip_list);
PhyTree* phytree_free(PhyTree** ptree);
void phytree_print(PhyTree* tree, FILE* fp);
void phytree_print_tips(PhyTree* tree, FILE* fp);
int phytree_compare_tips(PhyTree* tree, char** data, int size, int *perm, int *num_perm);

/* ****** PHYTREELIST interface ******* */
PhyTreeList* phytreelist_build(char* filename, PhyFormat format);
PhyTreeList* phytreelist_new(ListAny tree_list);

#endif
