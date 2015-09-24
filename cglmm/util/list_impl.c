/* 
Author: Igor Siveroni
 */
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "pm_error.h"
#include "pm_mem.h"

#define EMPTY(l) (l->size == 0)


#define LNODE(X) CAT(ListNode,X)


/* local typedef */
typedef struct LNODE(T) LNODE(T);

struct LNODE(T) {
  LNODE(T)* next;
  T e;
};

struct LT {
  int size;
  LNODE(T) *first,*last;
};

LT METHOD(list,T,new)() {
  LT l;
  /* LT l = (LT)malloc(sizeof(struct LT)); */
  l = (LT) PM_MEM_ALLOC(sizeof(struct LT));
  if (l==NULL) {
    pm_error_append_message(NULL,"Failed to allocate List node");
    return NULL;
  }
  l->first = l->last = NULL;
  l->size = 0;
  return l;
}

int METHOD(list,T,size)(LT l) {
  return l->size;
}


int METHOD(list,T,empty)(LT l) {
  return EMPTY(l);
}


int METHOD(list,T,front)(LT l,T* x) {
  /* if debug - insert assert */
  if (EMPTY(l)) {
    PM_ERROR(PM_EARGUMENTS,"invalid operation on empty list");
  }
  *x =  l->first->e;
  return PM_SUCCESS;
}

int METHOD(list,T,back)(LT l,T* x) {
  /* if debug - insert assert */
  if (EMPTY(l)) {
    PM_ERROR(PM_EARGUMENTS,"invalid operation on empty list");
  }
  *x = l->last->e;
  return PM_SUCCESS;
}


int METHOD(list,T,pop_front)(LT l,T* x) {
  LNODE(T)* node;
  
  /* if debug - insert assert */
  if (EMPTY(l)) {
    PM_ERROR(PM_EARGUMENTS,"invalid operation on empty list");
  }
  node = l->first;
  if (node == l->last) {
    l->last = NULL;
  }
  l->first = node->next;
  l->size--;
  *x = node->e;
  free(node);
  return PM_SUCCESS;
}

int METHOD(list,T,pop_back)(LT l,T* x) {
  LNODE(T) *node, *p;
  
  /* if debug - insert assert */
  if (EMPTY(l)) {
    PM_ERROR(PM_EARGUMENTS,"invalid operation on empty list");
  }
  node = l->last;
  if (l->first == node) {
    l->first = NULL;
    l->last = NULL;
  } else {
    p = l->first;
    while (p->next != node) {
      p = p->next;
    }
    p->next = NULL;
    l->last = p;
  }
  l->size--;
  *x = node->e;
  free(node);
  return PM_SUCCESS;
}

int METHOD(list,T,element)(LT l, int index, T* x) {
  LNODE(T)* node;

  if ((index < 0)||(index >= l->size)) {
    PM_ERROR(PM_EOUTOFBOUNDS,"Invalid index value for list element operation");
  }
  node = l->first;
  while (index > 0) {
    node = node->next;
    index--;
  }
  *x = node->e;
  return PM_SUCCESS;
}

void METHOD(list,T,push_front)(LT l, T e) {
  LNODE(T)* node = (LNODE(T)*)malloc(sizeof(struct LNODE(T)));
  node->e = e;
  node->next = l->first;
  l->first = node;
  if (l->last == NULL) {
    l->last = node;
  }
  l->size++;
  return;
}

void METHOD(list,T,push_back)(LT l, T e) {
  LNODE(T)* node = (LNODE(T)*)malloc(sizeof(struct LNODE(T)));
  node->e = e;
  node->next = NULL;
  if (l->first == NULL) {
    l->first = node;
  } else {
    l->last->next = node;
  }
  l->last = node;
  l->size++;
  return;
}


void METHOD(list,T,clear)(LT l) {
  LNODE(T) *pnode, *pnext;
  if EMPTY(l) {
      return;
  }
  pnode = l->first;
  while (pnode != NULL) {
    pnext = pnode->next;
    free(pnode);
    pnode = pnext;
  }
  l->size = 0;
  l->first  = l->last = NULL;
}

void METHOD(list,T,free)(LT* l) {
  METHOD(list,T,clear)(*l);
  free(*l);
  *l = NULL;
  return;
}


void METHOD(list,T,map)(LT l, void apply(T *x, void *cl), void *cl) {
  LNODE(T) *pnode = l->first;
  for( ; pnode; pnode = pnode->next) {
    apply(&pnode->e,cl);
  }
}



T* METHOD(list,T,toarray)(LT l) {
  LNODE(T) *pnode = l->first;
  T* array;
  int i;
  if (EMPTY(l)) {
    PM_ERROR_VAL(PM_EARGUMENTS,
		 NULL,"clear_toarray:invalid operation on empty array");
  }
  array = (T*)PM_MEM_ALLOC(l->size * sizeof(T));
  if (array==NULL) return NULL;
  for(i=0 ; pnode; pnode = pnode->next, i++) {
    array[i] = pnode->e;
  }
  return array;
}

T* METHOD(list,T,clear_toarray)(LT l) {
  LNODE(T) *pnext, *pnode;
  T* array;
  int i;
  if (EMPTY(l)) {
    PM_ERROR_VAL(PM_EARGUMENTS,
		 NULL,"clear_toarray:invalid operation on empty array");
  }
  array = (T*)PM_MEM_ALLOC(l->size * sizeof(T));
  if (array==NULL) return NULL;
  pnode = l->first;
  for(i=0 ; pnode; i++) {
    pnext = pnode->next;
    array[i] = pnode->e;
    free(pnode);
    pnode = pnext;
  }
  l->size = 0;
  l->first = l->last = NULL;
  return array;
}

#undef EMPTY
#undef LNODE



