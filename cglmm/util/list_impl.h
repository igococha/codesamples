/*
Author: Igor Siveroni
 */

/*
Important: 
Do not add header guards because this file must be able to be included several times.
 */

#include "templates.h"

typedef struct LT *LT;

LT METHOD(list,T,new)();
int METHOD(list,T,size)(LT l);
int METHOD(list,T,empty)(LT l);
int METHOD(list,T,front)(LT l,T* x);
int METHOD(list,T,back)(LT l, T* x);
int METHOD(list,T,pop_front)(LT l,T *x);
int METHOD(list,T,pop_back)(LT l, T* x);
int METHOD(list,T,element)(LT l, int index,T* x);
void METHOD(list,T,push_front)(LT l, T e);
void METHOD(list,T,push_back)(LT l, T e);
void METHOD(list,T,clear)(LT l);
void METHOD(list,T,free)(LT* l);
			
void METHOD(list,T,map)(LT l, void apply(T *x, void *cl), void *cl);
T* METHOD(list,T,toarray)(LT l);
T* METHOD(list,T,clear_toarray)(LT l);


