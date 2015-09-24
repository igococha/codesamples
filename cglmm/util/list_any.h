#ifndef UTIL_LIST_ANY_H
#define UTIL_LIST_ANY_H

#ifdef T
#undef T
#endif

#ifdef LT
#undef LT
#endif

typedef void* any;

#define T any
#define LT ListAny
#include "list_impl.h"
#undef LT
#undef T

#endif
