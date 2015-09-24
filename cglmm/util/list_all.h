#ifndef UTIL_LIST_ALL_H
#define UTIL_LIST_ALL_H


#ifdef T
#undef T
#endif

#ifdef LT
#undef LT
#endif

#define T int
#define LT ListInt
#include "list_impl.h"
#undef LT
#undef T

#define T double
#define LT ListDouble
#include "list_impl.h"
#undef LT
#undef T

typedef char* str;

#define T str
#define LT ListStr
#include "list_impl.h"
#undef LT
#undef T


#endif
