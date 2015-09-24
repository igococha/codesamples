#include "list_all.h"

#define T int
#define LT ListInt
#include "list_impl.c"
#undef LT
#undef T

#define T double
#define LT ListDouble
#include "list_impl.c"
#undef LT
#undef T

#define T str
#define LT ListStr
#include "list_impl.c"
#undef LT
#undef T
