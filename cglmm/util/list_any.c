#include "list_any.h"

/* 
list_imp.c must be together with the other header files.
list_imp.c must not be compiled. Do not add it to the list of source code. If sources and headers are all mixed up, make use not to use *.c.
 */

#define T any
#define LT ListAny
#include "list_impl.c"
#undef LT
#undef T

