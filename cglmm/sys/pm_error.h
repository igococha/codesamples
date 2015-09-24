#ifndef PM_ERROR_H
#define PM_ERROR_H

#include <stdio.h>
#include <stddef.h>


/* error messages */
#define PM_SUCCESS 0
#define PM_EUNDEFINED 1
#define PM_EARGUMENTS 2
#define PM_EOUTOFBOUNDS 3
#define PM_EFILE 4
#define PM_EPARSER 5
#define PM_EMEM 6
#define PM_EMODEL 7
#define PM_EDF 8
#define PM_EMATRIXSIZE 9
#define PM_ENP 10


typedef struct PmError PmError;

/* current malloc handling
  if (model==NULL) {
    pm_error_set_type(e,PM_EMEM);
    pm_error_set_info(e,__FILE__,__FUNCTION__,__LINE__);
  }

 */

/* Variadic macros */
#define PM_ERROR(errno,...) do { \
  pm_error_set_type(NULL,errno); \
  pm_error_append_message(NULL,__VA_ARGS__); \
  return errno; \
  } while(0)

#define PM_ERROR_NOMSG(errno,...) do { \
  pm_error_set_type(NULL,errno); \
  return errno; \
  } while(0)

#define PM_ERROR_VAL(errno,val,...) do {	\
  pm_error_set_type(NULL,errno); \
  pm_error_append_message(NULL,__VA_ARGS__); \
  return val; \
  } while(0)

#define PM_ERROR_VAL_NOMSG(errno,val,...) do {	\
  pm_error_set_type(NULL,errno); \
  return val; \
  } while(0)

const char* pm_error_str(PmError* e);
const char* pm_errorno_str(const int errno);

void pm_error_clear(PmError *e);
PmError* pm_error_default();

void pm_error_print(PmError *e, FILE* fp);

int pm_error_get_type(PmError* e);
void pm_error_set_type(PmError *e, short error_type);
void pm_error_set_info(PmError *e, const char* fname, const char* fn, int line);
int pm_error_set_message(PmError *e, const char* format, ...);
int pm_error_append_message(PmError *e, const char* format, ...);

#endif
