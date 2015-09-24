#include "pm_mem.h"
#include "pm_error.h"

#include <stdlib.h>
#include <assert.h>



void* pm_mem_alloc(long nbytes, const char* file, const char* fn, int line) {
  void* p;

  assert(nbytes>0);
  p = malloc(nbytes);
  if (p == NULL) { /* report error */
    printf("malloc returned null\n");
    pm_error_clear(NULL);
    pm_error_set_type(NULL,PM_EMEM);
    pm_error_set_info(NULL,file,fn,line);
    pm_error_set_message(NULL,"Couldn't allocate object\n");
  }
  
  return p;
}

void* pm_mem_calloc(long count, long nbytes, const char* file, const char* fn, int line) {
  void* p;

  assert(nbytes>0);
  p = calloc(count,nbytes);
  if (p == NULL) { /* report error */
    printf("calloc returned null\n");
    pm_error_clear(NULL);
    pm_error_set_type(NULL,PM_EMEM);
    pm_error_set_info(NULL,file,fn,line);
    pm_error_set_message(NULL,"Couldn't allocate object\n");
  }
  
  return p;
}

void pm_mem_free(void *p) {
  if (p) free(p);
}
