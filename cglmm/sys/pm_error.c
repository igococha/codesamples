#include "pm_error.h"

#include <stdarg.h>
#include <string.h>

#define PM_ERROR_MSGSIZE 256


struct PmError {
  short type;
  char msg[PM_ERROR_MSGSIZE];
  ptrdiff_t last_index;
  int line;
  char function[32];
  char file[32];
};

static PmError the_error;

const char* pm_error_str(PmError *e) {
  return pm_errorno_str(e->type);
}


void pm_error_clear(PmError *e) {
  if (e==NULL) e = &the_error;
  e->type = PM_SUCCESS;
  e->msg[0] = '\0';
  e->last_index = 0;
  e->line = -1;
  e->function[0] = '\0';
  e->file[0] = '\0';
  return;
}

PmError* pm_error_default() {
  return &the_error;
}

void pm_error_print(PmError *e, FILE* fp) {
  if (e==NULL) e = pm_error_default();
  fprintf(fp, "Error: %s\n", pm_errorno_str(e->type));
  if (e->line >= 0) {
    fprintf(fp,"<file: %s> <function: %s> <line: %d>\n",e->file,e->function,e->line);
  }
  if (e->type != PM_SUCCESS && e->last_index > 0) {
    fprintf(fp,"%s\n",e->msg);
  }
  return;
}

int pm_error_get_type(PmError* e) {
  if (e==NULL) e = pm_error_default();
  return e->type;
}

void pm_error_set_type(PmError *e, short error_type) {
  if (e==NULL) e = pm_error_default();
  e->type = error_type;
  return;
}

void pm_error_set_info(PmError *e, const char* file, const char* fn, int line) {
  if (e==NULL) e = pm_error_default();
  e->line = line;
  strcpy(e->file,file);
  strcpy(e->function,fn);
}

int pm_error_set_message(PmError *e, const char* format, ...) {
  va_list vparams;
  int n;
  if (e==NULL) e = pm_error_default();
  va_start(vparams,format); /* vparams start after format param */
  n = vsnprintf(&e->msg[0],PM_ERROR_MSGSIZE,format,vparams);
  if (PM_ERROR_MSGSIZE > n) {
    e->last_index = n;
    return n;
  } else {
    e->last_index = PM_ERROR_MSGSIZE-1;
    return 0;
  }
}

int pm_error_append_message(PmError *e, const char* format, ...) {
  va_list vparams;
  int space,n;
  if (e==NULL) e = pm_error_default();
  space = PM_ERROR_MSGSIZE - e->last_index;
  if (space <= 1) {
    return 0;
  }
  va_start(vparams,format); /* vparams start after format param */
  n = vsnprintf(&e->msg[e->last_index],space,format,vparams);
  if (space>n) {
    e->last_index += n;
    return n;
  } else {
    e->last_index = PM_ERROR_MSGSIZE-1;
    return (space-1);
  }
}
