#ifndef PM_MEM_H
#define PM_MEM_H

/* Macros */
#define PM_MEM_ALLOC(nbytes) pm_mem_alloc((nbytes),__FILE__,__FUNCTION__,__LINE__)

#define PM_MEM_CALLOC(n,nbytes) pm_mem_calloc((n),(nbytes),__FILE__,__FUNCTION__,__LINE__)

#define PM_MEM_FREE(p) ((void)(pm_mem_free((p)), (p) = 0))



/* Interface */
void* pm_mem_alloc(long nbytes, const char* file, const char* fn, int line);

void* pm_mem_calloc(long count, long nbytes, const char* file, const char* fn, int line);

void pm_mem_free(void* x);

#endif
