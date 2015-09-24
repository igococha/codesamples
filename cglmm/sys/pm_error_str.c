/* error/pm_error_str
* 
* Igor Siveroni
* 
*/

#include "pm_error.h"

const char* pm_errorno_str(const int errno) {
  switch (errno) {
  case PM_SUCCESS:
    return "success";
  case PM_EUNDEFINED:
    return "general error";
  case PM_EARGUMENTS:
    return "incorrent arguments";
  case PM_EOUTOFBOUNDS:
    return "incorrent arguments";
  case PM_EFILE:
    return "file error";
  case PM_EPARSER:
    return "parsing error";
  case PM_EMEM:
    return "out of memory";
  case PM_EMODEL:
    return "model load";
  case PM_EDF:
    return "data frame error";
  case PM_EMATRIXSIZE:
    return "wrong matrix dimensions";
  case PM_ENP:
    return "null pointer";
  default:
    return "unknown error code";
  }

}
