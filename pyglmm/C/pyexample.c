/* pyexample.c */

#include "Python.h"
#include "example.h"

static char py_gcd_doc[] = "Computes GCD of two integers";
static PyObject *
py_gcd (PyObject *self,PyObject *args) {
    int x,y,r;
    if (!PyArg_ParseTuple(args,"ii:gcd",&x,&y)) {
        return NULL;
    }
    r = gcd(x,y);
    return Py_BuildValue("i",r);
}

static PyMethodDef _examplemethods[] = {
    {"gcd",py_gcd,METH_VARARGS, py_gcd_doc},
    {NULL,NULL,0,NULL}
};

/* Python 2 module initialization */
/* Module name = _example */
void init_example(void) {
    PyObject *mod;
    mod = Py_InitModule("_example",_examplemethods);
    PyModule_AddIntMacro(mod,MAGIC);
}