#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>
#include "slicot.h"

/* Python wrapper for mb01qd */
static PyObject* py_mb01qd(PyObject* self, PyObject* args) {
    char type;
    i32 m, n, kl, ku, nbl, lda;
    f64 cfrom, cto;
    PyObject *a_obj, *nrows_obj = NULL;
    PyArrayObject *a_array, *nrows_array = NULL;
    i32 info;

    if (!PyArg_ParseTuple(args, "ciiiiddO|O",
                          &type, &m, &n, &kl, &ku, &cfrom, &cto,
                          &a_obj, &nrows_obj)) {
        return NULL;
    }

    /* Convert to NumPy arrays - preserve Fortran-order (column-major) */
    a_array = (PyArrayObject*)PyArray_FROM_OTF(a_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (a_array == NULL) {
        return NULL;
    }

    /* Extract leading dimension from array shape */
    npy_intp *a_dims = PyArray_DIMS(a_array);
    lda = (i32)a_dims[0];

    /* Handle optional nrows parameter */
    i32 *nrows_ptr = NULL;
    if (nrows_obj != NULL && nrows_obj != Py_None) {
        nrows_array = (PyArrayObject*)PyArray_FROM_OTF(nrows_obj, NPY_INT32,
                                                       NPY_ARRAY_IN_ARRAY);
        if (nrows_array == NULL) {
            Py_DECREF(a_array);
            return NULL;
        }
        nrows_ptr = (i32*)PyArray_DATA(nrows_array);
        nbl = (i32)PyArray_SIZE(nrows_array);
    } else {
        nbl = 0;
    }

    /* Call C function */
    f64 *a_data = (f64*)PyArray_DATA(a_array);
    mb01qd(type, m, n, kl, ku, cfrom, cto, nbl, nrows_ptr, a_data, lda, &info);

    /* Clean up and return */
    if (nrows_array != NULL) {
        Py_DECREF(nrows_array);
    }

    PyObject *result = Py_BuildValue("Oi", a_array, info);
    Py_DECREF(a_array);
    return result;
}

/* Python wrapper for mb03oy */
static PyObject* py_mb03oy(PyObject* self, PyObject* args) {
    i32 m, n, lda;
    f64 rcond, svlmax;
    PyObject *a_obj;
    PyArrayObject *a_array;
    i32 rank, info;

    if (!PyArg_ParseTuple(args, "iiOdd", &m, &n, &a_obj, &rcond, &svlmax)) {
        return NULL;
    }

    /* Validate dimensions before allocation */
    if (m < 0 || n < 0) {
        PyErr_Format(PyExc_ValueError, "Dimensions must be non-negative (m=%d, n=%d)", m, n);
        return NULL;
    }

    /* Convert to NumPy array - preserve Fortran-order (column-major) */
    a_array = (PyArrayObject*)PyArray_FROM_OTF(a_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (a_array == NULL) {
        return NULL;
    }

    /* Extract leading dimension (ensure lda >= 1 even if m=0) */
    npy_intp *a_dims = PyArray_DIMS(a_array);
    lda = (i32)a_dims[0];
    if (lda < 1) lda = 1;

    /* Allocate output arrays (handle n=0 edge case) */
    f64 *sval = (f64*)malloc(3 * sizeof(f64));
    i32 *jpvt = (n > 0) ? (i32*)malloc(n * sizeof(i32)) : NULL;
    i32 mn = (m < n) ? m : n;
    f64 *tau = (mn > 0) ? (f64*)malloc(mn * sizeof(f64)) : NULL;
    i32 dwork_size = (n > 0) ? (3*n - 1) : 1;
    f64 *dwork = (f64*)malloc(dwork_size * sizeof(f64));

    if (sval == NULL || dwork == NULL || (n > 0 && jpvt == NULL) || (mn > 0 && tau == NULL)) {
        free(sval); free(jpvt); free(tau); free(dwork);
        Py_DECREF(a_array);
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate work arrays");
        return NULL;
    }

    /* Call C function */
    f64 *a_data = (f64*)PyArray_DATA(a_array);
    mb03oy(m, n, a_data, lda, rcond, svlmax, &rank, sval, jpvt, tau, dwork, &info);

    /* Create output NumPy arrays */
    npy_intp sval_dims[1] = {3};
    npy_intp jpvt_dims[1] = {n > 0 ? n : 0};
    npy_intp tau_dims[1] = {mn};

    PyObject *sval_array = PyArray_SimpleNewFromData(1, sval_dims, NPY_DOUBLE, sval);
    PyObject *jpvt_array = (n > 0) ? PyArray_SimpleNewFromData(1, jpvt_dims, NPY_INT32, jpvt) : PyArray_EMPTY(1, jpvt_dims, NPY_INT32, 0);
    PyObject *tau_array = (mn > 0) ? PyArray_SimpleNewFromData(1, tau_dims, NPY_DOUBLE, tau) : PyArray_EMPTY(1, tau_dims, NPY_DOUBLE, 0);

    /* Transfer ownership to NumPy arrays */
    PyArray_ENABLEFLAGS((PyArrayObject*)sval_array, NPY_ARRAY_OWNDATA);
    if (n > 0) PyArray_ENABLEFLAGS((PyArrayObject*)jpvt_array, NPY_ARRAY_OWNDATA);
    if (mn > 0) PyArray_ENABLEFLAGS((PyArrayObject*)tau_array, NPY_ARRAY_OWNDATA);

    free(dwork);

    /* Resolve writebackifcopy before decref */
    PyArray_ResolveWritebackIfCopy(a_array);

    /* Build result tuple */
    PyObject *result = Py_BuildValue("(OiiOOO)", a_array, rank, info,
                                     sval_array, jpvt_array, tau_array);

    Py_DECREF(a_array);
    Py_DECREF(sval_array);
    Py_DECREF(jpvt_array);
    Py_DECREF(tau_array);

    return result;
}

/* Module method definitions */
static PyMethodDef SlicotMethods[] = {
    {"mb01qd", py_mb01qd, METH_VARARGS,
     "Multiply matrix by scalar CTO/CFROM without overflow/underflow.\n\n"
     "Parameters:\n"
     "  type (str): Matrix storage type ('G', 'L', 'U', 'H', 'B', 'Q', 'Z')\n"
     "  m (int): Number of rows\n"
     "  n (int): Number of columns\n"
     "  kl (int): Lower bandwidth\n"
     "  ku (int): Upper bandwidth\n"
     "  cfrom (float): Denominator scalar\n"
     "  cto (float): Numerator scalar\n"
     "  a (ndarray): Matrix array (column-major, shape (m,n))\n"
     "  nrows (ndarray, optional): Block sizes\n\n"
     "Returns:\n"
     "  (a, info): Modified matrix and exit code\n"},

    {"mb03oy", py_mb03oy, METH_VARARGS,
     "Matrix rank determination by incremental condition estimation.\n\n"
     "Parameters:\n"
     "  m (int): Number of rows\n"
     "  n (int): Number of columns\n"
     "  a (ndarray): Matrix array (column-major, shape (m,n))\n"
     "  rcond (float): Threshold for rank determination\n"
     "  svlmax (float): Estimate of largest singular value\n\n"
     "Returns:\n"
     "  (a, rank, info, sval, jpvt, tau): QR factorization results\n"},

    {NULL, NULL, 0, NULL}
};

/* Module definition */
static struct PyModuleDef slicotmodule = {
    PyModuleDef_HEAD_INIT,
    "_slicot",
    "SLICOT C library Python bindings",
    -1,
    SlicotMethods
};

/* Module initialization */
PyMODINIT_FUNC PyInit__slicot(void) {
    import_array();
    return PyModule_Create(&slicotmodule);
}
