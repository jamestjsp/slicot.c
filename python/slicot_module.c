/*
 * SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 1996-2025, The SLICOT Team (original Fortran77 code)
 * Copyright (c) 2025, slicot.c contributors (C11 translation)
 */

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdbool.h>
#include <stdlib.h>
#include <ctype.h>
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

    PyArray_ResolveWritebackIfCopy(a_array);
    PyObject *result = Py_BuildValue("Oi", a_array, info);
    Py_DECREF(a_array);
    return result;
}

/* Python wrapper for mb01rx */
static PyObject* py_mb01rx(PyObject* self, PyObject* args) {
    const char *side_str, *uplo_str, *trans_str;
    char side, uplo, trans;
    i32 m, n, ldr, lda, ldb;
    f64 alpha, beta;
    PyObject *r_obj, *a_obj, *b_obj;
    PyArrayObject *r_array, *a_array, *b_array;
    i32 info;

    if (!PyArg_ParseTuple(args, "sssiiddOOO",
                          &side_str, &uplo_str, &trans_str, &m, &n, &alpha, &beta,
                          &r_obj, &a_obj, &b_obj)) {
        return NULL;
    }

    side = side_str[0];
    uplo = uplo_str[0];
    trans = trans_str[0];

    /* Convert to NumPy arrays - preserve Fortran-order (column-major) */
    r_array = (PyArrayObject*)PyArray_FROM_OTF(r_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (r_array == NULL) {
        return NULL;
    }

    a_array = (PyArrayObject*)PyArray_FROM_OTF(a_obj, NPY_DOUBLE, NPY_ARRAY_FARRAY);
    if (a_array == NULL) {
        Py_DECREF(r_array);
        return NULL;
    }

    b_array = (PyArrayObject*)PyArray_FROM_OTF(b_obj, NPY_DOUBLE, NPY_ARRAY_FARRAY);
    if (b_array == NULL) {
        Py_DECREF(r_array);
        Py_DECREF(a_array);
        return NULL;
    }

    /* Extract leading dimensions from array shapes */
    npy_intp *r_dims = PyArray_DIMS(r_array);
    npy_intp *a_dims = PyArray_DIMS(a_array);
    npy_intp *b_dims = PyArray_DIMS(b_array);

    ldr = (i32)r_dims[0];
    lda = (i32)a_dims[0];
    ldb = (i32)b_dims[0];

    /* Call C function */
    f64 *r_data = (f64*)PyArray_DATA(r_array);
    const f64 *a_data = (const f64*)PyArray_DATA(a_array);
    const f64 *b_data = (const f64*)PyArray_DATA(b_array);

    info = slicot_mb01rx(side, uplo, trans, m, n, alpha, beta,
                         r_data, ldr, a_data, lda, b_data, ldb);

    /* Clean up and return */
    Py_DECREF(a_array);
    Py_DECREF(b_array);

    PyArray_ResolveWritebackIfCopy(r_array);
    PyObject *result = Py_BuildValue("Oi", r_array, info);
    Py_DECREF(r_array);
    return result;
}

/* Python wrapper for mb01rb */
static PyObject* py_mb01rb(PyObject* self, PyObject* args) {
    const char *side_str, *uplo_str, *trans_str;
    char side, uplo, trans;
    i32 m, n, ldr, lda, ldb;
    f64 alpha, beta;
    PyObject *r_obj, *a_obj, *b_obj;
    PyArrayObject *r_array, *a_array, *b_array;
    i32 info;

    if (!PyArg_ParseTuple(args, "sssiiddOOO",
                          &side_str, &uplo_str, &trans_str, &m, &n, &alpha, &beta,
                          &r_obj, &a_obj, &b_obj)) {
        return NULL;
    }

    side = side_str[0];
    uplo = uplo_str[0];
    trans = trans_str[0];

    r_array = (PyArrayObject*)PyArray_FROM_OTF(r_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (r_array == NULL) {
        return NULL;
    }

    a_array = (PyArrayObject*)PyArray_FROM_OTF(a_obj, NPY_DOUBLE, NPY_ARRAY_FARRAY);
    if (a_array == NULL) {
        Py_DECREF(r_array);
        return NULL;
    }

    b_array = (PyArrayObject*)PyArray_FROM_OTF(b_obj, NPY_DOUBLE, NPY_ARRAY_FARRAY);
    if (b_array == NULL) {
        Py_DECREF(r_array);
        Py_DECREF(a_array);
        return NULL;
    }

    npy_intp *r_dims = PyArray_DIMS(r_array);
    npy_intp *a_dims = PyArray_DIMS(a_array);
    npy_intp *b_dims = PyArray_DIMS(b_array);

    ldr = (i32)r_dims[0];
    lda = (i32)a_dims[0];
    ldb = (i32)b_dims[0];

    f64 *r_data = (f64*)PyArray_DATA(r_array);
    const f64 *a_data = (const f64*)PyArray_DATA(a_array);
    const f64 *b_data = (const f64*)PyArray_DATA(b_array);

    mb01rb(&side, &uplo, &trans, m, n, alpha, beta,
           r_data, ldr, a_data, lda, b_data, ldb, &info);

    Py_DECREF(a_array);
    Py_DECREF(b_array);

    if (info < 0) {
        PyArray_DiscardWritebackIfCopy(r_array);
        Py_DECREF(r_array);
        PyErr_Format(PyExc_ValueError, "Parameter %d had an illegal value", -info);
        return NULL;
    }

    PyArray_ResolveWritebackIfCopy(r_array);
    PyObject *result = Py_BuildValue("Oi", r_array, info);
    Py_DECREF(r_array);
    return result;
}

/* Python wrapper for mb01td */
static PyObject* py_mb01td(PyObject* self, PyObject* args) {
    i32 n, lda, ldb;
    PyObject *a_obj, *b_obj;
    PyArrayObject *a_array, *b_array;
    i32 info;

    if (!PyArg_ParseTuple(args, "OO", &a_obj, &b_obj)) {
        return NULL;
    }

    a_array = (PyArrayObject*)PyArray_FROM_OTF(a_obj, NPY_DOUBLE, NPY_ARRAY_FARRAY);
    if (a_array == NULL) {
        return NULL;
    }

    b_array = (PyArrayObject*)PyArray_FROM_OTF(b_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (b_array == NULL) {
        Py_DECREF(a_array);
        return NULL;
    }

    npy_intp *a_dims = PyArray_DIMS(a_array);

    n = (i32)a_dims[0];
    lda = n > 1 ? n : 1;
    ldb = n > 1 ? n : 1;

    const f64 *a_data = (const f64*)PyArray_DATA(a_array);
    f64 *b_data = (f64*)PyArray_DATA(b_array);

    f64 *dwork = NULL;
    if (n > 1) {
        dwork = (f64*)malloc((n - 1) * sizeof(f64));
        if (dwork == NULL) {
            Py_DECREF(a_array);
            Py_DECREF(b_array);
            return PyErr_NoMemory();
        }
    }

    mb01td(n, a_data, lda, b_data, ldb, dwork, &info);

    free(dwork);
    Py_DECREF(a_array);

    if (info < 0) {
        PyArray_DiscardWritebackIfCopy(b_array);
        Py_DECREF(b_array);
        PyErr_Format(PyExc_ValueError, "Parameter %d had an illegal value", -info);
        return NULL;
    }

    PyArray_ResolveWritebackIfCopy(b_array);
    PyObject *result = Py_BuildValue("Oi", b_array, info);
    Py_DECREF(b_array);
    return result;
}

/* Python wrapper for mb04od */
static PyObject* py_mb04od(PyObject* self, PyObject* args) {
    const char *uplo_str;
    char uplo;
    i32 n, m, p, ldr, lda, ldb, ldc;
    PyObject *r_obj, *a_obj, *b_obj, *c_obj;
    PyArrayObject *r_array, *a_array, *b_array, *c_array, *tau_array;
    f64 *dwork;
    i32 ldwork;

    if (!PyArg_ParseTuple(args, "siiiOOOO",
                          &uplo_str, &n, &m, &p,
                          &r_obj, &a_obj, &b_obj, &c_obj)) {
        return NULL;
    }

    uplo = uplo_str[0];

    r_array = (PyArrayObject*)PyArray_FROM_OTF(r_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (r_array == NULL) {
        return NULL;
    }

    a_array = (PyArrayObject*)PyArray_FROM_OTF(a_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (a_array == NULL) {
        Py_DECREF(r_array);
        return NULL;
    }

    b_array = (PyArrayObject*)PyArray_FROM_OTF(b_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (b_array == NULL) {
        Py_DECREF(r_array);
        Py_DECREF(a_array);
        return NULL;
    }

    c_array = (PyArrayObject*)PyArray_FROM_OTF(c_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (c_array == NULL) {
        Py_DECREF(r_array);
        Py_DECREF(a_array);
        Py_DECREF(b_array);
        return NULL;
    }

    npy_intp *r_dims = PyArray_DIMS(r_array);
    npy_intp *a_dims = PyArray_DIMS(a_array);
    npy_intp *b_dims = PyArray_DIMS(b_array);
    npy_intp *c_dims = PyArray_DIMS(c_array);

    ldr = (i32)r_dims[0];
    lda = (i32)a_dims[0];
    ldb = (i32)b_dims[0];
    ldc = (i32)c_dims[0];

    npy_intp tau_dims[1] = {n > 0 ? n : 1};
    tau_array = (PyArrayObject*)PyArray_SimpleNew(1, tau_dims, NPY_DOUBLE);
    if (tau_array == NULL) {
        Py_DECREF(r_array);
        Py_DECREF(a_array);
        Py_DECREF(b_array);
        Py_DECREF(c_array);
        return NULL;
    }

    ldwork = (n - 1) > m ? (n - 1) : m;
    if (ldwork < 1) ldwork = 1;
    dwork = (f64*)malloc(ldwork * sizeof(f64));
    if (dwork == NULL) {
        Py_DECREF(r_array);
        Py_DECREF(a_array);
        Py_DECREF(b_array);
        Py_DECREF(c_array);
        Py_DECREF(tau_array);
        PyErr_NoMemory();
        return NULL;
    }

    f64 *r_data = (f64*)PyArray_DATA(r_array);
    f64 *a_data = (f64*)PyArray_DATA(a_array);
    f64 *b_data = (f64*)PyArray_DATA(b_array);
    f64 *c_data = (f64*)PyArray_DATA(c_array);
    f64 *tau_data = (f64*)PyArray_DATA(tau_array);

    mb04od(&uplo, n, m, p, r_data, ldr, a_data, lda,
           b_data, ldb, c_data, ldc, tau_data, dwork);

    free(dwork);

    PyArray_ResolveWritebackIfCopy(r_array);
    PyArray_ResolveWritebackIfCopy(a_array);
    PyArray_ResolveWritebackIfCopy(b_array);
    PyArray_ResolveWritebackIfCopy(c_array);
    PyObject *result = Py_BuildValue("OOOOO", r_array, a_array, b_array, c_array, tau_array);
    Py_DECREF(r_array);
    Py_DECREF(a_array);
    Py_DECREF(b_array);
    Py_DECREF(c_array);
    Py_DECREF(tau_array);
    return result;
}

/* Python wrapper for mb03oy */
static PyObject* py_mb03oy(PyObject* self, PyObject* args) {
    i32 m, n, lda;
    f64 rcond, svlmax;
    PyObject *a_obj;
    PyArrayObject *a_array;
    i32 rank = 0, info = 0;

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

/* Python wrapper for mb03od */
static PyObject* py_mb03od(PyObject* self, PyObject* args, PyObject* kwargs) {
    const char *jobqr = "Q";
    i32 m, n, lda;
    f64 rcond, svlmax;
    PyObject *a_obj;
    PyArrayObject *a_array;
    i32 rank = 0, info = 0;

    static char *kwlist[] = {"m", "n", "a", "rcond", "svlmax", "jobqr", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "iiOdd|s", kwlist,
                                     &m, &n, &a_obj, &rcond, &svlmax, &jobqr)) {
        return NULL;
    }

    if (m < 0 || n < 0) {
        PyErr_Format(PyExc_ValueError, "Dimensions must be non-negative (m=%d, n=%d)", m, n);
        return NULL;
    }

    a_array = (PyArrayObject*)PyArray_FROM_OTF(a_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (a_array == NULL) {
        return NULL;
    }

    npy_intp *a_dims = PyArray_DIMS(a_array);
    lda = (i32)a_dims[0];
    if (lda < 1) lda = 1;

    i32 mn = (m < n) ? m : n;
    i32 *jpvt = (n > 0) ? (i32*)calloc(n, sizeof(i32)) : NULL;
    f64 *tau = (mn > 0) ? (f64*)malloc(mn * sizeof(f64)) : NULL;
    f64 *sval = (f64*)malloc(3 * sizeof(f64));

    bool ljobqr = (*jobqr == 'Q' || *jobqr == 'q');
    i32 ldwork = ljobqr ? (3*n + 1) : ((2*mn > 1) ? 2*mn : 1);
    f64 *dwork = (f64*)malloc(ldwork * sizeof(f64));

    if (sval == NULL || dwork == NULL || (n > 0 && jpvt == NULL) || (mn > 0 && tau == NULL)) {
        free(sval); free(jpvt); free(tau); free(dwork);
        Py_DECREF(a_array);
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate work arrays");
        return NULL;
    }

    f64 *a_data = (f64*)PyArray_DATA(a_array);
    mb03od(jobqr, m, n, a_data, lda, jpvt, rcond, svlmax, tau, &rank, sval, dwork, ldwork, &info);

    npy_intp sval_dims[1] = {3};
    npy_intp jpvt_dims[1] = {n > 0 ? n : 0};

    PyObject *sval_array = PyArray_SimpleNewFromData(1, sval_dims, NPY_DOUBLE, sval);
    PyObject *jpvt_array = (n > 0) ? PyArray_SimpleNewFromData(1, jpvt_dims, NPY_INT32, jpvt) : PyArray_EMPTY(1, jpvt_dims, NPY_INT32, 0);

    PyArray_ENABLEFLAGS((PyArrayObject*)sval_array, NPY_ARRAY_OWNDATA);
    if (n > 0) PyArray_ENABLEFLAGS((PyArrayObject*)jpvt_array, NPY_ARRAY_OWNDATA);

    free(tau);
    free(dwork);

    PyArray_ResolveWritebackIfCopy(a_array);

    PyObject *result = Py_BuildValue("OiOi", jpvt_array, rank, sval_array, info);

    Py_DECREF(a_array);
    Py_DECREF(sval_array);
    Py_DECREF(jpvt_array);

    return result;
}

static PyObject* py_tg01fd(PyObject* self, PyObject* args) {
    const char *compq, *compz, *joba;
    i32 l, n, m, p;
    PyObject *a_obj, *e_obj, *b_obj, *c_obj;
    f64 tol;
    PyArrayObject *a_array, *e_array, *b_array, *c_array;
    i32 ranke = 0, rnka22 = 0, info = 0;
    i32 lda, lde, ldb, ldc, ldq, ldz;

    if (!PyArg_ParseTuple(args, "sssiiiiOOOOd", &compq, &compz, &joba,
                          &l, &n, &m, &p, &a_obj, &e_obj, &b_obj, &c_obj, &tol)) {
        return NULL;
    }

    if (l < 0 || n < 0 || m < 0 || p < 0) {
        PyErr_Format(PyExc_ValueError, "Dimensions must be non-negative");
        return NULL;
    }

    a_array = (PyArrayObject*)PyArray_FROM_OTF(a_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    e_array = (PyArrayObject*)PyArray_FROM_OTF(e_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    b_array = (PyArrayObject*)PyArray_FROM_OTF(b_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    c_array = (PyArrayObject*)PyArray_FROM_OTF(c_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);

    if (a_array == NULL || e_array == NULL || b_array == NULL || c_array == NULL) {
        Py_XDECREF(a_array);
        Py_XDECREF(e_array);
        Py_XDECREF(b_array);
        Py_XDECREF(c_array);
        return NULL;
    }

    npy_intp *a_dims = PyArray_DIMS(a_array);
    npy_intp *e_dims = PyArray_DIMS(e_array);
    npy_intp *b_dims = PyArray_DIMS(b_array);
    npy_intp *c_dims = PyArray_DIMS(c_array);

    lda = (l > 0) ? (i32)a_dims[0] : 1;
    lde = (l > 0) ? (i32)e_dims[0] : 1;
    ldb = (l > 0) ? (i32)b_dims[0] : 1;
    ldc = (p > 0) ? (i32)c_dims[0] : 1;
    ldq = (l > 0) ? l : 1;
    ldz = (n > 0) ? n : 1;

    i32 ln = (l < n) ? l : n;
    i32 temp1 = n + p;
    i32 temp2 = (3 * n - 1 > m) ? 3 * n - 1 : m;
    temp2 = (temp2 > l) ? temp2 : l;
    temp2 = ln + temp2;
    i32 ldwork = (temp1 > temp2) ? temp1 : temp2;
    ldwork = (ldwork > 1) ? ldwork : 1;

    i32 *iwork = (n > 0) ? (i32*)malloc(n * sizeof(i32)) : NULL;
    f64 *dwork = (f64*)malloc(ldwork * sizeof(f64));
    f64 *q = NULL;
    f64 *z = NULL;

    bool compq_needed = (compq[0] == 'I' || compq[0] == 'i' || compq[0] == 'U' || compq[0] == 'u');
    bool compz_needed = (compz[0] == 'I' || compz[0] == 'i' || compz[0] == 'U' || compz[0] == 'u');

    if (compq_needed && l > 0) {
        q = (f64*)calloc(l * l, sizeof(f64));
    }
    if (compz_needed && n > 0) {
        z = (f64*)calloc(n * n, sizeof(f64));
    }

    if (dwork == NULL || (n > 0 && iwork == NULL) ||
        (compq_needed && l > 0 && q == NULL) ||
        (compz_needed && n > 0 && z == NULL)) {
        free(iwork);
        free(dwork);
        free(q);
        free(z);
        Py_DECREF(a_array);
        Py_DECREF(e_array);
        Py_DECREF(b_array);
        Py_DECREF(c_array);
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate work arrays");
        return NULL;
    }

    f64 *a_data = (f64*)PyArray_DATA(a_array);
    f64 *e_data = (f64*)PyArray_DATA(e_array);
    f64 *b_data = (f64*)PyArray_DATA(b_array);
    f64 *c_data = (f64*)PyArray_DATA(c_array);

    tg01fd(compq, compz, joba, l, n, m, p, a_data, lda, e_data, lde,
           b_data, ldb, c_data, ldc, q, ldq, z, ldz, &ranke, &rnka22,
           tol, iwork, dwork, ldwork, &info);

    npy_intp q_dims[2] = {l, l};
    npy_intp z_dims[2] = {n, n};
    npy_intp q_strides[2] = {sizeof(f64), l * sizeof(f64)};
    npy_intp z_strides[2] = {sizeof(f64), n * sizeof(f64)};

    PyObject *q_array, *z_array;
    if (compq_needed && l > 0) {
        q_array = PyArray_New(&PyArray_Type, 2, q_dims, NPY_DOUBLE, q_strides, q, 0, NPY_ARRAY_FARRAY, NULL);
        PyArray_ENABLEFLAGS((PyArrayObject*)q_array, NPY_ARRAY_OWNDATA);
    } else {
        q_array = PyArray_EMPTY(2, q_dims, NPY_DOUBLE, 1);
    }

    if (compz_needed && n > 0) {
        z_array = PyArray_New(&PyArray_Type, 2, z_dims, NPY_DOUBLE, z_strides, z, 0, NPY_ARRAY_FARRAY, NULL);
        PyArray_ENABLEFLAGS((PyArrayObject*)z_array, NPY_ARRAY_OWNDATA);
    } else {
        z_array = PyArray_EMPTY(2, z_dims, NPY_DOUBLE, 1);
    }

    free(iwork);
    free(dwork);

    PyArray_ResolveWritebackIfCopy(a_array);
    PyArray_ResolveWritebackIfCopy(e_array);
    PyArray_ResolveWritebackIfCopy(b_array);
    PyArray_ResolveWritebackIfCopy(c_array);

    PyObject *result = Py_BuildValue("(OOOOOOiii)", a_array, e_array, b_array, c_array,
                                     q_array, z_array, ranke, rnka22, info);

    Py_DECREF(a_array);
    Py_DECREF(e_array);
    Py_DECREF(b_array);
    Py_DECREF(c_array);
    Py_DECREF(q_array);
    Py_DECREF(z_array);

    return result;
}

typedef struct {
    PyObject *fcn_callable;
    PyObject *jac_callable;
    i32 m;
    i32 n;
} md03bd_callback_data;

static md03bd_callback_data *g_cb_data = NULL;

static void md03bd_fcn_wrapper(
    i32* iflag, i32 m, i32 n, i32* ipar, i32 lipar,
    const f64* dpar1, i32 ldpar1, const f64* dpar2, i32 ldpar2,
    const f64* x, i32* nfevl, f64* e, f64* j, i32* ldj,
    f64* dwork, i32 ldwork, i32* info
)
{
    *info = 0;
    *nfevl = 0;

    if (*iflag == 3) {
        *ldj = m;
        ipar[0] = m * n;
        ipar[1] = 0;
        ipar[2] = 0;
        ipar[3] = 4*n + 1;
        ipar[4] = 4*n;
        return;
    }

    if (*iflag == 0 || g_cb_data == NULL) {
        return;
    }

    npy_intp x_dims[1] = {n};
    PyObject *x_array = PyArray_SimpleNewFromData(1, x_dims, NPY_DOUBLE, (void*)x);
    if (x_array == NULL) {
        *info = -1;
        return;
    }
    PyArray_CLEARFLAGS((PyArrayObject*)x_array, NPY_ARRAY_WRITEABLE);

    if (*iflag == 1) {
        PyObject *result = PyObject_CallFunctionObjArgs(g_cb_data->fcn_callable, x_array, NULL);
        Py_DECREF(x_array);

        if (result == NULL) {
            PyErr_Print();
            *info = -1;
            return;
        }

        PyArrayObject *e_result = (PyArrayObject*)PyArray_FROM_OTF(result, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
        Py_DECREF(result);

        if (e_result == NULL || PyArray_SIZE(e_result) != m) {
            Py_XDECREF(e_result);
            *info = -1;
            return;
        }

        f64 *e_data = (f64*)PyArray_DATA(e_result);
        for (i32 i = 0; i < m; i++) {
            e[i] = e_data[i];
        }

        Py_DECREF(e_result);

    } else if (*iflag == 2) {
        PyObject *result = PyObject_CallFunctionObjArgs(g_cb_data->jac_callable, x_array, NULL);
        Py_DECREF(x_array);

        if (result == NULL) {
            PyErr_Print();
            *info = -1;
            return;
        }

        PyArrayObject *j_result = (PyArrayObject*)PyArray_FROM_OTF(result, NPY_DOUBLE, NPY_ARRAY_FARRAY);
        Py_DECREF(result);

        if (j_result == NULL || PyArray_DIM(j_result, 0) != m || PyArray_DIM(j_result, 1) != n) {
            Py_XDECREF(j_result);
            *info = -1;
            return;
        }

        f64 *j_data = (f64*)PyArray_DATA(j_result);
        for (i32 col = 0; col < n; col++) {
            for (i32 row = 0; row < m; row++) {
                j[row + col * m] = j_data[row + col * m];
            }
        }

        Py_DECREF(j_result);
    }
}

static void md03bd_qrfact_wrapper(
    i32 n, const i32* ipar, i32 lipar, f64 fnorm,
    f64* j, i32* ldj, f64* e, f64* jnorms, f64* gnorm,
    i32* ipvt, f64* dwork, i32 ldwork, i32* info
)
{
    i32 m = g_cb_data->m;
    md03bx(m, n, fnorm, j, ldj, e, jnorms, gnorm, ipvt, dwork, ldwork, info);
}

static void md03bd_lmparm_wrapper(
    const char* cond, i32 n, const i32* ipar, i32 lipar,
    f64* r, i32 ldr, const i32* ipvt, const f64* diag,
    const f64* qtb, f64 delta, f64* par, i32* rank,
    f64* x, f64* rx, f64 tol, f64* dwork, i32 ldwork, i32* info
)
{
    md03by(cond, n, r, ldr, ipvt, diag, qtb, delta, par, rank, x, rx, tol, dwork, ldwork, info);
}

static PyObject* py_md03bd(PyObject* self, PyObject* args, PyObject* kwargs) {
    i32 m, n, itmax = 100;
    f64 ftol = -1.0, xtol = -1.0, gtol = -1.0;
    PyObject *x_obj, *fcn_obj, *jac_obj;

    static char *kwlist[] = {"m", "n", "x", "fcn", "jac", "itmax", "ftol", "xtol", "gtol", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "iiOOO|iddd", kwlist,
                                     &m, &n, &x_obj, &fcn_obj, &jac_obj,
                                     &itmax, &ftol, &xtol, &gtol)) {
        return NULL;
    }

    if (m < 0 || n < 0 || n > m) {
        PyErr_SetString(PyExc_ValueError, "Invalid dimensions: m >= n >= 0 required");
        return NULL;
    }

    if (!PyCallable_Check(fcn_obj) || !PyCallable_Check(jac_obj)) {
        PyErr_SetString(PyExc_TypeError, "fcn and jac must be callable");
        return NULL;
    }

    PyArrayObject *x_array = (PyArrayObject*)PyArray_FROM_OTF(x_obj, NPY_DOUBLE,
                                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (x_array == NULL) {
        return NULL;
    }

    if (PyArray_SIZE(x_array) != n) {
        PyErr_SetString(PyExc_ValueError, "x must have length n");
        Py_DECREF(x_array);
        return NULL;
    }

    f64 *x_data = (f64*)PyArray_DATA(x_array);
    f64 *diag = (f64*)calloc(n, sizeof(f64));
    i32 *iwork = (i32*)calloc(n + 1, sizeof(i32));

    i32 sizej = m * n;
    i32 lfcn1 = 0;
    i32 lfcn2 = 0;
    i32 lqrf = 4*n + 1;
    i32 llmp = 4*n;

    i32 max1 = (lfcn1 > lfcn2) ? lfcn1 : lfcn2;
    max1 = (max1 > (n + lqrf)) ? max1 : (n + lqrf);
    i32 max2 = (m + lfcn1 > n + llmp) ? (m + lfcn1) : (n + llmp);
    i32 max3 = (n*n + n + max2 > sizej + max1) ? (n*n + n + max2) : (sizej + max1);
    i32 ldwork = (m + max3 > 4) ? (m + max3) : 4;
    f64 *dwork = (f64*)malloc(ldwork * sizeof(f64));

    if (diag == NULL || iwork == NULL || dwork == NULL) {
        free(diag);
        free(iwork);
        free(dwork);
        Py_DECREF(x_array);
        return PyErr_NoMemory();
    }

    md03bd_callback_data cb_data;
    cb_data.fcn_callable = fcn_obj;
    cb_data.jac_callable = jac_obj;
    cb_data.m = m;
    cb_data.n = n;

    g_cb_data = &cb_data;

    i32 lipar = 5;
    i32 ipar_storage[5];
    for (i32 k = 0; k < 5; k++) ipar_storage[k] = 0;

    i32 nfev, njev, iwarn, info;
    const char xinit = 'G';
    const char scale = 'I';
    const char cond = 'N';
    f64 factor = 100.0;
    i32 nprint = 0;
    f64 tol_rank = -1.0;

    md03bd(&xinit, &scale, &cond,
           md03bd_fcn_wrapper, md03bd_qrfact_wrapper, md03bd_lmparm_wrapper,
           m, n, itmax, factor, nprint,
           ipar_storage, lipar,
           NULL, 0,
           NULL, 0,
           x_data, diag, &nfev, &njev,
           ftol, xtol, gtol, tol_rank,
           iwork, dwork, ldwork, &iwarn, &info);

    g_cb_data = NULL;
    f64 fnorm = dwork[1];

    npy_intp x_dims[1] = {n};
    PyObject *x_out = PyArray_EMPTY(1, x_dims, NPY_DOUBLE, 0);
    if (x_out == NULL) {
        free(diag);
        free(iwork);
        free(dwork);
        Py_DECREF(x_array);
        return NULL;
    }

    f64 *x_out_data = (f64*)PyArray_DATA((PyArrayObject*)x_out);
    for (i32 i = 0; i < n; i++) {
        x_out_data[i] = x_data[i];
    }

    free(diag);
    free(iwork);
    free(dwork);

    PyArray_ResolveWritebackIfCopy(x_array);
    PyObject *result = Py_BuildValue("(Oiidii)", x_out, nfev, njev, fnorm, iwarn, info);
    Py_DECREF(x_array);
    Py_DECREF(x_out);

    return result;
}

/* Python wrapper for mb01uy */
static PyObject* py_mb01uy(PyObject* self, PyObject* args) {
    char *side, *uplo, *trans;
    i32 m, n;
    f64 alpha;
    PyObject *t_obj, *a_obj;
    PyArrayObject *t_array, *a_array;
    i32 info;

    if (!PyArg_ParseTuple(args, "sssiidOO", &side, &uplo, &trans, &m, &n, &alpha, &t_obj, &a_obj)) {
        return NULL;
    }

    t_array = (PyArrayObject*)PyArray_FROM_OTF(t_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
    if (t_array == NULL) return NULL;

    a_array = (PyArrayObject*)PyArray_FROM_OTF(a_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
    if (a_array == NULL) {
        Py_DECREF(t_array);
        return NULL;
    }

    i32 ldt = (m > n) ? m : n;
    if ((*side == 'R' || *side == 'r') && ldt < n) ldt = n;
    i32 lda = (i32)PyArray_DIM(a_array, 0);

    f64 *t_in = (f64*)PyArray_DATA(t_array);
    f64 *a_data = (f64*)PyArray_DATA(a_array);

    npy_intp t_dims[2] = {ldt, (ldt > n) ? ldt : n};
    npy_intp t_strides[2] = {sizeof(f64), ldt * sizeof(f64)};
    PyObject *t_out = PyArray_New(&PyArray_Type, 2, t_dims, NPY_DOUBLE, t_strides,
                                   NULL, 0, NPY_ARRAY_FARRAY, NULL);
    if (t_out == NULL) {
        Py_DECREF(t_array);
        Py_DECREF(a_array);
        return NULL;
    }

    f64 *t_data = (f64*)PyArray_DATA((PyArrayObject*)t_out);
    i32 k = (*side == 'L' || *side == 'l') ? m : n;
    for (i32 j = 0; j < k; j++) {
        for (i32 i = 0; i < k; i++) {
            t_data[i + j * ldt] = t_in[i + j * PyArray_DIM(t_array, 0)];
        }
    }

    i32 ldwork = m * n;
    f64 *dwork = (f64*)malloc(ldwork * sizeof(f64));
    if (dwork == NULL) {
        Py_DECREF(t_array);
        Py_DECREF(a_array);
        Py_DECREF(t_out);
        return PyErr_NoMemory();
    }

    mb01uy(side, uplo, trans, m, n, alpha, t_data, ldt, a_data, lda, dwork, ldwork, &info);

    free(dwork);

    npy_intp result_dims[2] = {m, n};
    npy_intp result_strides[2] = {sizeof(f64), m * sizeof(f64)};
    PyObject *result_array = PyArray_New(&PyArray_Type, 2, result_dims, NPY_DOUBLE, result_strides,
                                          NULL, 0, NPY_ARRAY_FARRAY, NULL);
    if (result_array == NULL) {
        Py_DECREF(t_array);
        Py_DECREF(a_array);
        Py_DECREF(t_out);
        return NULL;
    }

    f64 *result_data = (f64*)PyArray_DATA((PyArrayObject*)result_array);
    for (i32 j = 0; j < n; j++) {
        for (i32 i = 0; i < m; i++) {
            result_data[i + j * m] = t_data[i + j * ldt];
        }
    }

    PyObject *result = Py_BuildValue("Oi", result_array, info);
    Py_DECREF(t_array);
    Py_DECREF(a_array);
    Py_DECREF(t_out);
    Py_DECREF(result_array);

    return result;
}

/* Python wrapper for mb02yd */
static PyObject* py_mb02yd(PyObject* self, PyObject* args) {
    char* cond;
    i32 n, rank_in;
    f64 tol;
    PyObject *r_obj, *ipvt_obj, *diag_obj, *qtb_obj;
    PyArrayObject *r_array, *ipvt_array, *diag_array, *qtb_array;
    i32 info;

    if (!PyArg_ParseTuple(args, "siOOOOid", &cond, &n, &r_obj, &ipvt_obj, &diag_obj, &qtb_obj, &rank_in, &tol)) {
        return NULL;
    }

    if (n < 0) {
        PyErr_SetString(PyExc_ValueError, "n must be non-negative");
        return NULL;
    }

    r_array = (PyArrayObject*)PyArray_FROM_OTF(r_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (r_array == NULL) return NULL;

    ipvt_array = (PyArrayObject*)PyArray_FROM_OTF(ipvt_obj, NPY_INT32, NPY_ARRAY_IN_FARRAY);
    if (ipvt_array == NULL) {
        Py_DECREF(r_array);
        return NULL;
    }

    diag_array = (PyArrayObject*)PyArray_FROM_OTF(diag_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
    if (diag_array == NULL) {
        Py_DECREF(r_array);
        Py_DECREF(ipvt_array);
        return NULL;
    }

    qtb_array = (PyArrayObject*)PyArray_FROM_OTF(qtb_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
    if (qtb_array == NULL) {
        Py_DECREF(r_array);
        Py_DECREF(ipvt_array);
        Py_DECREF(diag_array);
        return NULL;
    }

    i32 ldr = (i32)PyArray_DIM(r_array, 0);
    f64* r_data = (f64*)PyArray_DATA(r_array);
    i32* ipvt_data = (i32*)PyArray_DATA(ipvt_array);
    f64* diag_data = (f64*)PyArray_DATA(diag_array);
    f64* qtb_data = (f64*)PyArray_DATA(qtb_array);

    bool econd = (*cond == 'E' || *cond == 'e');
    i32 ldwork = econd ? 4*n : 2*n;
    if (ldwork < 1) ldwork = 1;  /* Ensure at least 1 element for malloc */
    f64* dwork = (f64*)malloc(ldwork * sizeof(f64));
    if (dwork == NULL) {
        Py_DECREF(r_array);
        Py_DECREF(ipvt_array);
        Py_DECREF(diag_array);
        Py_DECREF(qtb_array);
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate workspace");
        return NULL;
    }

    /* Allocate output array x */
    f64* x_data = (n > 0) ? (f64*)calloc(n, sizeof(f64)) : NULL;
    if (n > 0 && x_data == NULL) {
        free(dwork);
        Py_DECREF(r_array);
        Py_DECREF(ipvt_array);
        Py_DECREF(diag_array);
        Py_DECREF(qtb_array);
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate output array");
        return NULL;
    }

    i32 rank = rank_in;

    mb02yd(cond, n, r_data, ldr, ipvt_data, diag_data, qtb_data, &rank, x_data, tol, dwork, ldwork, &info);

    free(dwork);

    /* Resolve writebackifcopy before decref */
    PyArray_ResolveWritebackIfCopy(r_array);

    /* Create NumPy array from x_data */
    npy_intp x_dims[1] = {n > 0 ? n : 0};
    PyObject* x_array = (n > 0) ? PyArray_SimpleNewFromData(1, x_dims, NPY_DOUBLE, x_data) 
                                 : PyArray_EMPTY(1, x_dims, NPY_DOUBLE, 0);
    if (x_array == NULL) {
        free(x_data);
        Py_DECREF(r_array);
        Py_DECREF(ipvt_array);
        Py_DECREF(diag_array);
        Py_DECREF(qtb_array);
        PyErr_SetString(PyExc_MemoryError, "Failed to create output array");
        return NULL;
    }

    /* Transfer ownership to NumPy */
    if (n > 0) PyArray_ENABLEFLAGS((PyArrayObject*)x_array, NPY_ARRAY_OWNDATA);

    if (info < 0) {
        Py_DECREF(r_array);
        Py_DECREF(ipvt_array);
        Py_DECREF(diag_array);
        Py_DECREF(qtb_array);
        Py_DECREF(x_array);
        PyErr_Format(PyExc_ValueError, "mb02yd: parameter %d is invalid", -info);
        return NULL;
    }

    PyObject* result = Py_BuildValue("Oii", x_array, rank, info);

    Py_DECREF(r_array);
    Py_DECREF(ipvt_array);
    Py_DECREF(diag_array);
    Py_DECREF(qtb_array);
    Py_DECREF(x_array);

    return result;
}

static PyObject* py_mb02ud(PyObject* self, PyObject* args, PyObject* kwargs) {
    static char* kwlist[] = {"fact", "side", "trans", "jobp", "m", "n", "alpha", "rcond",
                             "r", "b", "q", "sv", "rank", "rp", "ldwork", NULL};

    char *fact_str, *side_str, *trans_str, *jobp_str;
    i32 m_in, n_in, rank_in = 0;
    f64 alpha, rcond;
    i32 ldwork = 0;
    PyObject *r_obj, *b_obj;
    PyObject *q_obj = NULL, *sv_obj = NULL, *rp_obj = NULL;
    PyArrayObject *r_array, *b_array;
    PyArrayObject *q_array = NULL, *sv_array = NULL;
    i32 info;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ssssiiddOO|OOiOi", kwlist,
                                     &fact_str, &side_str, &trans_str, &jobp_str,
                                     &m_in, &n_in, &alpha, &rcond, &r_obj, &b_obj,
                                     &q_obj, &sv_obj, &rank_in, &rp_obj, &ldwork)) {
        return NULL;
    }

    char fact = fact_str[0];
    char side = side_str[0];
    char jobp = jobp_str[0];

    bool nfct = (fact == 'N' || fact == 'n');
    bool left = (side == 'L' || side == 'l');
    bool pinv = (jobp == 'P' || jobp == 'p');
    (void)rp_obj;

    i32 l = left ? m_in : n_in;

    if (m_in < 0) {
        PyErr_SetString(PyExc_ValueError, "m must be non-negative");
        return NULL;
    }
    if (n_in < 0) {
        PyErr_SetString(PyExc_ValueError, "n must be non-negative");
        return NULL;
    }

    r_array = (PyArrayObject*)PyArray_FROM_OTF(r_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (r_array == NULL) return NULL;

    b_array = (PyArrayObject*)PyArray_FROM_OTF(b_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (b_array == NULL) {
        Py_DECREF(r_array);
        return NULL;
    }

    i32 ldr = (i32)PyArray_DIM(r_array, 0);
    i32 ldb = (i32)PyArray_DIM(b_array, 0);

    f64* r_data = (f64*)PyArray_DATA(r_array);
    f64* b_data = (f64*)PyArray_DATA(b_array);

    f64* q_data = NULL;
    f64* sv_data = NULL;
    i32 ldq = l > 0 ? l : 1;

    if (nfct) {
        q_data = (f64*)calloc(l > 0 ? l * l : 1, sizeof(f64));
        sv_data = (f64*)calloc(l > 0 ? l : 1, sizeof(f64));
        if ((l > 0 && q_data == NULL) || (l > 0 && sv_data == NULL)) {
            free(q_data);
            free(sv_data);
            Py_DECREF(r_array);
            Py_DECREF(b_array);
            PyErr_SetString(PyExc_MemoryError, "Failed to allocate Q or SV");
            return NULL;
        }
    } else {
        if (q_obj == NULL || sv_obj == NULL) {
            Py_DECREF(r_array);
            Py_DECREF(b_array);
            PyErr_SetString(PyExc_ValueError, "q and sv are required when fact='F'");
            return NULL;
        }
        q_array = (PyArrayObject*)PyArray_FROM_OTF(q_obj, NPY_DOUBLE,
                                                   NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
        if (q_array == NULL) {
            Py_DECREF(r_array);
            Py_DECREF(b_array);
            return NULL;
        }
        sv_array = (PyArrayObject*)PyArray_FROM_OTF(sv_obj, NPY_DOUBLE,
                                                    NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
        if (sv_array == NULL) {
            Py_DECREF(r_array);
            Py_DECREF(b_array);
            Py_DECREF(q_array);
            return NULL;
        }
        q_data = (f64*)PyArray_DATA(q_array);
        sv_data = (f64*)PyArray_DATA(sv_array);
        ldq = (i32)PyArray_DIM(q_array, 0);
    }

    f64* rp_data = NULL;
    i32 ldrp = 1;
    if (pinv) {
        ldrp = l > 0 ? l : 1;
        rp_data = (f64*)calloc(ldrp * l > 0 ? ldrp * l : 1, sizeof(f64));
        if (l > 0 && rp_data == NULL) {
            if (nfct) { free(q_data); free(sv_data); }
            else { Py_XDECREF(q_array); Py_XDECREF(sv_array); }
            Py_DECREF(r_array);
            Py_DECREF(b_array);
            PyErr_SetString(PyExc_MemoryError, "Failed to allocate RP");
            return NULL;
        }
    }

    i32 minwrk = nfct ? (5 * l > 1 ? 5 * l : 1) : (l > 1 ? l : 1);
    i32 mn = m_in * n_in;
    i32 optwork = minwrk > mn ? minwrk : mn;
    if (ldwork == 0) ldwork = optwork;
    if (ldwork < minwrk) ldwork = minwrk;

    f64* dwork = (f64*)calloc(ldwork > 0 ? ldwork : 1, sizeof(f64));
    if (dwork == NULL) {
        if (pinv) free(rp_data);
        if (nfct) { free(q_data); free(sv_data); }
        else { Py_XDECREF(q_array); Py_XDECREF(sv_array); }
        Py_DECREF(r_array);
        Py_DECREF(b_array);
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate workspace");
        return NULL;
    }

    i32 rank = nfct ? 0 : rank_in;

    mb02ud(fact_str, side_str, trans_str, jobp_str, m_in, n_in, alpha, rcond,
           &rank, r_data, ldr, q_data, ldq, sv_data, b_data, ldb,
           rp_data, ldrp, dwork, ldwork, &info);

    free(dwork);

    PyArray_ResolveWritebackIfCopy(r_array);
    PyArray_ResolveWritebackIfCopy(b_array);
    if (!nfct) {
        PyArray_ResolveWritebackIfCopy(q_array);
        PyArray_ResolveWritebackIfCopy(sv_array);
    }

    if (info < 0) {
        if (pinv) free(rp_data);
        if (nfct) { free(q_data); free(sv_data); }
        else { Py_XDECREF(q_array); Py_XDECREF(sv_array); }
        Py_DECREF(r_array);
        Py_DECREF(b_array);
        PyErr_Format(PyExc_ValueError, "mb02ud: parameter %d is invalid", -info);
        return NULL;
    }

    npy_intp q_dims[2] = {l, l};
    npy_intp q_strides[2] = {sizeof(f64), l * sizeof(f64)};
    PyObject* q_out;
    if (nfct) {
        q_out = PyArray_New(&PyArray_Type, 2, q_dims, NPY_DOUBLE,
                            q_strides, q_data, 0, NPY_ARRAY_FARRAY, NULL);
        if (q_out == NULL) {
            free(q_data);
            free(sv_data);
            if (pinv) free(rp_data);
            Py_DECREF(r_array);
            Py_DECREF(b_array);
            return NULL;
        }
        PyArray_ENABLEFLAGS((PyArrayObject*)q_out, NPY_ARRAY_OWNDATA);
    } else {
        q_out = (PyObject*)q_array;
    }

    npy_intp sv_dims[1] = {l};
    PyObject* sv_out;
    if (nfct) {
        sv_out = PyArray_SimpleNewFromData(1, sv_dims, NPY_DOUBLE, sv_data);
        if (sv_out == NULL) {
            free(sv_data);
            Py_DECREF(q_out);
            if (pinv) free(rp_data);
            Py_DECREF(r_array);
            Py_DECREF(b_array);
            return NULL;
        }
        PyArray_ENABLEFLAGS((PyArrayObject*)sv_out, NPY_ARRAY_OWNDATA);
    } else {
        sv_out = (PyObject*)sv_array;
    }

    PyObject* rp_out;
    if (pinv && rank > 0) {
        npy_intp rp_dims[2] = {l, l};
        npy_intp rp_strides[2] = {sizeof(f64), l * sizeof(f64)};
        rp_out = PyArray_New(&PyArray_Type, 2, rp_dims, NPY_DOUBLE,
                             rp_strides, rp_data, 0, NPY_ARRAY_FARRAY, NULL);
        if (rp_out == NULL) {
            free(rp_data);
            Py_DECREF(q_out);
            Py_DECREF(sv_out);
            Py_DECREF(r_array);
            Py_DECREF(b_array);
            return NULL;
        }
        PyArray_ENABLEFLAGS((PyArrayObject*)rp_out, NPY_ARRAY_OWNDATA);
    } else {
        if (pinv) free(rp_data);
        Py_INCREF(Py_None);
        rp_out = Py_None;
    }

    PyObject* result = Py_BuildValue("OOOiOi", b_array, q_out, sv_out, rank, rp_out, info);

    Py_DECREF(r_array);
    Py_DECREF(b_array);
    Py_DECREF(q_out);
    Py_DECREF(sv_out);
    Py_DECREF(rp_out);

    return result;
}

static PyObject* py_md03by(PyObject* self, PyObject* args, PyObject* kwargs) {
    static char* kwlist[] = {"cond", "n", "r", "ipvt", "diag", "qtb", "delta", "par", "rank", "tol", NULL};

    char* cond;
    i32 n, rank_in;
    f64 delta, par, tol;
    PyObject *r_obj, *ipvt_obj, *diag_obj, *qtb_obj;
    PyArrayObject *r_array, *ipvt_array, *diag_array, *qtb_array;
    i32 info;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "siOOOOddid", kwlist,
                                     &cond, &n, &r_obj, &ipvt_obj, &diag_obj, &qtb_obj,
                                     &delta, &par, &rank_in, &tol)) {
        return NULL;
    }

    if (n < 0) {
        PyErr_SetString(PyExc_ValueError, "n must be non-negative");
        return NULL;
    }

    if (delta <= 0.0) {
        PyErr_SetString(PyExc_ValueError, "delta must be positive");
        return NULL;
    }

    if (par < 0.0) {
        PyErr_SetString(PyExc_ValueError, "par must be non-negative");
        return NULL;
    }

    r_array = (PyArrayObject*)PyArray_FROM_OTF(r_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (r_array == NULL) return NULL;

    ipvt_array = (PyArrayObject*)PyArray_FROM_OTF(ipvt_obj, NPY_INT32, NPY_ARRAY_IN_FARRAY);
    if (ipvt_array == NULL) {
        Py_DECREF(r_array);
        return NULL;
    }

    diag_array = (PyArrayObject*)PyArray_FROM_OTF(diag_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
    if (diag_array == NULL) {
        Py_DECREF(r_array);
        Py_DECREF(ipvt_array);
        return NULL;
    }

    qtb_array = (PyArrayObject*)PyArray_FROM_OTF(qtb_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
    if (qtb_array == NULL) {
        Py_DECREF(r_array);
        Py_DECREF(ipvt_array);
        Py_DECREF(diag_array);
        return NULL;
    }

    i32 ldr = (i32)PyArray_DIM(r_array, 0);
    f64* r_data = (f64*)PyArray_DATA(r_array);
    i32* ipvt_data = (i32*)PyArray_DATA(ipvt_array);
    f64* diag_data = (f64*)PyArray_DATA(diag_array);
    f64* qtb_data = (f64*)PyArray_DATA(qtb_array);

    bool econd = (*cond == 'E' || *cond == 'e');
    i32 ldwork = econd ? 4*n : 2*n;
    if (ldwork < 1) ldwork = 1;
    f64* dwork = (f64*)malloc(ldwork * sizeof(f64));
    if (dwork == NULL) {
        Py_DECREF(r_array);
        Py_DECREF(ipvt_array);
        Py_DECREF(diag_array);
        Py_DECREF(qtb_array);
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate workspace");
        return NULL;
    }

    f64* x_data = (n > 0) ? (f64*)calloc(n, sizeof(f64)) : NULL;
    if (n > 0 && x_data == NULL) {
        free(dwork);
        Py_DECREF(r_array);
        Py_DECREF(ipvt_array);
        Py_DECREF(diag_array);
        Py_DECREF(qtb_array);
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate x array");
        return NULL;
    }

    f64* rx_data = (n > 0) ? (f64*)calloc(n, sizeof(f64)) : NULL;
    if (n > 0 && rx_data == NULL) {
        free(x_data);
        free(dwork);
        Py_DECREF(r_array);
        Py_DECREF(ipvt_array);
        Py_DECREF(diag_array);
        Py_DECREF(qtb_array);
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate rx array");
        return NULL;
    }

    i32 rank = rank_in;

    md03by(cond, n, r_data, ldr, ipvt_data, diag_data, qtb_data, delta,
           &par, &rank, x_data, rx_data, tol, dwork, ldwork, &info);

    free(dwork);

    PyArray_ResolveWritebackIfCopy(r_array);

    npy_intp x_dims[1] = {n > 0 ? n : 0};
    PyObject* x_array = (n > 0) ? PyArray_SimpleNewFromData(1, x_dims, NPY_DOUBLE, x_data)
                                 : PyArray_EMPTY(1, x_dims, NPY_DOUBLE, 0);
    if (x_array == NULL) {
        free(x_data);
        free(rx_data);
        Py_DECREF(r_array);
        Py_DECREF(ipvt_array);
        Py_DECREF(diag_array);
        Py_DECREF(qtb_array);
        PyErr_SetString(PyExc_MemoryError, "Failed to create x array");
        return NULL;
    }
    if (n > 0) PyArray_ENABLEFLAGS((PyArrayObject*)x_array, NPY_ARRAY_OWNDATA);

    PyObject* rx_array = (n > 0) ? PyArray_SimpleNewFromData(1, x_dims, NPY_DOUBLE, rx_data)
                                  : PyArray_EMPTY(1, x_dims, NPY_DOUBLE, 0);
    if (rx_array == NULL) {
        free(rx_data);
        Py_DECREF(r_array);
        Py_DECREF(ipvt_array);
        Py_DECREF(diag_array);
        Py_DECREF(qtb_array);
        Py_DECREF(x_array);
        PyErr_SetString(PyExc_MemoryError, "Failed to create rx array");
        return NULL;
    }
    if (n > 0) PyArray_ENABLEFLAGS((PyArrayObject*)rx_array, NPY_ARRAY_OWNDATA);

    if (info < 0) {
        Py_DECREF(r_array);
        Py_DECREF(ipvt_array);
        Py_DECREF(diag_array);
        Py_DECREF(qtb_array);
        Py_DECREF(x_array);
        Py_DECREF(rx_array);
        PyErr_Format(PyExc_ValueError, "md03by: parameter %d is invalid", -info);
        return NULL;
    }

    PyObject* result = Py_BuildValue("OdiOOi", r_array, par, rank, x_array, rx_array, info);

    Py_DECREF(r_array);
    Py_DECREF(ipvt_array);
    Py_DECREF(diag_array);
    Py_DECREF(qtb_array);
    Py_DECREF(x_array);
    Py_DECREF(rx_array);

    return result;
}

/* Python wrapper for md03bb */
static PyObject* py_md03bb(PyObject* self, PyObject* args, PyObject* kwargs) {
    static char* kwlist[] = {"cond", "n", "ipar", "r", "ipvt", "diag", "qtb", "delta", "par", "ranks", "tol", NULL};

    char* cond;
    i32 n;
    f64 delta, par, tol;
    PyObject *ipar_obj, *r_obj, *ipvt_obj, *diag_obj, *qtb_obj, *ranks_obj;
    PyArrayObject *ipar_array, *r_array, *ipvt_array, *diag_array, *qtb_array, *ranks_array;
    i32 info;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "siOOOOOddOd", kwlist,
                                     &cond, &n, &ipar_obj, &r_obj, &ipvt_obj, &diag_obj, &qtb_obj,
                                     &delta, &par, &ranks_obj, &tol)) {
        return NULL;
    }

    if (n < 0) {
        PyErr_SetString(PyExc_ValueError, "n must be non-negative");
        return NULL;
    }

    if (delta <= 0.0) {
        PyErr_SetString(PyExc_ValueError, "delta must be positive");
        return NULL;
    }

    if (par < 0.0) {
        PyErr_SetString(PyExc_ValueError, "par must be non-negative");
        return NULL;
    }
    
    ipar_array = (PyArrayObject*)PyArray_FROM_OTF(ipar_obj, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    if (ipar_array == NULL) return NULL;

    r_array = (PyArrayObject*)PyArray_FROM_OTF(r_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (r_array == NULL) {
        Py_DECREF(ipar_array);
        return NULL;
    }

    ipvt_array = (PyArrayObject*)PyArray_FROM_OTF(ipvt_obj, NPY_INT32, NPY_ARRAY_IN_FARRAY);
    if (ipvt_array == NULL) {
        Py_DECREF(ipar_array);
        Py_DECREF(r_array);
        return NULL;
    }

    diag_array = (PyArrayObject*)PyArray_FROM_OTF(diag_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
    if (diag_array == NULL) {
        Py_DECREF(ipar_array);
        Py_DECREF(r_array);
        Py_DECREF(ipvt_array);
        return NULL;
    }

    qtb_array = (PyArrayObject*)PyArray_FROM_OTF(qtb_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
    if (qtb_array == NULL) {
        Py_DECREF(ipar_array);
        Py_DECREF(r_array);
        Py_DECREF(ipvt_array);
        Py_DECREF(diag_array);
        return NULL;
    }

    ranks_array = (PyArrayObject*)PyArray_FROM_OTF(ranks_obj, NPY_INT32, NPY_ARRAY_IN_ARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (ranks_array == NULL) {
        Py_DECREF(ipar_array);
        Py_DECREF(r_array);
        Py_DECREF(ipvt_array);
        Py_DECREF(diag_array);
        Py_DECREF(qtb_array);
        return NULL;
    }

    i32 lipar = (i32)PyArray_SIZE(ipar_array);
    i32 *ipar_data = (i32*)PyArray_DATA(ipar_array);
    
    i32 ldr = (i32)PyArray_DIM(r_array, 0);
    f64* r_data = (f64*)PyArray_DATA(r_array);
    i32* ipvt_data = (i32*)PyArray_DATA(ipvt_array);
    f64* diag_data = (f64*)PyArray_DATA(diag_array);
    f64* qtb_data = (f64*)PyArray_DATA(qtb_array);
    i32* ranks_data = (i32*)PyArray_DATA(ranks_array);

    bool econd = (*cond == 'E' || *cond == 'e');
    i32 ldwork = econd ? 4*n : 2*n;
    if (ldwork < 1) ldwork = 1;
    f64* dwork = (f64*)malloc(ldwork * sizeof(f64));
    if (dwork == NULL) {
        Py_DECREF(ipar_array);
        Py_DECREF(r_array);
        Py_DECREF(ipvt_array);
        Py_DECREF(diag_array);
        Py_DECREF(qtb_array);
        Py_DECREF(ranks_array);
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate workspace");
        return NULL;
    }

    f64* x_data = (n > 0) ? (f64*)calloc(n, sizeof(f64)) : NULL;
    if (n > 0 && x_data == NULL) {
        free(dwork);
        Py_DECREF(ipar_array);
        Py_DECREF(r_array);
        Py_DECREF(ipvt_array);
        Py_DECREF(diag_array);
        Py_DECREF(qtb_array);
        Py_DECREF(ranks_array);
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate x array");
        return NULL;
    }

    f64* rx_data = (n > 0) ? (f64*)calloc(n, sizeof(f64)) : NULL;
    if (n > 0 && rx_data == NULL) {
        free(x_data);
        free(dwork);
        Py_DECREF(ipar_array);
        Py_DECREF(r_array);
        Py_DECREF(ipvt_array);
        Py_DECREF(diag_array);
        Py_DECREF(qtb_array);
        Py_DECREF(ranks_array);
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate rx array");
        return NULL;
    }

    md03bb(cond, n, ipar_data, lipar, r_data, ldr, ipvt_data, diag_data, qtb_data, delta,
           &par, ranks_data, x_data, rx_data, tol, dwork, ldwork, &info);

    free(dwork);

    PyArray_ResolveWritebackIfCopy(r_array);
    PyArray_ResolveWritebackIfCopy(ranks_array);
    Py_DECREF(ipar_array);

    npy_intp x_dims[1] = {n > 0 ? n : 0};
    PyObject* x_array = (n > 0) ? PyArray_SimpleNewFromData(1, x_dims, NPY_DOUBLE, x_data)
                                 : PyArray_EMPTY(1, x_dims, NPY_DOUBLE, 0);
    if (x_array == NULL) {
        free(x_data);
        free(rx_data);
        Py_DECREF(r_array);
        Py_DECREF(ipvt_array);
        Py_DECREF(diag_array);
        Py_DECREF(qtb_array);
        Py_DECREF(ranks_array);
        PyErr_SetString(PyExc_MemoryError, "Failed to create x array");
        return NULL;
    }
    if (n > 0) PyArray_ENABLEFLAGS((PyArrayObject*)x_array, NPY_ARRAY_OWNDATA);

    PyObject* rx_array = (n > 0) ? PyArray_SimpleNewFromData(1, x_dims, NPY_DOUBLE, rx_data)
                                  : PyArray_EMPTY(1, x_dims, NPY_DOUBLE, 0);
    if (rx_array == NULL) {
        free(rx_data);
        Py_DECREF(r_array);
        Py_DECREF(ipvt_array);
        Py_DECREF(diag_array);
        Py_DECREF(qtb_array);
        Py_DECREF(x_array);
        Py_DECREF(ranks_array);
        PyErr_SetString(PyExc_MemoryError, "Failed to create rx array");
        return NULL;
    }
    if (n > 0) PyArray_ENABLEFLAGS((PyArrayObject*)rx_array, NPY_ARRAY_OWNDATA);

    if (info < 0) {
        Py_DECREF(r_array);
        Py_DECREF(ipvt_array);
        Py_DECREF(diag_array);
        Py_DECREF(qtb_array);
        Py_DECREF(x_array);
        Py_DECREF(rx_array);
        Py_DECREF(ranks_array);
        PyErr_Format(PyExc_ValueError, "md03bb: parameter %d is invalid", -info);
        return NULL;
    }

    PyObject* result = Py_BuildValue("OdOOOi", r_array, par, ranks_array, x_array, rx_array, info);

    Py_DECREF(r_array);
    Py_DECREF(ipvt_array);
    Py_DECREF(diag_array);
    Py_DECREF(qtb_array);
    Py_DECREF(x_array);
    Py_DECREF(rx_array);
    Py_DECREF(ranks_array);

    return result;
}
static PyObject* py_sg03br(PyObject* self, PyObject* args) {
    f64 xr, xi, yr, yi;
    f64 c, sr, si, zr, zi;

    if (!PyArg_ParseTuple(args, "dddd", &xr, &xi, &yr, &yi)) {
        return NULL;
    }

    sg03br(xr, xi, yr, yi, &c, &sr, &si, &zr, &zi);

    return Py_BuildValue("(ddddd)", c, sr, si, zr, zi);
}

/* Python wrapper for sg03bw */
static PyObject* py_sg03bw(PyObject* self, PyObject* args) {
    char* trans;
    PyObject *a_obj, *c_obj, *e_obj, *d_obj, *x_obj;
    PyArrayObject *a_array, *c_array, *e_array, *d_array, *x_array;
    f64 scale;
    i32 info;

    if (!PyArg_ParseTuple(args, "sOOOOO", &trans, &a_obj, &e_obj, &c_obj, &d_obj, &x_obj)) {
        return NULL;
    }

    a_array = (PyArrayObject*)PyArray_FROM_OTF(a_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
    if (a_array == NULL) return NULL;

    e_array = (PyArrayObject*)PyArray_FROM_OTF(e_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
    if (e_array == NULL) {
        Py_DECREF(a_array);
        return NULL;
    }

    c_array = (PyArrayObject*)PyArray_FROM_OTF(c_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
    if (c_array == NULL) {
        Py_DECREF(a_array);
        Py_DECREF(e_array);
        return NULL;
    }

    d_array = (PyArrayObject*)PyArray_FROM_OTF(d_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
    if (d_array == NULL) {
        Py_DECREF(a_array);
        Py_DECREF(e_array);
        Py_DECREF(c_array);
        return NULL;
    }

    x_array = (PyArrayObject*)PyArray_FROM_OTF(x_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (x_array == NULL) {
        Py_DECREF(a_array);
        Py_DECREF(e_array);
        Py_DECREF(c_array);
        Py_DECREF(d_array);
        return NULL;
    }

    npy_intp *a_dims = PyArray_DIMS(a_array);
    npy_intp *c_dims = PyArray_DIMS(c_array);
    npy_intp *x_dims = PyArray_DIMS(x_array);

    i32 m = (i32)a_dims[0];
    i32 n = (i32)c_dims[0];
    i32 lda = (i32)a_dims[0];
    i32 ldc = (i32)c_dims[0];
    i32 lde = (i32)PyArray_DIM(e_array, 0);
    i32 ldd = (i32)PyArray_DIM(d_array, 0);
    i32 ldx = (i32)x_dims[0];

    f64 *a_data = (f64*)PyArray_DATA(a_array);
    f64 *e_data = (f64*)PyArray_DATA(e_array);
    f64 *c_data = (f64*)PyArray_DATA(c_array);
    f64 *d_data = (f64*)PyArray_DATA(d_array);
    f64 *x_data = (f64*)PyArray_DATA(x_array);

    sg03bw(trans, m, n, a_data, lda, c_data, ldc, e_data, lde, d_data, ldd, x_data, ldx, &scale, &info);

    PyObject *result = Py_BuildValue("Odi", x_array, scale, info);

    Py_DECREF(a_array);
    Py_DECREF(e_array);
    Py_DECREF(c_array);
    Py_DECREF(d_array);
    Py_DECREF(x_array);

    return result;
}

/* Python wrapper for sg03bu */
static PyObject* py_sg03bu(PyObject* self, PyObject* args) {
    char *trans;
    PyObject *a_obj, *e_obj, *b_obj;
    PyArrayObject *a_array, *e_array, *b_array;
    f64 scale;
    i32 info;

    if (!PyArg_ParseTuple(args, "sOOO", &trans, &a_obj, &e_obj, &b_obj)) {
        return NULL;
    }

    a_array = (PyArrayObject*)PyArray_FROM_OTF(a_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
    if (a_array == NULL) return NULL;

    e_array = (PyArrayObject*)PyArray_FROM_OTF(e_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
    if (e_array == NULL) {
        Py_DECREF(a_array);
        return NULL;
    }

    b_array = (PyArrayObject*)PyArray_FROM_OTF(b_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (b_array == NULL) {
        Py_DECREF(a_array);
        Py_DECREF(e_array);
        return NULL;
    }

    npy_intp n_rows = PyArray_DIM(b_array, 0);
    npy_intp n_cols = PyArray_DIM(b_array, 1);
    i32 n = (i32)(n_rows < n_cols ? n_rows : n_cols);
    i32 lda = (i32)PyArray_DIM(a_array, 0);
    i32 lde = (i32)PyArray_DIM(e_array, 0);
    i32 ldb = (i32)PyArray_DIM(b_array, 0);

    f64 *a_data = (f64*)PyArray_DATA(a_array);
    f64 *e_data = (f64*)PyArray_DATA(e_array);
    f64 *b_data = (f64*)PyArray_DATA(b_array);

    i32 ldwork = (n > 1) ? 6 * n - 6 : 1;
    f64 *dwork = (f64*)malloc(ldwork * sizeof(f64));

    if (dwork == NULL) {
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate workspace");
        Py_DECREF(a_array);
        Py_DECREF(e_array);
        Py_DECREF(b_array);
        return NULL;
    }

    sg03bu(trans, n, a_data, lda, e_data, lde, b_data, ldb, &scale, dwork, &info);

    free(dwork);

    PyObject *result = Py_BuildValue("Odi", b_array, scale, info);

    Py_DECREF(a_array);
    Py_DECREF(e_array);
    Py_DECREF(b_array);

    return result;
}

/* Python wrapper for sg03bv */
static PyObject* py_sg03bv(PyObject* self, PyObject* args) {
    char *trans;
    PyObject *a_obj, *e_obj, *b_obj;
    PyArrayObject *a_array, *e_array, *b_array;
    f64 scale;
    i32 info;

    if (!PyArg_ParseTuple(args, "sOOO", &trans, &a_obj, &e_obj, &b_obj)) {
        return NULL;
    }

    a_array = (PyArrayObject*)PyArray_FROM_OTF(a_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
    if (a_array == NULL) return NULL;

    e_array = (PyArrayObject*)PyArray_FROM_OTF(e_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
    if (e_array == NULL) {
        Py_DECREF(a_array);
        return NULL;
    }

    b_array = (PyArrayObject*)PyArray_FROM_OTF(b_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (b_array == NULL) {
        Py_DECREF(a_array);
        Py_DECREF(e_array);
        return NULL;
    }

    npy_intp n_rows = PyArray_DIM(b_array, 0);
    npy_intp n_cols = PyArray_DIM(b_array, 1);
    i32 n = (i32)(n_rows < n_cols ? n_rows : n_cols);
    i32 lda = (i32)PyArray_DIM(a_array, 0);
    i32 lde = (i32)PyArray_DIM(e_array, 0);
    i32 ldb = (i32)PyArray_DIM(b_array, 0);

    f64 *a_data = (f64*)PyArray_DATA(a_array);
    f64 *e_data = (f64*)PyArray_DATA(e_array);
    f64 *b_data = (f64*)PyArray_DATA(b_array);

    i32 ldwork = (n > 1) ? 6 * n - 6 : 1;
    f64 *dwork = (f64*)malloc(ldwork * sizeof(f64));

    if (dwork == NULL) {
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate workspace");
        Py_DECREF(a_array);
        Py_DECREF(e_array);
        Py_DECREF(b_array);
        return NULL;
    }

    sg03bv(trans, n, a_data, lda, e_data, lde, b_data, ldb, &scale, dwork, &info);

    free(dwork);

    PyObject *result = Py_BuildValue("Odi", b_array, scale, info);

    Py_DECREF(a_array);
    Py_DECREF(e_array);
    Py_DECREF(b_array);

    return result;
}

/* Python wrapper for sg03bx */
static PyObject* py_sg03bx(PyObject* self, PyObject* args) {
    char *dico, *trans;
    PyObject *a_obj, *e_obj, *b_obj;
    PyArrayObject *a_array, *e_array, *b_array;
    f64 scale;
    i32 info;

    if (!PyArg_ParseTuple(args, "ssOOO", &dico, &trans, &a_obj, &e_obj, &b_obj)) {
        return NULL;
    }

    a_array = (PyArrayObject*)PyArray_FROM_OTF(a_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
    if (a_array == NULL) return NULL;

    e_array = (PyArrayObject*)PyArray_FROM_OTF(e_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
    if (e_array == NULL) {
        Py_DECREF(a_array);
        return NULL;
    }

    b_array = (PyArrayObject*)PyArray_FROM_OTF(b_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
    if (b_array == NULL) {
        Py_DECREF(a_array);
        Py_DECREF(e_array);
        return NULL;
    }

    i32 lda = (i32)PyArray_DIM(a_array, 0);
    i32 lde = (i32)PyArray_DIM(e_array, 0);
    i32 ldb = (i32)PyArray_DIM(b_array, 0);

    f64 *a_data = (f64*)PyArray_DATA(a_array);
    f64 *e_data = (f64*)PyArray_DATA(e_array);
    f64 *b_data = (f64*)PyArray_DATA(b_array);

    npy_intp dims[2] = {2, 2};
    npy_intp strides[2] = {sizeof(f64), 2 * sizeof(f64)};

    PyObject *u_array = PyArray_New(&PyArray_Type, 2, dims, NPY_DOUBLE, strides,
                                     NULL, 0, NPY_ARRAY_FARRAY, NULL);
    PyObject *m1_array = PyArray_New(&PyArray_Type, 2, dims, NPY_DOUBLE, strides,
                                      NULL, 0, NPY_ARRAY_FARRAY, NULL);
    PyObject *m2_array = PyArray_New(&PyArray_Type, 2, dims, NPY_DOUBLE, strides,
                                      NULL, 0, NPY_ARRAY_FARRAY, NULL);

    if (u_array == NULL || m1_array == NULL || m2_array == NULL) {
        Py_XDECREF(u_array);
        Py_XDECREF(m1_array);
        Py_XDECREF(m2_array);
        Py_DECREF(a_array);
        Py_DECREF(e_array);
        Py_DECREF(b_array);
        return NULL;
    }

    f64 *u_data = (f64*)PyArray_DATA((PyArrayObject*)u_array);
    f64 *m1_data = (f64*)PyArray_DATA((PyArrayObject*)m1_array);
    f64 *m2_data = (f64*)PyArray_DATA((PyArrayObject*)m2_array);

    sg03bx(dico, trans, a_data, lda, e_data, lde, b_data, ldb,
           u_data, 2, &scale, m1_data, 2, m2_data, 2, &info);

    PyObject *result = Py_BuildValue("OdOOi", u_array, scale, m1_array, m2_array, info);

    Py_DECREF(a_array);
    Py_DECREF(e_array);
    Py_DECREF(b_array);
    Py_DECREF(u_array);
    Py_DECREF(m1_array);
    Py_DECREF(m2_array);

    return result;
}

/* Python wrapper for sg03bd */
static PyObject* py_sg03bd(PyObject* self, PyObject* args) {
    char *dico, *fact, *trans;
    i32 n, m;
    PyObject *a_obj, *e_obj, *b_obj;
    PyArrayObject *a_array, *e_array, *b_array;
    f64 scale;
    i32 info;

    if (!PyArg_ParseTuple(args, "sssiiOOO", &dico, &fact, &trans, &n, &m, &a_obj, &e_obj, &b_obj)) {
        return NULL;
    }

    a_array = (PyArrayObject*)PyArray_FROM_OTF(a_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (a_array == NULL) return NULL;

    e_array = (PyArrayObject*)PyArray_FROM_OTF(e_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (e_array == NULL) {
        Py_DECREF(a_array);
        return NULL;
    }

    /* B array is modified in place and may need more space than input provides.
     * Allocate workspace of correct size and copy input B into it. */
    b_array = (PyArrayObject*)PyArray_FROM_OTF(b_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
    if (b_array == NULL) {
        Py_DECREF(a_array);
        Py_DECREF(e_array);
        return NULL;
    }
    
    /* Validate B has 2 dimensions */
    if (PyArray_NDIM(b_array) != 2) {
        PyErr_SetString(PyExc_ValueError, "B must be a 2D array");
        Py_DECREF(a_array);
        Py_DECREF(e_array);
        Py_DECREF(b_array);
        return NULL;
    }

    i32 lda = (i32)PyArray_DIM(a_array, 0);
    i32 lde = (i32)PyArray_DIM(e_array, 0);
    
    /* Determine required LDB and allocate workspace for B */
    bool istran = (trans[0] == 'T' || trans[0] == 't');
    i32 ldb_required = istran ? (n > 1 ? n : 1) : ((m > n ? m : n) > 1 ? (m > n ? m : n) : 1);
    i32 b_cols = n;
    
    /* Get input B dimensions */
    i32 b_in_rows = (i32)PyArray_DIM(b_array, 0);
    i32 b_in_cols = (i32)PyArray_DIM(b_array, 1);
    i32 b_in_ld = (i32)PyArray_DIM(b_array, 0);  /* Leading dimension of input */
    
    /* Allocate B workspace */
    npy_intp dims_b[2] = {ldb_required, b_cols};
    npy_intp strides_b[2] = {sizeof(f64), ldb_required * sizeof(f64)};
    PyObject *b_work = PyArray_New(&PyArray_Type, 2, dims_b, NPY_DOUBLE, strides_b,
                                    NULL, 0, NPY_ARRAY_FARRAY, NULL);
    if (b_work == NULL) {
        Py_DECREF(a_array);
        Py_DECREF(e_array);
        Py_DECREF(b_array);
        return NULL;
    }
    
    /* Copy input B into workspace */
    f64 *b_work_data = (f64*)PyArray_DATA((PyArrayObject*)b_work);
    f64 *b_data = (f64*)PyArray_DATA(b_array);
    
    /* Zero out workspace first */
    for (i32 j = 0; j < b_cols; j++) {
        for (i32 i = 0; i < ldb_required; i++) {
            b_work_data[i + j*ldb_required] = 0.0;
        }
    }
    
    /* Copy input B (column by column) */
    for (i32 j = 0; j < b_in_cols && j < b_cols; j++) {
        for (i32 i = 0; i < b_in_rows && i < ldb_required; i++) {
            b_work_data[i + j*ldb_required] = b_data[i + j*b_in_ld];
        }
    }
    
    Py_DECREF(b_array);  /* Done with input B */
    i32 ldb = ldb_required;

    f64 *a_data = (f64*)PyArray_DATA(a_array);
    f64 *e_data = (f64*)PyArray_DATA(e_array);

    /* Allocate workspace for Q and Z matrices (not returned) */
    i32 ldq = (n > 1) ? n : 1;
    i32 ldz = (n > 1) ? n : 1;
    f64 *q_data = (f64*)malloc(ldq * n * sizeof(f64));
    f64 *z_data = (f64*)malloc(ldz * n * sizeof(f64));

    if (q_data == NULL || z_data == NULL) {
        free(q_data);
        free(z_data);
        Py_DECREF(a_array);
        Py_DECREF(e_array);
        Py_DECREF(b_work);
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate Q/Z workspace");
        return NULL;
    }

    npy_intp dims_eig[1] = {n};
    PyObject *alphar_array = PyArray_New(&PyArray_Type, 1, dims_eig, NPY_DOUBLE, NULL,
                                          NULL, 0, NPY_ARRAY_FARRAY, NULL);
    PyObject *alphai_array = PyArray_New(&PyArray_Type, 1, dims_eig, NPY_DOUBLE, NULL,
                                          NULL, 0, NPY_ARRAY_FARRAY, NULL);
    PyObject *beta_array = PyArray_New(&PyArray_Type, 1, dims_eig, NPY_DOUBLE, NULL,
                                        NULL, 0, NPY_ARRAY_FARRAY, NULL);

    if (alphar_array == NULL || alphai_array == NULL || beta_array == NULL) {
        Py_XDECREF(alphar_array);
        Py_XDECREF(alphai_array);
        Py_XDECREF(beta_array);
        free(q_data);
        free(z_data);
        Py_DECREF(a_array);
        Py_DECREF(e_array);
        Py_DECREF(b_work);
        return NULL;
    }

    f64 *alphar_data = (f64*)PyArray_DATA((PyArrayObject*)alphar_array);
    f64 *alphai_data = (f64*)PyArray_DATA((PyArrayObject*)alphai_array);
    f64 *beta_data = (f64*)PyArray_DATA((PyArrayObject*)beta_array);

    i32 ldwork;
    if (fact[0] == 'F' || fact[0] == 'f') {
        ldwork = (2 * n > 6 * n - 6) ? 2 * n : (6 * n - 6);
        ldwork = (ldwork > 1) ? ldwork : 1;
    } else {
        ldwork = (4 * n > 6 * n - 6) ? 4 * n : (6 * n - 6);
        ldwork = (ldwork > 1) ? ldwork : 1;
        i32 mingg = (ldwork > 8 * n + 16) ? ldwork : (8 * n + 16);
        ldwork = mingg;
    }

    f64 *dwork = (f64*)malloc(ldwork * sizeof(f64));
    if (dwork == NULL) {
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate workspace");
        Py_DECREF(alphar_array);
        Py_DECREF(alphai_array);
        Py_DECREF(beta_array);
        free(q_data);
        free(z_data);
        Py_DECREF(a_array);
        Py_DECREF(e_array);
        Py_DECREF(b_work);
        return NULL;
    }

    sg03bd(dico, fact, trans, n, m, a_data, lda, e_data, lde,
           q_data, ldq, z_data, ldz, b_work_data, ldb, &scale,
           alphar_data, alphai_data, beta_data, dwork, ldwork, &info);

    free(dwork);
    free(q_data);
    free(z_data);

    /* Resolve writebackifcopy before decref */
    PyArray_ResolveWritebackIfCopy(a_array);
    PyArray_ResolveWritebackIfCopy(e_array);

    /* sg03bd modifies B in place to produce U (n x n upper triangular).
     * Create a view of the n x n submatrix from b_work.
     * See CLAUDE.md: "CRITICAL: In-place modification - return input array directly"
     */
    npy_intp dims_u[2] = {n, n};
    npy_intp strides_u[2] = {sizeof(f64), ldb * sizeof(f64)};
    PyObject *u_array = PyArray_New(&PyArray_Type, 2, dims_u, NPY_DOUBLE, strides_u,
                                     b_work_data, 0, NPY_ARRAY_FARRAY, NULL);
    if (u_array == NULL) {
        Py_DECREF(a_array);
        Py_DECREF(e_array);
        Py_DECREF(b_work);
        Py_DECREF(alphar_array);
        Py_DECREF(alphai_array);
        Py_DECREF(beta_array);
        return NULL;
    }
    
    /* Set b_work as the base object so the memory stays alive */
    if (PyArray_SetBaseObject((PyArrayObject*)u_array, (PyObject*)b_work) < 0) {
        Py_DECREF(u_array);
        Py_DECREF(a_array);
        Py_DECREF(e_array);
        Py_DECREF(b_work);
        Py_DECREF(alphar_array);
        Py_DECREF(alphai_array);
        Py_DECREF(beta_array);
        return NULL;
    }

    PyObject *result = Py_BuildValue("OdOOOi", u_array, scale,
                                     alphar_array, alphai_array, beta_array, info);

    /* Py_BuildValue with "O" increments refcounts, so we need to DECREF all arrays
     * that were passed to it to avoid leaks */
    Py_DECREF(a_array);
    Py_DECREF(e_array);
    /* Don't DECREF b_work here - it's now owned by u_array as base object */
    Py_DECREF(u_array);
    Py_DECREF(alphar_array);
    Py_DECREF(alphai_array);
    Py_DECREF(beta_array);

    return result;
}

/* Python wrapper for tb01vy */
static PyObject* py_tb01vy(PyObject* self, PyObject* args, PyObject* kwargs) {
    i32 n, m, l, ltheta;
    const char *apply = "N";
    PyObject *theta_obj;
    PyArrayObject *theta_array;

    static char *kwlist[] = {"n", "m", "l", "theta", "apply", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "iiiO|s", kwlist,
                                     &n, &m, &l, &theta_obj, &apply)) {
        return NULL;
    }

    if (n < 0 || m < 0 || l < 0) {
        PyErr_Format(PyExc_ValueError, "Dimensions must be non-negative (n=%d, m=%d, l=%d)", n, m, l);
        return NULL;
    }

    theta_array = (PyArrayObject*)PyArray_FROM_OTF(theta_obj, NPY_DOUBLE,
                                                   NPY_ARRAY_IN_FARRAY);
    if (theta_array == NULL) return NULL;

    ltheta = (i32)PyArray_SIZE(theta_array);
    i32 required_ltheta = n * (l + m + 1) + l * m;

    if (ltheta < required_ltheta) {
        Py_DECREF(theta_array);
        PyErr_Format(PyExc_ValueError, "ltheta=%d is too small (need >= %d)", ltheta, required_ltheta);
        return NULL;
    }

    if (!(apply[0] == 'A' || apply[0] == 'a' || apply[0] == 'N' || apply[0] == 'n')) {
        Py_DECREF(theta_array);
        PyErr_Format(PyExc_ValueError, "apply must be 'A' or 'N', got '%s'", apply);
        return NULL;
    }

    i32 lda = n > 0 ? n : 1;
    i32 ldb = n > 0 ? n : 1;
    i32 ldc = l > 0 ? l : 1;
    i32 ldd = l > 0 ? l : 1;
    i32 ldwork = n * (n + l + 1);

    npy_intp a_dims[2] = {n, n};
    npy_intp a_strides[2] = {sizeof(f64), lda * sizeof(f64)};
    f64 *a_data = (f64*)calloc(lda * n, sizeof(f64));

    npy_intp b_dims[2] = {n, m};
    npy_intp b_strides[2] = {sizeof(f64), ldb * sizeof(f64)};
    f64 *b_data = (f64*)calloc(ldb * m, sizeof(f64));

    npy_intp c_dims[2] = {l, n};
    npy_intp c_strides[2] = {sizeof(f64), ldc * sizeof(f64)};
    f64 *c_data = (f64*)calloc(ldc * n, sizeof(f64));

    npy_intp d_dims[2] = {l, m};
    npy_intp d_strides[2] = {sizeof(f64), ldd * sizeof(f64)};
    f64 *d_data = (f64*)calloc(ldd * m, sizeof(f64));

    f64 *x0 = (f64*)calloc(n > 0 ? n : 1, sizeof(f64));
    f64 *dwork = (f64*)malloc(ldwork > 0 ? ldwork * sizeof(f64) : sizeof(f64));

    if (a_data == NULL || b_data == NULL || c_data == NULL ||
        d_data == NULL || x0 == NULL || dwork == NULL) {
        free(a_data); free(b_data); free(c_data); free(d_data); free(x0); free(dwork);
        Py_DECREF(theta_array);
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate memory");
        return NULL;
    }

    i32 info;
    f64 *theta_data = (f64*)PyArray_DATA(theta_array);

    tb01vy(apply, n, m, l, theta_data, ltheta, a_data, lda, b_data, ldb,
           c_data, ldc, d_data, ldd, x0, dwork, ldwork, &info);

    free(dwork);

    PyObject *a_array = PyArray_New(&PyArray_Type, 2, a_dims, NPY_DOUBLE,
                                    a_strides, a_data, 0, NPY_ARRAY_FARRAY, NULL);
    PyArray_ENABLEFLAGS((PyArrayObject*)a_array, NPY_ARRAY_OWNDATA);

    PyObject *b_array = PyArray_New(&PyArray_Type, 2, b_dims, NPY_DOUBLE,
                                    b_strides, b_data, 0, NPY_ARRAY_FARRAY, NULL);
    PyArray_ENABLEFLAGS((PyArrayObject*)b_array, NPY_ARRAY_OWNDATA);

    PyObject *c_array = PyArray_New(&PyArray_Type, 2, c_dims, NPY_DOUBLE,
                                    c_strides, c_data, 0, NPY_ARRAY_FARRAY, NULL);
    PyArray_ENABLEFLAGS((PyArrayObject*)c_array, NPY_ARRAY_OWNDATA);

    PyObject *d_array = PyArray_New(&PyArray_Type, 2, d_dims, NPY_DOUBLE,
                                    d_strides, d_data, 0, NPY_ARRAY_FARRAY, NULL);
    PyArray_ENABLEFLAGS((PyArrayObject*)d_array, NPY_ARRAY_OWNDATA);

    npy_intp x0_dims[1] = {n};
    PyObject *x0_array = PyArray_SimpleNewFromData(1, x0_dims, NPY_DOUBLE, x0);
    PyArray_ENABLEFLAGS((PyArrayObject*)x0_array, NPY_ARRAY_OWNDATA);

    PyObject *result = Py_BuildValue("OOOOOi", a_array, b_array, c_array, d_array, x0_array, info);

    Py_DECREF(theta_array);
    Py_DECREF(a_array);
    Py_DECREF(b_array);
    Py_DECREF(c_array);
    Py_DECREF(d_array);
    Py_DECREF(x0_array);

    return result;
}

/* Python wrapper for ma02ad */
static PyObject* py_tb01wd(PyObject* self, PyObject* args) {
    i32 n, m, p;
    PyObject *a_obj, *b_obj, *c_obj;
    PyArrayObject *a_array, *b_array, *c_array;

    if (!PyArg_ParseTuple(args, "iiiOOO", &n, &m, &p, &a_obj, &b_obj, &c_obj)) {
        return NULL;
    }

    /* Validate dimensions */
    if (n < 0 || m < 0 || p < 0) {
        PyErr_Format(PyExc_ValueError, "Dimensions must be non-negative (n=%d, m=%d, p=%d)", n, m, p);
        return NULL;
    }

    /* Convert inputs to NumPy arrays - preserve Fortran-order */
    a_array = (PyArrayObject*)PyArray_FROM_OTF(a_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (a_array == NULL) return NULL;

    b_array = (PyArrayObject*)PyArray_FROM_OTF(b_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (b_array == NULL) {
        Py_DECREF(a_array);
        return NULL;
    }

    c_array = (PyArrayObject*)PyArray_FROM_OTF(c_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (c_array == NULL) {
        Py_DECREF(a_array);
        Py_DECREF(b_array);
        return NULL;
    }

    /* Get leading dimensions */
    i32 lda = n > 0 ? n : 1;
    i32 ldb = n > 0 ? n : 1;
    i32 ldc = p > 0 ? p : 1;
    i32 ldu = n > 0 ? n : 1;

    /* Allocate output arrays */
    npy_intp u_dims[2] = {n, n};
    npy_intp u_strides[2] = {sizeof(f64), ldu * sizeof(f64)};
    f64 *u_data = (f64*)calloc(ldu * n, sizeof(f64));

    f64 *wr = (f64*)calloc(n > 0 ? n : 1, sizeof(f64));
    f64 *wi = (f64*)calloc(n > 0 ? n : 1, sizeof(f64));

    i32 ldwork = n > 0 ? 3 * n : 1;
    f64 *dwork = (f64*)malloc(ldwork * sizeof(f64));

    if (u_data == NULL || wr == NULL || wi == NULL || dwork == NULL) {
        free(u_data);
        free(wr);
        free(wi);
        free(dwork);
        Py_DECREF(a_array);
        Py_DECREF(b_array);
        Py_DECREF(c_array);
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate workspace");
        return NULL;
    }

    /* Call C function */
    i32 info;
    f64 *a_data = (f64*)PyArray_DATA(a_array);
    f64 *b_data = (f64*)PyArray_DATA(b_array);
    f64 *c_data = (f64*)PyArray_DATA(c_array);

    tb01wd(n, m, p, a_data, lda, b_data, ldb, c_data, ldc,
           u_data, ldu, wr, wi, dwork, ldwork, &info);

    free(dwork);

    /* Create output arrays */
    PyObject *u_array = PyArray_New(&PyArray_Type, 2, u_dims, NPY_DOUBLE,
                                    u_strides, u_data, 0, NPY_ARRAY_FARRAY, NULL);
    PyArray_ENABLEFLAGS((PyArrayObject*)u_array, NPY_ARRAY_OWNDATA);

    npy_intp wr_dims[1] = {n};
    PyObject *wr_array = PyArray_SimpleNewFromData(1, wr_dims, NPY_DOUBLE, wr);
    PyArray_ENABLEFLAGS((PyArrayObject*)wr_array, NPY_ARRAY_OWNDATA);

    npy_intp wi_dims[1] = {n};
    PyObject *wi_array = PyArray_SimpleNewFromData(1, wi_dims, NPY_DOUBLE, wi);
    PyArray_ENABLEFLAGS((PyArrayObject*)wi_array, NPY_ARRAY_OWNDATA);

    /* Build result tuple */
    PyObject *result = Py_BuildValue("OOOOOOi", a_array, b_array, c_array,
                                     u_array, wr_array, wi_array, info);

    Py_DECREF(a_array);
    Py_DECREF(b_array);
    Py_DECREF(c_array);
    Py_DECREF(u_array);
    Py_DECREF(wr_array);
    Py_DECREF(wi_array);

    return result;
}

static PyObject* py_ma02ed(PyObject* self, PyObject* args) {
    const char *uplo_str;
    PyObject *a_obj;
    PyArrayObject *a_array;

    if (!PyArg_ParseTuple(args, "sO", &uplo_str, &a_obj)) {
        return NULL;
    }

    if (uplo_str == NULL || uplo_str[0] == '\0') {
        PyErr_SetString(PyExc_ValueError, "uplo must be a non-empty string");
        return NULL;
    }
    char uplo = uplo_str[0];

    /* Convert to NumPy array - preserve Fortran-order, in-place modification */
    a_array = (PyArrayObject*)PyArray_FROM_OTF(a_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (a_array == NULL) {
        return NULL;
    }

    /* Get dimensions */
    npy_intp *a_dims = PyArray_DIMS(a_array);
    i32 n = (i32)a_dims[0];
    i32 lda = n > 0 ? n : 1;

    /* Call C function - modifies a in place */
    f64 *a_data = (f64*)PyArray_DATA(a_array);
    ma02ed(uplo, n, a_data, lda);

    /* Resolve writebackifcopy before returning */
    PyArray_ResolveWritebackIfCopy(a_array);

    /* Return modified array */
    PyObject *result = Py_BuildValue("O", a_array);
    Py_DECREF(a_array);
    return result;
}

static PyObject* py_mb04kd(PyObject* self, PyObject* args) {
    const char *uplo_str;
    i32 n, m, p;
    PyObject *r_obj, *a_obj, *b_obj;
    PyArrayObject *r_array, *a_array, *b_array;

    if (!PyArg_ParseTuple(args, "siiiOOO", &uplo_str, &n, &m, &p,
                          &r_obj, &a_obj, &b_obj)) {
        return NULL;
    }

    if (uplo_str == NULL || uplo_str[0] == '\0') {
        PyErr_SetString(PyExc_ValueError, "uplo must be a non-empty string");
        return NULL;
    }
    char uplo = uplo_str[0];

    /* Convert input arrays - in-place modification */
    r_array = (PyArrayObject*)PyArray_FROM_OTF(r_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (r_array == NULL) return NULL;

    a_array = (PyArrayObject*)PyArray_FROM_OTF(a_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (a_array == NULL) {
        Py_DECREF(r_array);
        return NULL;
    }

    b_array = (PyArrayObject*)PyArray_FROM_OTF(b_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (b_array == NULL) {
        Py_DECREF(r_array);
        Py_DECREF(a_array);
        return NULL;
    }

    /* Get dimensions */
    npy_intp *r_dims = PyArray_DIMS(r_array);
    npy_intp *a_dims = PyArray_DIMS(a_array);
    npy_intp *b_dims = PyArray_DIMS(b_array);

    i32 ldr = (i32)r_dims[0];
    i32 lda = (i32)a_dims[0];
    i32 ldb = (i32)b_dims[0];
    i32 ldc = n > 0 ? n : 1;

    /* Allocate output arrays */
    npy_intp c_dims[2] = {n, m};
    npy_intp c_strides[2] = {sizeof(f64), ldc * sizeof(f64)};
    f64 *c_data = (f64*)calloc(ldc * m, sizeof(f64));
    if (c_data == NULL) {
        Py_DECREF(r_array);
        Py_DECREF(a_array);
        Py_DECREF(b_array);
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate C matrix");
        return NULL;
    }

    f64 *tau = (f64*)calloc(n > 0 ? n : 1, sizeof(f64));
    if (tau == NULL) {
        free(c_data);
        Py_DECREF(r_array);
        Py_DECREF(a_array);
        Py_DECREF(b_array);
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate tau");
        return NULL;
    }

    f64 *dwork = (f64*)calloc(n > 0 ? n : 1, sizeof(f64));
    if (dwork == NULL) {
        free(c_data);
        free(tau);
        Py_DECREF(r_array);
        Py_DECREF(a_array);
        Py_DECREF(b_array);
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate workspace");
        return NULL;
    }

    /* Call C function */
    f64 *r_data = (f64*)PyArray_DATA(r_array);
    f64 *a_data = (f64*)PyArray_DATA(a_array);
    f64 *b_data = (f64*)PyArray_DATA(b_array);

    mb04kd(uplo, n, m, p, r_data, ldr, a_data, lda, b_data, ldb,
           c_data, ldc, tau, dwork);

    /* Resolve writeback */
    PyArray_ResolveWritebackIfCopy(r_array);
    PyArray_ResolveWritebackIfCopy(a_array);
    PyArray_ResolveWritebackIfCopy(b_array);

    /* Create output arrays */
    PyObject *c_array = PyArray_New(&PyArray_Type, 2, c_dims, NPY_DOUBLE,
                                    c_strides, c_data, 0, NPY_ARRAY_FARRAY, NULL);
    if (c_array == NULL) {
        free(c_data);
        free(tau);
        free(dwork);
        Py_DECREF(r_array);
        Py_DECREF(a_array);
        Py_DECREF(b_array);
        return NULL;
    }
    PyArray_ENABLEFLAGS((PyArrayObject*)c_array, NPY_ARRAY_OWNDATA);

    npy_intp tau_dims[1] = {n > 0 ? n : 1};
    PyObject *tau_array = PyArray_SimpleNewFromData(1, tau_dims, NPY_DOUBLE, tau);
    if (tau_array == NULL) {
        free(tau);
        free(dwork);
        Py_DECREF(c_array);
        Py_DECREF(r_array);
        Py_DECREF(a_array);
        Py_DECREF(b_array);
        return NULL;
    }
    PyArray_ENABLEFLAGS((PyArrayObject*)tau_array, NPY_ARRAY_OWNDATA);

    free(dwork);

    /* Return (R_bar, A_out, D, C, tau) */
    PyObject *result = Py_BuildValue("OOOOO", r_array, a_array, b_array, c_array, tau_array);
    Py_DECREF(r_array);
    Py_DECREF(a_array);
    Py_DECREF(b_array);
    Py_DECREF(c_array);
    Py_DECREF(tau_array);

    return result;
}

static PyObject* py_ma02ad(PyObject* self, PyObject* args) {
    const char* job;
    PyObject *a_obj;
    PyArrayObject *a_array;

    if (!PyArg_ParseTuple(args, "sO", &job, &a_obj)) {
        return NULL;
    }

    /* Convert to NumPy array - preserve Fortran-order */
    a_array = (PyArrayObject*)PyArray_FROM_OTF(a_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY);
    if (a_array == NULL) {
        return NULL;
    }

    /* Get dimensions */
    npy_intp *a_dims = PyArray_DIMS(a_array);
    i32 m = (i32)a_dims[0];
    i32 n = (i32)a_dims[1];
    i32 lda = m > 0 ? m : 1;
    i32 ldb = n > 0 ? n : 1;

    /* Allocate output array B with shape (n, m) */
    npy_intp b_dims[2] = {n, m};
    npy_intp b_strides[2] = {sizeof(f64), ldb * sizeof(f64)};
    f64 *b_data = (f64*)calloc(ldb * m, sizeof(f64));
    if (b_data == NULL) {
        Py_DECREF(a_array);
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate output array");
        return NULL;
    }

    /* Call C function */
    f64 *a_data = (f64*)PyArray_DATA(a_array);
    ma02ad(job, m, n, a_data, lda, b_data, ldb);

    /* Create output NumPy array */
    PyObject *b_array = PyArray_New(&PyArray_Type, 2, b_dims, NPY_DOUBLE,
                                    b_strides, b_data, 0, NPY_ARRAY_FARRAY, NULL);
    PyArray_ENABLEFLAGS((PyArrayObject*)b_array, NPY_ARRAY_OWNDATA);

    Py_DECREF(a_array);
    return b_array;
}

static PyObject* py_tf01mx(PyObject* self, PyObject* args) {
    i32 n, m, p, ny;
    PyObject *s_obj, *u_obj, *x_obj;
    PyArrayObject *s_array, *u_array, *x_array;
    i32 info;

    if (!PyArg_ParseTuple(args, "iiiiOOO", &n, &m, &p, &ny, &s_obj, &u_obj, &x_obj)) {
        return NULL;
    }

    if (n < 0) {
        PyErr_SetString(PyExc_ValueError, "N must be >= 0");
        return NULL;
    }
    if (m < 0) {
        PyErr_SetString(PyExc_ValueError, "M must be >= 0");
        return NULL;
    }
    if (p < 0) {
        PyErr_SetString(PyExc_ValueError, "P must be >= 0");
        return NULL;
    }
    if (ny < 0) {
        PyErr_SetString(PyExc_ValueError, "NY must be >= 0");
        return NULL;
    }

    s_array = (PyArrayObject*)PyArray_FROM_OTF(s_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
    if (s_array == NULL) return NULL;

    u_array = (PyArrayObject*)PyArray_FROM_OTF(u_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
    if (u_array == NULL) {
        Py_DECREF(s_array);
        return NULL;
    }

    x_array = (PyArrayObject*)PyArray_FROM_OTF(x_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (x_array == NULL) {
        Py_DECREF(s_array);
        Py_DECREF(u_array);
        return NULL;
    }

    i32 lds = (i32)PyArray_DIM(s_array, 0);
    i32 ldu = (ny > 1) ? (i32)PyArray_DIM(u_array, 0) : 1;
    i32 ldy = (ny > 1) ? ny : 1;

    f64 *s_data = (f64*)PyArray_DATA(s_array);
    f64 *u_data = (f64*)PyArray_DATA(u_array);
    f64 *x_data = (f64*)PyArray_DATA(x_array);

    npy_intp dims_y[2] = {ny, p};
    npy_intp strides_y[2] = {sizeof(f64), ldy * sizeof(f64)};
    PyObject *y_array = PyArray_New(&PyArray_Type, 2, dims_y, NPY_DOUBLE, strides_y,
                                     NULL, 0, NPY_ARRAY_FARRAY, NULL);
    if (y_array == NULL) {
        Py_DECREF(s_array);
        Py_DECREF(u_array);
        Py_DECREF(x_array);
        return NULL;
    }

    f64 *y_data = (f64*)PyArray_DATA((PyArrayObject*)y_array);

    i32 ldwork;
    if (n == 0 || p == 0 || ny == 0) {
        ldwork = 0;
    } else if (m == 0) {
        ldwork = n + p;
    } else {
        ldwork = 2 * n + m + p;
    }

    f64 *dwork = NULL;
    if (ldwork > 0) {
        dwork = (f64*)malloc(ldwork * sizeof(f64));
        if (dwork == NULL) {
            PyErr_SetString(PyExc_MemoryError, "Failed to allocate workspace");
            Py_DECREF(s_array);
            Py_DECREF(u_array);
            Py_DECREF(x_array);
            Py_DECREF(y_array);
            return NULL;
        }
    }

    tf01mx(n, m, p, ny, s_data, lds, u_data, ldu, x_data, y_data, ldy, dwork, ldwork, &info);

    if (dwork != NULL) {
        free(dwork);
    }

    PyObject *result = Py_BuildValue("OOi", y_array, x_array, info);

    Py_DECREF(s_array);
    Py_DECREF(u_array);
    Py_DECREF(x_array);
    Py_DECREF(y_array);

    return result;
}

/* Python wrapper for md03ba */
static PyObject* py_md03ba(PyObject* self, PyObject* args) {
    i32 n;
    PyObject *ipar_obj, *j_obj, *e_obj;
    f64 fnorm;
    PyArrayObject *ipar_array, *j_array, *e_array;
    i32 info;

    if (!PyArg_ParseTuple(args, "iOdOO", &n, &ipar_obj, &fnorm, &j_obj, &e_obj)) {
        return NULL;
    }
    
    if (n < 0 || fnorm < 0) {
        PyErr_Format(PyExc_ValueError, "Invalid arguments");
        return NULL;
    }

    ipar_array = (PyArrayObject*)PyArray_FROM_OTF(ipar_obj, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    j_array = (PyArrayObject*)PyArray_FROM_OTF(j_obj, NPY_DOUBLE, NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    e_array = (PyArrayObject*)PyArray_FROM_OTF(e_obj, NPY_DOUBLE, NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    
    if (!ipar_array || !j_array || !e_array) {
         Py_XDECREF(ipar_array); Py_XDECREF(j_array); Py_XDECREF(e_array);
         return NULL;
    }
    
    i32 lipar = (i32)PyArray_SIZE(ipar_array);
    i32 *ipar_data = (i32*)PyArray_DATA(ipar_array);
    
    if (lipar < 1) {
         Py_DECREF(ipar_array); Py_DECREF(j_array); Py_DECREF(e_array);
         return NULL;
    }
    
    i32 ldj = (i32)PyArray_DIM(j_array, 0);
    f64 *j_data = (f64*)PyArray_DATA(j_array);
    f64 *e_data = (f64*)PyArray_DATA(e_array);
    
    npy_intp jnorms_dims[1] = {n};
    PyObject *jnorms_array = PyArray_SimpleNew(1, jnorms_dims, NPY_DOUBLE);
    
    npy_intp ipvt_dims[1] = {n};
    PyObject *ipvt_array = PyArray_SimpleNew(1, ipvt_dims, NPY_INT32);
    
    f64 gnorm;
    
    if (!jnorms_array || !ipvt_array) {
         Py_DECREF(ipar_array); Py_DECREF(j_array); Py_DECREF(e_array);
         Py_XDECREF(jnorms_array); Py_XDECREF(ipvt_array);
         return NULL;
    }
    
    f64 *jnorms_data = (f64*)PyArray_DATA((PyArrayObject*)jnorms_array);
    i32 *ipvt_data = (i32*)PyArray_DATA((PyArrayObject*)ipvt_array);
    
    i32 ldwork = (n > 1) ? (4*n + 1) : 1;
    f64 *dwork = (f64*)malloc(ldwork * sizeof(f64));
    if (!dwork) {
         Py_DECREF(ipar_array); Py_DECREF(j_array); Py_DECREF(e_array);
         Py_DECREF(jnorms_array); Py_DECREF(ipvt_array);
         return PyErr_NoMemory();
    }
    
    md03ba(n, ipar_data, lipar, fnorm, j_data, &ldj, e_data, 
           jnorms_data, &gnorm, ipvt_data, dwork, ldwork, &info);
           
    free(dwork);
    
    PyArray_ResolveWritebackIfCopy(j_array);
    PyArray_ResolveWritebackIfCopy(e_array);
    Py_DECREF(ipar_array);
    
    PyObject *result = Py_BuildValue("OOOdOi", j_array, e_array, jnorms_array, gnorm, ipvt_array, info);
    
    Py_DECREF(j_array);
    Py_DECREF(e_array);
    Py_DECREF(jnorms_array);
    Py_DECREF(ipvt_array);
    
    return result;
}

/* Python wrapper for nf01br */
static PyObject* py_nf01br(PyObject* self, PyObject* args) {
    const char *cond, *uplo, *trans;
    i32 n;
    PyObject *ipar_obj, *r_obj, *sdiag_obj, *s_obj, *b_obj, *ranks_obj;
    f64 tol;
    PyArrayObject *ipar_array, *r_array, *sdiag_array, *s_array, *b_array, *ranks_array;
    i32 info;

    if (!PyArg_ParseTuple(args, "sssiOOOOOOd", &cond, &uplo, &trans, &n, 
                          &ipar_obj, &r_obj, &sdiag_obj, &s_obj, &b_obj, &ranks_obj, &tol)) {
        return NULL;
    }
    
    ipar_array = (PyArrayObject*)PyArray_FROM_OTF(ipar_obj, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    r_array = (PyArrayObject*)PyArray_FROM_OTF(r_obj, NPY_DOUBLE, NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    sdiag_array = (PyArrayObject*)PyArray_FROM_OTF(sdiag_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    s_array = (PyArrayObject*)PyArray_FROM_OTF(s_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    b_array = (PyArrayObject*)PyArray_FROM_OTF(b_obj, NPY_DOUBLE, NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    ranks_array = (PyArrayObject*)PyArray_FROM_OTF(ranks_obj, NPY_INT32, NPY_ARRAY_IN_ARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    
    if (!ipar_array || !r_array || !sdiag_array || !s_array || !b_array || !ranks_array) {
         Py_XDECREF(ipar_array); Py_XDECREF(r_array); Py_XDECREF(sdiag_array);
         Py_XDECREF(s_array); Py_XDECREF(b_array); Py_XDECREF(ranks_array);
         return NULL;
    }
    
    i32 lipar = (i32)PyArray_SIZE(ipar_array);
    if (lipar < 4) {
         Py_DECREF(ipar_array); Py_DECREF(r_array); Py_DECREF(sdiag_array);
         Py_DECREF(s_array); Py_DECREF(b_array); Py_DECREF(ranks_array);
         PyErr_SetString(PyExc_ValueError, "ipar must have length >= 4");
         return NULL;
    }
    
    i32 *ipar_data = (i32*)PyArray_DATA(ipar_array);
    i32 st = ipar_data[0];
    i32 bn = ipar_data[1];
    i32 bsn = ipar_data[3];
    bool full = (bn <= 1 || bsn == 0);
    bool econd = (cond[0] == 'E' || cond[0] == 'e');
    
    i32 lwork_size;
    if (econd) {
        if (full) lwork_size = 2 * n;
        else lwork_size = 2 * ((bsn > st) ? bsn : st);
    } else {
        lwork_size = 1;
    }
    
    f64 *dwork = (f64*)malloc((lwork_size > 0 ? lwork_size : 1) * sizeof(f64));
    if (!dwork) {
         Py_DECREF(ipar_array); Py_DECREF(r_array); Py_DECREF(sdiag_array);
         Py_DECREF(s_array); Py_DECREF(b_array); Py_DECREF(ranks_array);
         return PyErr_NoMemory();
    }
    
    i32 ldr = (i32)PyArray_DIM(r_array, 0);
    i32 lds = (i32)PyArray_DIM(s_array, 0);
    if (lds < 1) lds = 1; /* Avoid zero dim */
    
    f64 *r_data = (f64*)PyArray_DATA(r_array);
    f64 *sdiag_data = (f64*)PyArray_DATA(sdiag_array);
    f64 *s_data = (f64*)PyArray_DATA(s_array);
    f64 *b_data = (f64*)PyArray_DATA(b_array);
    i32 *ranks_data = (i32*)PyArray_DATA(ranks_array);
    
    nf01br(cond, uplo, trans, n, ipar_data, lipar, r_data, ldr, sdiag_data, 
           s_data, lds, b_data, ranks_data, tol, dwork, lwork_size, &info);
           
    free(dwork);
    
    PyArray_ResolveWritebackIfCopy(b_array);
    PyArray_ResolveWritebackIfCopy(ranks_array);
    PyArray_ResolveWritebackIfCopy(r_array); // Potentially modified if UPLO=L
    PyArray_ResolveWritebackIfCopy(sdiag_array);
    PyArray_ResolveWritebackIfCopy(s_array);
    
    Py_DECREF(ipar_array); 
    Py_DECREF(r_array); 
    Py_DECREF(sdiag_array); 
    Py_DECREF(s_array); 
    
    /* Return modified b and ranks */
    PyObject *result = Py_BuildValue("OOi", b_array, ranks_array, info);
    Py_DECREF(b_array); Py_DECREF(ranks_array);
    
    return result;
}

/* Python wrapper for nf01bs */
static PyObject* py_nf01bs(PyObject* self, PyObject* args) {
    i32 n;
    PyObject *ipar_obj, *j_obj, *e_obj;
    f64 fnorm;
    PyArrayObject *ipar_array, *j_array, *e_array;
    i32 info;

    if (!PyArg_ParseTuple(args, "iOdOO", &n, &ipar_obj, &fnorm, &j_obj, &e_obj)) {
        return NULL;
    }
    
    if (n < 0 || fnorm < 0) {
        PyErr_Format(PyExc_ValueError, "Invalid arguments");
        return NULL;
    }

    ipar_array = (PyArrayObject*)PyArray_FROM_OTF(ipar_obj, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    j_array = (PyArrayObject*)PyArray_FROM_OTF(j_obj, NPY_DOUBLE, NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    e_array = (PyArrayObject*)PyArray_FROM_OTF(e_obj, NPY_DOUBLE, NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    
    if (!ipar_array || !j_array || !e_array) {
         Py_XDECREF(ipar_array); Py_XDECREF(j_array); Py_XDECREF(e_array);
         return NULL;
    }
    
    i32 lipar = (i32)PyArray_SIZE(ipar_array);
    if (lipar < 4) {
         Py_DECREF(ipar_array); Py_DECREF(j_array); Py_DECREF(e_array);
         PyErr_SetString(PyExc_ValueError, "ipar must have length >= 4");
         return NULL;
    }
    
    i32 *ipar_data = (i32*)PyArray_DATA(ipar_array);
    i32 bn = ipar_data[1];
    i32 bsn = ipar_data[3];
    
    i32 ldj = (i32)PyArray_DIM(j_array, 0);
    f64 *j_data = (f64*)PyArray_DATA(j_array);
    f64 *e_data = (f64*)PyArray_DATA(e_array);
    
    npy_intp jnorms_dims[1] = {n};
    PyObject *jnorms_array = PyArray_SimpleNew(1, jnorms_dims, NPY_DOUBLE);
    
    npy_intp ipvt_dims[1] = {n};
    PyObject *ipvt_array = PyArray_SimpleNew(1, ipvt_dims, NPY_INT32);
    
    f64 gnorm;
    
    if (!jnorms_array || !ipvt_array) {
         Py_DECREF(ipar_array); Py_DECREF(j_array); Py_DECREF(e_array);
         Py_XDECREF(jnorms_array); Py_XDECREF(ipvt_array);
         return NULL;
    }
    
    f64 *jnorms_data = (f64*)PyArray_DATA((PyArrayObject*)jnorms_array);
    i32 *ipvt_data = (i32*)PyArray_DATA((PyArrayObject*)ipvt_array);
    
    i32 lwork_size;
    if (n == 0) lwork_size = 1;
    else if (bn <= 1 || bsn == 0) lwork_size = 4*n + 1;
    else {
        /* Conservative estimate for general case */
        i32 st = ipar_data[0];
        i32 bsm = ipar_data[2];
        i32 jwork = bsn + (3*bsn + 1 > st ? 3*bsn + 1 : st);
        if (bsm > bsn) jwork = (jwork > 4*st + 1) ? jwork : 4*st + 1;
        lwork_size = jwork;
    }
    
    f64 *dwork = (f64*)malloc((lwork_size > 0 ? lwork_size : 1) * sizeof(f64));
    if (!dwork) {
         Py_DECREF(ipar_array); Py_DECREF(j_array); Py_DECREF(e_array);
         Py_DECREF(jnorms_array); Py_DECREF(ipvt_array);
         return PyErr_NoMemory();
    }
    
    nf01bs(n, ipar_data, lipar, fnorm, j_data, &ldj, e_data, 
           jnorms_data, &gnorm, ipvt_data, dwork, lwork_size, &info);
           
    free(dwork);
    
    PyArray_ResolveWritebackIfCopy(j_array);
    PyArray_ResolveWritebackIfCopy(e_array);
    Py_DECREF(ipar_array);
    
    PyObject *result = Py_BuildValue("OOOdOi", j_array, e_array, jnorms_array, gnorm, ipvt_array, info);
    
    Py_DECREF(j_array);
    Py_DECREF(e_array);
    Py_DECREF(jnorms_array);
    Py_DECREF(ipvt_array);
    
    return result;
}

/* Python wrapper for nf01ay */
static PyObject* py_nf01ay(PyObject* self, PyObject* args) {
    i32 nsmp, nz, l;
    PyObject *ipar_obj, *wb_obj, *z_obj;
    PyArrayObject *ipar_array, *wb_array, *z_array;
    i32 info;

    if (!PyArg_ParseTuple(args, "iiiOOO", &nsmp, &nz, &l, &ipar_obj, &wb_obj, &z_obj)) {
        return NULL;
    }
    
    if (nsmp < 0 || nz < 0 || l < 0) {
        PyErr_Format(PyExc_ValueError, "Dimensions must be non-negative");
        return NULL;
    }

    ipar_array = (PyArrayObject*)PyArray_FROM_OTF(ipar_obj, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    wb_array = (PyArrayObject*)PyArray_FROM_OTF(wb_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    z_array = (PyArrayObject*)PyArray_FROM_OTF(z_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
    
    if (!ipar_array || !wb_array || !z_array) {
         Py_XDECREF(ipar_array); Py_XDECREF(wb_array); Py_XDECREF(z_array);
         return NULL;
    }
    
    i32 lipar = (i32)PyArray_SIZE(ipar_array);
    i32 lwb = (i32)PyArray_SIZE(wb_array);
    i32 ldz = (i32)PyArray_DIM(z_array, 0);
    
    if (ldz < 1) ldz = 1;
    
    i32 *ipar_data = (i32*)PyArray_DATA(ipar_array);
    if (lipar < 1) {
         Py_DECREF(ipar_array); Py_DECREF(wb_array); Py_DECREF(z_array);
         PyErr_SetString(PyExc_ValueError, "ipar must have length >= 1");
         return NULL;
    }
    i32 nn = ipar_data[0];
    
    i32 block_size = 64;
    if (block_size > nsmp) block_size = nsmp;
    if (block_size < 2) block_size = 2;
    
    i32 ldwork_target = nn * (block_size + 1);
    if (ldwork_target < 2*nn) ldwork_target = 2*nn;
    
    f64 *dwork = (f64*)malloc(ldwork_target * sizeof(f64));
    if (!dwork) {
         Py_DECREF(ipar_array); Py_DECREF(wb_array); Py_DECREF(z_array);
         return PyErr_NoMemory();
    }
    
    npy_intp y_dims[2] = {nsmp, l};
    npy_intp y_strides[2] = {sizeof(f64), (nsmp > 0 ? nsmp : 1) * sizeof(f64)};
    PyObject *y_array = PyArray_New(&PyArray_Type, 2, y_dims, NPY_DOUBLE, y_strides, NULL, 0, NPY_ARRAY_FARRAY, NULL);
    
    if (!y_array) {
         free(dwork);
         Py_DECREF(ipar_array); Py_DECREF(wb_array); Py_DECREF(z_array);
         return NULL;
    }
    
    i32 ldy = (nsmp > 0) ? nsmp : 1;
    f64 *y_data = (f64*)PyArray_DATA((PyArrayObject*)y_array);
    f64 *wb_data = (f64*)PyArray_DATA(wb_array);
    f64 *z_data = (f64*)PyArray_DATA(z_array);
    
    nf01ay(nsmp, nz, l, ipar_data, lipar, wb_data, lwb, z_data, ldz, y_data, ldy, dwork, ldwork_target, &info);
    
    free(dwork);
    Py_DECREF(ipar_array); 
    Py_DECREF(wb_array); 
    Py_DECREF(z_array);
    
    PyObject *result = Py_BuildValue("Oi", y_array, info);
    Py_DECREF(y_array);
    
    return result;
}

/* Python wrapper for nf01by */
static PyObject* py_nf01by(PyObject* self, PyObject* args) {
    const char *cjte;
    i32 nsmp, nz, l;
    PyObject *ipar_obj, *wb_obj, *z_obj, *e_obj;
    PyArrayObject *ipar_array, *wb_array, *z_array, *e_array;
    i32 info;

    if (!PyArg_ParseTuple(args, "siiiOOOO", &cjte, &nsmp, &nz, &l, &ipar_obj, &wb_obj, &z_obj, &e_obj)) {
        return NULL;
    }
    
    if (nsmp < 0 || nz < 0 || l != 1) {
        PyErr_Format(PyExc_ValueError, "Dimensions invalid");
        return NULL;
    }

    ipar_array = (PyArrayObject*)PyArray_FROM_OTF(ipar_obj, NPY_INT32, NPY_ARRAY_IN_ARRAY);
    wb_array = (PyArrayObject*)PyArray_FROM_OTF(wb_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    z_array = (PyArrayObject*)PyArray_FROM_OTF(z_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
    e_array = (PyArrayObject*)PyArray_FROM_OTF(e_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    
    if (!ipar_array || !wb_array || !z_array || !e_array) {
         Py_XDECREF(ipar_array); Py_XDECREF(wb_array); Py_XDECREF(z_array); Py_XDECREF(e_array);
         return NULL;
    }
    
    i32 lipar = (i32)PyArray_SIZE(ipar_array);
    i32 lwb = (i32)PyArray_SIZE(wb_array);
    i32 ldz = (i32)PyArray_DIM(z_array, 0);
    
    if (ldz < 1) ldz = 1;
    
    i32 *ipar_data = (i32*)PyArray_DATA(ipar_array);
    if (lipar < 1) {
         Py_DECREF(ipar_array); Py_DECREF(wb_array); Py_DECREF(z_array); Py_DECREF(e_array);
         PyErr_SetString(PyExc_ValueError, "ipar must have length >= 1");
         return NULL;
    }
    i32 nn = ipar_data[0];
    i32 nwb = nn * (nz + 2) + 1;
    
    i32 ldj = (nsmp > 0) ? nsmp : 1;
    
    npy_intp j_dims[2] = {nsmp, nwb};
    npy_intp j_strides[2] = {sizeof(f64), ldj * sizeof(f64)};
    PyObject *j_array = PyArray_New(&PyArray_Type, 2, j_dims, NPY_DOUBLE, j_strides, NULL, 0, NPY_ARRAY_FARRAY, NULL);
    
    npy_intp jte_dims[1] = {nwb};
    PyObject *jte_array = PyArray_SimpleNew(1, jte_dims, NPY_DOUBLE);
    
    if (!j_array || !jte_array) {
         Py_XDECREF(j_array); Py_XDECREF(jte_array);
         Py_DECREF(ipar_array); Py_DECREF(wb_array); Py_DECREF(z_array); Py_DECREF(e_array);
         return NULL;
    }
    
    i32 ldwork = 2 * nn;
    if (ldwork < 1) ldwork = 1;
    
    f64 *dwork = (f64*)malloc(ldwork * sizeof(f64));
    if (!dwork) {
         Py_DECREF(j_array); Py_DECREF(jte_array);
         Py_DECREF(ipar_array); Py_DECREF(wb_array); Py_DECREF(z_array); Py_DECREF(e_array);
         return PyErr_NoMemory();
    }
    
    f64 *wb_data = (f64*)PyArray_DATA(wb_array);
    f64 *z_data = (f64*)PyArray_DATA(z_array);
    f64 *e_data = (f64*)PyArray_DATA(e_array);
    f64 *j_data = (f64*)PyArray_DATA((PyArrayObject*)j_array);
    f64 *jte_data = (f64*)PyArray_DATA((PyArrayObject*)jte_array);
    
    nf01by(cjte, nsmp, nz, l, ipar_data, lipar, wb_data, lwb, z_data, ldz, e_data, 
           j_data, ldj, jte_data, dwork, ldwork, &info);
           
    free(dwork);
    Py_DECREF(ipar_array); Py_DECREF(wb_array); Py_DECREF(z_array); Py_DECREF(e_array);
    
    PyObject *result = Py_BuildValue("OOi", j_array, jte_array, info);
    Py_DECREF(j_array); Py_DECREF(jte_array);
    
    return result;
}

/* Python wrapper for mb04ow */
static PyObject* py_mb04ow(PyObject* self, PyObject* args, PyObject* kwargs) {
    i32 m, n, p, incx = 1, incd = 1;
    PyObject *a_obj, *t_obj, *x_obj, *b_obj, *c_obj, *d_obj;
    PyArrayObject *a_array, *t_array, *x_array, *b_array, *c_array, *d_array;

    static char *kwlist[] = {"m", "n", "p", "a", "t", "x", "b", "c", "d", "incx", "incd", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "iiiOOOOOO|ii", kwlist,
                                     &m, &n, &p, &a_obj, &t_obj, &x_obj, 
                                     &b_obj, &c_obj, &d_obj, &incx, &incd)) {
        return NULL;
    }
    
    if (m < 0 || n < 0 || p < 0) {
        PyErr_Format(PyExc_ValueError, "Dimensions must be non-negative");
        return NULL;
    }
    
    a_array = (PyArrayObject*)PyArray_FROM_OTF(a_obj, NPY_DOUBLE, NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    t_array = (PyArrayObject*)PyArray_FROM_OTF(t_obj, NPY_DOUBLE, NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    x_array = (PyArrayObject*)PyArray_FROM_OTF(x_obj, NPY_DOUBLE, NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    b_array = (PyArrayObject*)PyArray_FROM_OTF(b_obj, NPY_DOUBLE, NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    c_array = (PyArrayObject*)PyArray_FROM_OTF(c_obj, NPY_DOUBLE, NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    d_array = (PyArrayObject*)PyArray_FROM_OTF(d_obj, NPY_DOUBLE, NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    
    if (!a_array || !t_array || !x_array || !b_array || !c_array || !d_array) {
         Py_XDECREF(a_array); Py_XDECREF(t_array); Py_XDECREF(x_array);
         Py_XDECREF(b_array); Py_XDECREF(c_array); Py_XDECREF(d_array);
         return NULL;
    }
    
    i32 lda = (i32)PyArray_DIM(a_array, 0);
    i32 ldt = (i32)PyArray_DIM(t_array, 0);
    i32 ldb = (i32)PyArray_DIM(b_array, 0);
    i32 ldc = (i32)PyArray_DIM(c_array, 0);
    
    if (lda < 1) lda = 1;
    if (ldt < 1) ldt = 1;
    if (ldb < 1) ldb = 1;
    if (ldc < 1) ldc = 1;
    
    f64 *a_data = (f64*)PyArray_DATA(a_array);
    f64 *t_data = (f64*)PyArray_DATA(t_array);
    f64 *x_data = (f64*)PyArray_DATA(x_array);
    f64 *b_data = (f64*)PyArray_DATA(b_array);
    f64 *c_data = (f64*)PyArray_DATA(c_array);
    f64 *d_data = (f64*)PyArray_DATA(d_array);
    
    mb04ow(m, n, p, a_data, lda, t_data, ldt, x_data, incx, b_data, ldb, c_data, ldc, d_data, incd);
    
    PyArray_ResolveWritebackIfCopy(a_array);
    PyArray_ResolveWritebackIfCopy(t_array);
    PyArray_ResolveWritebackIfCopy(x_array);
    PyArray_ResolveWritebackIfCopy(b_array);
    PyArray_ResolveWritebackIfCopy(c_array);
    PyArray_ResolveWritebackIfCopy(d_array);
    
    PyObject *result = Py_BuildValue("OOOOOO", a_array, t_array, x_array, b_array, c_array, d_array);
    
    Py_DECREF(a_array);
    Py_DECREF(t_array);
    Py_DECREF(x_array);
    Py_DECREF(b_array);
    Py_DECREF(c_array);
    Py_DECREF(d_array);

    return result;
}

static PyObject* py_mb03ud(PyObject* self, PyObject* args, PyObject* kwargs) {
    static char *kwlist[] = {"n", "a", "jobq", "jobp", "ldwork", NULL};
    i32 n;
    PyObject *a_obj;
    char *jobq_str = "N";
    char *jobp_str = "N";
    i32 ldwork = 0;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "iO|ssi", kwlist,
                                     &n, &a_obj, &jobq_str, &jobp_str, &ldwork)) {
        return NULL;
    }

    char jobq = jobq_str[0];
    char jobp = jobp_str[0];

    PyArrayObject *a_array = (PyArrayObject*)PyArray_FROM_OTF(a_obj, NPY_DOUBLE,
                                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (a_array == NULL) return NULL;

    npy_intp *a_dims = PyArray_DIMS(a_array);
    i32 lda = (i32)a_dims[0];

    bool wantq = (jobq == 'V' || jobq == 'v');
    bool wantp = (jobp == 'V' || jobp == 'v');

    i32 ldq = wantq ? (n > 0 ? n : 1) : 1;

    f64 *q = NULL;
    PyArrayObject *q_array = NULL;
    if (wantq) {
        npy_intp q_dims[2] = {n, n};
        npy_intp q_strides[2] = {sizeof(f64), n * sizeof(f64)};
        q = (f64*)calloc(n * n, sizeof(f64));
        if (q == NULL) {
            Py_DECREF(a_array);
            PyErr_SetString(PyExc_MemoryError, "Failed to allocate Q");
            return NULL;
        }
        q_array = (PyArrayObject*)PyArray_New(&PyArray_Type, 2, q_dims, NPY_DOUBLE,
                                              q_strides, q, 0, NPY_ARRAY_FARRAY, NULL);
        if (q_array == NULL) {
            free(q);
            Py_DECREF(a_array);
            return NULL;
        }
        PyArray_ENABLEFLAGS(q_array, NPY_ARRAY_OWNDATA);
    }

    f64 *sv = (f64*)calloc(n > 0 ? n : 1, sizeof(f64));
    if (sv == NULL) {
        Py_XDECREF(q_array);
        Py_DECREF(a_array);
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate sv");
        return NULL;
    }

    i32 minwork = n > 0 ? 5 * n : 1;
    if (ldwork == 0) ldwork = minwork;

    f64 *dwork = (f64*)calloc(ldwork, sizeof(f64));
    if (dwork == NULL) {
        free(sv);
        Py_XDECREF(q_array);
        Py_DECREF(a_array);
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate workspace");
        return NULL;
    }

    f64 *a_data = (f64*)PyArray_DATA(a_array);
    f64 *q_data = wantq ? q : NULL;

    i32 info;
    mb03ud(jobq, jobp, n, a_data, lda, q_data, ldq, sv, dwork, ldwork, &info);

    free(dwork);

    PyArray_ResolveWritebackIfCopy(a_array);

    if (info < 0) {
        free(sv);
        Py_XDECREF(q_array);
        Py_DECREF(a_array);
        PyErr_Format(PyExc_ValueError, "MB03UD: invalid parameter at position %d", -info);
        return NULL;
    }

    npy_intp sv_dims[1] = {n};
    PyArrayObject *sv_array = (PyArrayObject*)PyArray_SimpleNewFromData(1, sv_dims, NPY_DOUBLE, sv);
    if (sv_array == NULL) {
        free(sv);
        Py_XDECREF(q_array);
        Py_DECREF(a_array);
        return NULL;
    }
    PyArray_ENABLEFLAGS(sv_array, NPY_ARRAY_OWNDATA);

    PyObject *p_array = wantp ? (PyObject*)a_array : Py_None;
    PyObject *q_result = wantq ? (PyObject*)q_array : Py_None;

    if (!wantp) Py_INCREF(Py_None);
    if (!wantq) Py_INCREF(Py_None);

    PyObject *result = Py_BuildValue("OOOi", sv_array, p_array, q_result, info);

    Py_DECREF(sv_array);
    if (wantq) Py_DECREF(q_array);
    if (wantp) Py_DECREF(a_array);

    return result;
}

static PyObject* py_mb04id(PyObject* self, PyObject* args, PyObject* kwargs) {
    static char *kwlist[] = {"n", "m", "p", "a", "b", "l", "ldwork", NULL};
    i32 n, m, p, l = 0;
    i32 ldwork = 0;
    PyObject *a_obj, *b_obj = NULL;
    PyArrayObject *a_array = NULL, *b_array = NULL;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "iiiO|Oii", kwlist,
                                     &n, &m, &p, &a_obj, &b_obj, &l, &ldwork)) {
        return NULL;
    }

    a_array = (PyArrayObject*)PyArray_FROM_OTF(a_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (a_array == NULL) return NULL;

    npy_intp *a_dims = PyArray_DIMS(a_array);
    i32 lda = (i32)a_dims[0];

    i32 ldb = n > 0 ? n : 1;
    bool has_b = (b_obj != NULL && l > 0);

    if (has_b) {
        b_array = (PyArrayObject*)PyArray_FROM_OTF(b_obj, NPY_DOUBLE,
                                                   NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
        if (b_array == NULL) {
            Py_DECREF(a_array);
            return NULL;
        }
        npy_intp *b_dims = PyArray_DIMS(b_array);
        ldb = (i32)b_dims[0];
    }

    i32 minwork = 1;
    if (m > 1 && m - 1 > minwork) minwork = m - 1;
    if (m > p && m - p > minwork) minwork = m - p;
    if (l > minwork) minwork = l;

    if (ldwork == 0) ldwork = minwork;

    f64 *dwork = (f64*)calloc(ldwork > 0 ? ldwork : 1, sizeof(f64));
    if (dwork == NULL) {
        Py_XDECREF(b_array);
        Py_DECREF(a_array);
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate workspace");
        return NULL;
    }

    i32 minval = (n < m ? n : m);
    f64 *tau = (f64*)calloc(minval > 0 ? minval : 1, sizeof(f64));
    if (tau == NULL) {
        free(dwork);
        Py_XDECREF(b_array);
        Py_DECREF(a_array);
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate tau");
        return NULL;
    }

    f64 *a_data = (f64*)PyArray_DATA(a_array);
    f64 *b_data = has_b ? (f64*)PyArray_DATA(b_array) : NULL;
    f64 dummy_b = 0.0;
    if (!has_b) {
        b_data = &dummy_b;
        ldb = 1;
    }

    i32 info;
    mb04id(n, m, p, l, a_data, lda, b_data, ldb, tau, dwork, ldwork, &info);

    f64 optimal_work = dwork[0];
    free(dwork);

    PyArray_ResolveWritebackIfCopy(a_array);
    if (has_b) {
        PyArray_ResolveWritebackIfCopy(b_array);
    }

    if (info != 0) {
        free(tau);
        Py_XDECREF(b_array);
        Py_DECREF(a_array);
        PyErr_Format(PyExc_ValueError, "mb04id failed with info=%d", info);
        return NULL;
    }

    npy_intp tau_dims[1] = {minval > 0 ? minval : 1};
    PyObject *tau_array = PyArray_SimpleNewFromData(1, tau_dims, NPY_DOUBLE, tau);
    if (tau_array == NULL) {
        free(tau);
        Py_XDECREF(b_array);
        Py_DECREF(a_array);
        return NULL;
    }
    PyArray_ENABLEFLAGS((PyArrayObject*)tau_array, NPY_ARRAY_OWNDATA);

    PyObject *result;
    if (ldwork == -1) {
        if (has_b) {
            result = Py_BuildValue("OOOid", a_array, b_array, tau_array, info, optimal_work);
        } else {
            result = Py_BuildValue("OOid", a_array, tau_array, info, optimal_work);
        }
    } else {
        if (has_b) {
            result = Py_BuildValue("OOOi", a_array, b_array, tau_array, info);
        } else {
            result = Py_BuildValue("OOi", a_array, tau_array, info);
        }
    }

    Py_XDECREF(b_array);
    Py_DECREF(a_array);
    Py_DECREF(tau_array);

    return result;
}

static PyObject* py_mb04iy(PyObject* self, PyObject* args, PyObject* kwargs) {
    static char *kwlist[] = {"side", "trans", "a", "tau", "c", "p", NULL};
    const char *side, *trans;
    i32 p = 0;
    PyObject *a_obj, *tau_obj, *c_obj;
    PyArrayObject *a_array, *tau_array, *c_array;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ssOOO|i", kwlist,
                                     &side, &trans, &a_obj, &tau_obj, &c_obj, &p)) {
        return NULL;
    }

    a_array = (PyArrayObject*)PyArray_FROM_OTF(a_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (a_array == NULL) return NULL;

    tau_array = (PyArrayObject*)PyArray_FROM_OTF(tau_obj, NPY_DOUBLE, NPY_ARRAY_FARRAY);
    if (tau_array == NULL) {
        Py_DECREF(a_array);
        return NULL;
    }

    c_array = (PyArrayObject*)PyArray_FROM_OTF(c_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (c_array == NULL) {
        Py_DECREF(a_array);
        Py_DECREF(tau_array);
        return NULL;
    }

    npy_intp *a_dims = PyArray_DIMS(a_array);
    npy_intp *c_dims = PyArray_DIMS(c_array);
    npy_intp *tau_dims = PyArray_DIMS(tau_array);

    i32 lda = (i32)a_dims[0];
    i32 ldc = (i32)c_dims[0];
    i32 n = (i32)c_dims[0];
    i32 m = (i32)c_dims[1];
    i32 k = (i32)tau_dims[0];

    bool left = (side[0] == 'L' || side[0] == 'l');
    i32 ldwork = left ? (m > 0 ? m : 1) : (n > 0 ? n : 1);

    f64 *dwork = (f64*)calloc(ldwork, sizeof(f64));
    if (dwork == NULL) {
        Py_DECREF(a_array);
        Py_DECREF(tau_array);
        Py_DECREF(c_array);
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate workspace");
        return NULL;
    }

    f64 *a_data = (f64*)PyArray_DATA(a_array);
    f64 *tau_data = (f64*)PyArray_DATA(tau_array);
    f64 *c_data = (f64*)PyArray_DATA(c_array);

    i32 info;
    mb04iy(side, trans, n, m, k, p, a_data, lda, tau_data, c_data, ldc, dwork, ldwork, &info);

    free(dwork);

    // Resolve WRITEBACKIFCOPY for 'a' array before DECREF
    if (PyArray_ResolveWritebackIfCopy(a_array) < 0) {
        Py_DECREF(a_array);
        Py_DECREF(tau_array);
        Py_DECREF(c_array);
        PyErr_SetString(PyExc_RuntimeError, "Failed to resolve WRITEBACKIFCOPY for array 'a'");
        return NULL;
    }

    // Resolve WRITEBACKIFCOPY for 'c' array
    if (PyArray_ResolveWritebackIfCopy(c_array) < 0) {
        Py_DECREF(a_array);
        Py_DECREF(tau_array);
        Py_DECREF(c_array);
        PyErr_SetString(PyExc_RuntimeError, "Failed to resolve WRITEBACKIFCOPY for array 'c'");
        return NULL;
    }

    if (info != 0) {
        Py_DECREF(a_array);
        Py_DECREF(tau_array);
        Py_DECREF(c_array);
        PyErr_Format(PyExc_ValueError, "mb04iy failed with info=%d", info);
        return NULL;
    }

    Py_DECREF(a_array);
    Py_DECREF(tau_array);

    PyObject *result = Py_BuildValue("Oi", c_array, info);
    Py_DECREF(c_array);
    return result;
}

static PyObject* py_mb04oy(PyObject* self, PyObject* args) {
    i32 m, n;
    f64 tau;
    PyObject *v_obj, *a_obj, *b_obj;
    PyArrayObject *v_array, *a_array, *b_array;

    if (!PyArg_ParseTuple(args, "iiOdOO", &m, &n, &v_obj, &tau, &a_obj, &b_obj)) {
        return NULL;
    }

    v_array = (PyArrayObject*)PyArray_FROM_OTF(v_obj, NPY_DOUBLE, NPY_ARRAY_FARRAY);
    if (v_array == NULL) return NULL;

    a_array = (PyArrayObject*)PyArray_FROM_OTF(a_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (a_array == NULL) {
        Py_DECREF(v_array);
        return NULL;
    }

    b_array = (PyArrayObject*)PyArray_FROM_OTF(b_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (b_array == NULL) {
        Py_DECREF(v_array);
        Py_DECREF(a_array);
        return NULL;
    }

    npy_intp *a_dims = PyArray_DIMS(a_array);
    npy_intp *b_dims = PyArray_DIMS(b_array);

    i32 lda = (i32)a_dims[0];
    i32 ldb = (i32)b_dims[0];

    f64 *dwork = NULL;
    if (m + 1 >= 11) {
        dwork = (f64*)calloc(n > 0 ? n : 1, sizeof(f64));
        if (dwork == NULL) {
            Py_DECREF(v_array);
            Py_DECREF(a_array);
            Py_DECREF(b_array);
            PyErr_SetString(PyExc_MemoryError, "Failed to allocate workspace");
            return NULL;
        }
    }

    f64 *v_data = (f64*)PyArray_DATA(v_array);
    f64 *a_data = (f64*)PyArray_DATA(a_array);
    f64 *b_data = (f64*)PyArray_DATA(b_array);

    SLC_MB04OY(m, n, v_data, tau, a_data, lda, b_data, ldb, dwork);

    if (dwork != NULL) {
        free(dwork);
    }

    Py_DECREF(v_array);

    PyObject *result = Py_BuildValue("OO", a_array, b_array);
    Py_DECREF(a_array);
    Py_DECREF(b_array);
    return result;
}

/* Python wrapper for mb04ny - Apply Householder reflector to [A B] from right */
static PyObject* py_mb04ny(PyObject* self, PyObject* args) {
    i32 m, n, incv;
    f64 tau;
    PyObject *v_obj, *a_obj, *b_obj;
    PyArrayObject *v_array, *a_array, *b_array;

    if (!PyArg_ParseTuple(args, "iiOidOO", &m, &n, &v_obj, &incv, &tau, &a_obj, &b_obj)) {
        return NULL;
    }

    v_array = (PyArrayObject*)PyArray_FROM_OTF(v_obj, NPY_DOUBLE, NPY_ARRAY_FARRAY);
    if (v_array == NULL) return NULL;

    a_array = (PyArrayObject*)PyArray_FROM_OTF(a_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (a_array == NULL) {
        Py_DECREF(v_array);
        return NULL;
    }

    b_array = (PyArrayObject*)PyArray_FROM_OTF(b_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (b_array == NULL) {
        Py_DECREF(v_array);
        Py_DECREF(a_array);
        return NULL;
    }

    npy_intp *a_dims = PyArray_DIMS(a_array);
    npy_intp *b_dims = PyArray_DIMS(b_array);

    i32 lda = (m > 0) ? (i32)a_dims[0] : 1;
    i32 ldb = (m > 0) ? (i32)b_dims[0] : 1;

    f64 *dwork = NULL;
    if (n + 1 >= 11) {
        dwork = (f64*)calloc(m > 0 ? m : 1, sizeof(f64));
        if (dwork == NULL) {
            Py_DECREF(v_array);
            Py_DECREF(a_array);
            Py_DECREF(b_array);
            PyErr_SetString(PyExc_MemoryError, "Failed to allocate workspace");
            return NULL;
        }
    }

    f64 *v_data = (f64*)PyArray_DATA(v_array);
    f64 *a_data = (f64*)PyArray_DATA(a_array);
    f64 *b_data = (f64*)PyArray_DATA(b_array);

    SLC_MB04NY(m, n, v_data, incv, tau, a_data, lda, b_data, ldb, dwork);

    if (dwork != NULL) {
        free(dwork);
    }

    Py_DECREF(v_array);

    PyObject *result = Py_BuildValue("OO", a_array, b_array);
    Py_DECREF(a_array);
    Py_DECREF(b_array);
    return result;
}

/* Python wrapper for ib01od */
static PyObject* py_ib01od(PyObject* self, PyObject* args) {
    const char *ctrl_str;
    i32 nobr, l;
    f64 tol;
    PyObject *sv_obj;
    PyArrayObject *sv_array;
    i32 n, iwarn, info;

    if (!PyArg_ParseTuple(args, "siiOd", &ctrl_str, &nobr, &l, &sv_obj, &tol)) {
        return NULL;
    }

    char ctrl = ctrl_str[0];

    /* Validate CTRL parameter */
    if (ctrl != 'C' && ctrl != 'c' && ctrl != 'N' && ctrl != 'n') {
        PyErr_SetString(PyExc_ValueError, "CTRL must be 'C' or 'N'");
        return NULL;
    }

    /* Validate NOBR */
    if (nobr <= 0) {
        PyErr_SetString(PyExc_ValueError, "NOBR must be positive");
        return NULL;
    }

    /* Validate L */
    if (l <= 0) {
        PyErr_SetString(PyExc_ValueError, "L must be positive");
        return NULL;
    }

    /* Convert SV array */
    sv_array = (PyArrayObject*)PyArray_FROM_OTF(sv_obj, NPY_DOUBLE,
                                                NPY_ARRAY_IN_FARRAY);
    if (sv_array == NULL) {
        return NULL;
    }

    /* Validate SV size */
    npy_intp sv_size = PyArray_SIZE(sv_array);
    i32 lnobr = l * nobr;
    if (sv_size < lnobr) {
        Py_DECREF(sv_array);
        PyErr_SetString(PyExc_ValueError, "SV must have at least L*NOBR elements");
        return NULL;
    }

    f64 *sv = (f64*)PyArray_DATA(sv_array);

    /* Call C function */
    SLC_IB01OD(ctrl, nobr, l, sv, &n, tol, &iwarn, &info);

    Py_DECREF(sv_array);

    if (info < 0) {
        PyErr_Format(PyExc_ValueError, "Parameter %d had an illegal value", -info);
        return NULL;
    }

    return Py_BuildValue("iii", n, iwarn, info);
}

/* Python wrapper for ib01nd */
static PyObject* py_ib01nd(PyObject* self, PyObject* args) {
    const char *meth_str, *jobd_str;
    i32 nobr, m, l;
    f64 tol;
    PyObject *r_obj;
    PyArrayObject *r_array;
    i32 iwarn, info;

    if (!PyArg_ParseTuple(args, "ssiiiOd", &meth_str, &jobd_str, &nobr, &m, &l,
                          &r_obj, &tol)) {
        return NULL;
    }

    char meth = meth_str[0];
    char jobd = jobd_str[0];

    /* Validate METH */
    if (meth != 'M' && meth != 'm' && meth != 'N' && meth != 'n') {
        PyErr_SetString(PyExc_ValueError, "METH must be 'M' or 'N'");
        return NULL;
    }

    /* Validate JOBD for MOESP */
    if ((meth == 'M' || meth == 'm') &&
        jobd != 'M' && jobd != 'm' && jobd != 'N' && jobd != 'n') {
        PyErr_SetString(PyExc_ValueError, "JOBD must be 'M' or 'N'");
        return NULL;
    }

    /* Validate dimensions */
    if (nobr <= 0) {
        PyErr_SetString(PyExc_ValueError, "NOBR must be positive");
        return NULL;
    }
    if (m < 0) {
        PyErr_SetString(PyExc_ValueError, "M must be non-negative");
        return NULL;
    }
    if (l <= 0) {
        PyErr_SetString(PyExc_ValueError, "L must be positive");
        return NULL;
    }

    /* Convert R array */
    r_array = (PyArrayObject*)PyArray_FROM_OTF(r_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (r_array == NULL) {
        return NULL;
    }

    /* Get dimensions */
    int ndim = PyArray_NDIM(r_array);
    if (ndim != 2) {
        Py_DECREF(r_array);
        PyErr_SetString(PyExc_ValueError, "R must be a 2D array");
        return NULL;
    }

    npy_intp *r_dims = PyArray_DIMS(r_array);
    i32 ldr = (i32)r_dims[0];
    i32 nr = 2 * (m + l) * nobr;

    /* Validate R size */
    if (ldr < nr || r_dims[1] < nr) {
        Py_DECREF(r_array);
        PyErr_SetString(PyExc_ValueError, "R must be at least 2*(m+l)*nobr x 2*(m+l)*nobr");
        return NULL;
    }

    f64 *r_data = (f64*)PyArray_DATA(r_array);

    /* Allocate output singular values */
    i32 lnobr = l * nobr;
    npy_intp sv_dims[1] = {lnobr};
    PyArrayObject *sv_array = (PyArrayObject*)PyArray_SimpleNew(1, sv_dims, NPY_DOUBLE);
    if (sv_array == NULL) {
        Py_DECREF(r_array);
        return NULL;
    }
    f64 *sv = (f64*)PyArray_DATA(sv_array);

    /* Allocate workspace */
    i32 lmnobr = lnobr + m * nobr;
    i32 ldwork;
    bool moesp = (meth == 'M' || meth == 'm');
    bool jobdm = (jobd == 'M' || jobd == 'm');

    if (moesp) {
        if (jobdm) {
            i32 t1 = (2 * m - 1) * nobr;
            if (t1 < 0) t1 = 1;
            i32 t2 = lmnobr;
            i32 t3 = 5 * lnobr;
            ldwork = t1 > t2 ? t1 : t2;
            ldwork = ldwork > t3 ? ldwork : t3;
        } else {
            ldwork = 5 * lnobr;
        }
    } else {
        ldwork = 5 * lmnobr + 1;
    }
    ldwork = ldwork > 1 ? ldwork : 1;

    f64 *dwork = (f64*)malloc(ldwork * sizeof(f64));
    i32 *iwork = (i32*)malloc(lmnobr * sizeof(i32));
    if (dwork == NULL || iwork == NULL) {
        free(dwork);
        free(iwork);
        Py_DECREF(r_array);
        Py_DECREF(sv_array);
        PyErr_NoMemory();
        return NULL;
    }

    /* Call C function */
    SLC_IB01ND(meth, jobd, nobr, m, l, r_data, ldr, sv, tol,
               iwork, dwork, ldwork, &iwarn, &info);

    /* Get rcond values for N4SID */
    f64 rcond1 = 0.0, rcond2 = 0.0;
    if (!moesp) {
        rcond1 = dwork[1];
        rcond2 = dwork[2];
    }

    free(dwork);
    free(iwork);

    if (info < 0) {
        Py_DECREF(r_array);
        Py_DECREF(sv_array);
        PyErr_Format(PyExc_ValueError, "Parameter %d had an illegal value", -info);
        return NULL;
    }

    /* Build result: (r, sv, rcond1, rcond2, iwarn, info) */
    PyObject *result = Py_BuildValue("OOddii", r_array, sv_array,
                                      rcond1, rcond2, iwarn, info);
    Py_DECREF(r_array);
    Py_DECREF(sv_array);
    return result;
}

/* Python wrapper for ib01qd */
static PyObject* py_ib01qd(PyObject* self, PyObject* args) {
    const char *jobx0_str, *job_str;
    i32 n, m, l;
    f64 tol;
    PyObject *a_obj, *c_obj, *u_obj, *y_obj;
    PyArrayObject *a_array, *c_array, *u_array, *y_array;
    i32 iwarn, info;

    if (!PyArg_ParseTuple(args, "ssiiiOOOOd",
                          &jobx0_str, &job_str, &n, &m, &l,
                          &a_obj, &c_obj, &u_obj, &y_obj, &tol)) {
        return NULL;
    }

    char jobx0 = toupper((unsigned char)jobx0_str[0]);
    char job = toupper((unsigned char)job_str[0]);

    if (jobx0 != 'X' && jobx0 != 'N') {
        PyErr_SetString(PyExc_ValueError, "JOBX0 must be 'X' or 'N'");
        return NULL;
    }
    if (job != 'B' && job != 'D') {
        PyErr_SetString(PyExc_ValueError, "JOB must be 'B' or 'D'");
        return NULL;
    }
    if (n < 0) {
        PyErr_SetString(PyExc_ValueError, "N must be non-negative");
        return NULL;
    }
    if (m < 0) {
        PyErr_SetString(PyExc_ValueError, "M must be non-negative");
        return NULL;
    }
    if (l <= 0) {
        PyErr_SetString(PyExc_ValueError, "L must be positive");
        return NULL;
    }

    a_array = (PyArrayObject*)PyArray_FROM_OTF(a_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
    c_array = (PyArrayObject*)PyArray_FROM_OTF(c_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
    u_array = (PyArrayObject*)PyArray_FROM_OTF(u_obj, NPY_DOUBLE, NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    y_array = (PyArrayObject*)PyArray_FROM_OTF(y_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);

    if (!a_array || !c_array || !u_array || !y_array) {
        Py_XDECREF(a_array);
        Py_XDECREF(c_array);
        Py_XDECREF(u_array);
        Py_XDECREF(y_array);
        return NULL;
    }

    npy_intp *u_dims = PyArray_DIMS(u_array);
    npy_intp *y_dims = PyArray_DIMS(y_array);
    i32 nsmp = (i32)u_dims[0];
    i32 ldu = nsmp > 0 ? nsmp : 1;
    i32 ldy = nsmp > 0 ? nsmp : 1;
    i32 lda = n > 0 ? n : 1;
    i32 ldc = l;
    i32 ldb = (n > 0 && m > 0) ? n : 1;
    i32 ldd = (m > 0 && job == 'D') ? l : 1;

    bool withx0 = (jobx0 == 'X');
    bool withd = (job == 'D');
    i32 ncol = n * m + (withx0 ? n : 0);
    i32 minsmp = ncol;
    if (withd) {
        minsmp += m;
    } else if (!withx0) {
        minsmp += 1;
    }

    if (nsmp < minsmp) {
        Py_DECREF(a_array);
        Py_DECREF(c_array);
        Py_DECREF(u_array);
        Py_DECREF(y_array);
        PyErr_SetString(PyExc_ValueError, "NSMP too small for the problem dimensions");
        return NULL;
    }

    const f64 *a_data = (const f64*)PyArray_DATA(a_array);
    const f64 *c_data = (const f64*)PyArray_DATA(c_array);
    f64 *u_data = (f64*)PyArray_DATA(u_array);
    const f64 *y_data = (const f64*)PyArray_DATA(y_array);

    npy_intp x0_dims[1] = {n};
    PyArrayObject *x0_array = (PyArrayObject*)PyArray_SimpleNew(1, x0_dims, NPY_DOUBLE);

    npy_intp b_dims[2] = {n > 0 ? n : 1, m > 0 ? m : 1};
    npy_intp b_strides[2] = {sizeof(f64), ldb * sizeof(f64)};
    f64 *b_data = NULL;
    PyArrayObject *b_array = NULL;
    if (n > 0 && m > 0) {
        b_data = (f64*)calloc(n * m, sizeof(f64));
        if (!b_data) {
            Py_DECREF(a_array);
            Py_DECREF(c_array);
            Py_DECREF(u_array);
            Py_DECREF(y_array);
            Py_DECREF(x0_array);
            return PyErr_NoMemory();
        }
        b_array = (PyArrayObject*)PyArray_New(&PyArray_Type, 2, b_dims, NPY_DOUBLE,
                                               b_strides, b_data, 0, NPY_ARRAY_FARRAY, NULL);
        PyArray_ENABLEFLAGS(b_array, NPY_ARRAY_OWNDATA);
    } else {
        b_array = (PyArrayObject*)PyArray_SimpleNew(2, b_dims, NPY_DOUBLE);
        b_data = (f64*)PyArray_DATA(b_array);
    }

    npy_intp d_dims[2] = {l, m > 0 ? m : 1};
    npy_intp d_strides[2] = {sizeof(f64), ldd * sizeof(f64)};
    f64 *d_data = NULL;
    PyArrayObject *d_array = NULL;
    if (m > 0 && withd) {
        d_data = (f64*)calloc(l * m, sizeof(f64));
        if (!d_data) {
            Py_DECREF(a_array);
            Py_DECREF(c_array);
            Py_DECREF(u_array);
            Py_DECREF(y_array);
            Py_DECREF(x0_array);
            Py_DECREF(b_array);
            return PyErr_NoMemory();
        }
        d_array = (PyArrayObject*)PyArray_New(&PyArray_Type, 2, d_dims, NPY_DOUBLE,
                                               d_strides, d_data, 0, NPY_ARRAY_FARRAY, NULL);
        PyArray_ENABLEFLAGS(d_array, NPY_ARRAY_OWNDATA);
    } else {
        d_array = (PyArrayObject*)PyArray_SimpleNew(2, d_dims, NPY_DOUBLE);
        d_data = (f64*)PyArray_DATA(d_array);
    }

    if (!x0_array || !b_array || !d_array) {
        Py_DECREF(a_array);
        Py_DECREF(c_array);
        Py_DECREF(u_array);
        Py_DECREF(y_array);
        Py_XDECREF(x0_array);
        Py_XDECREF(b_array);
        Py_XDECREF(d_array);
        return NULL;
    }

    f64 *x0_data = (f64*)PyArray_DATA(x0_array);

    i32 nsmpl = nsmp * l;
    i32 iq = ncol + (withd ? m : (withx0 ? 0 : 1));
    iq = iq * l;
    i32 ncp1 = ncol + 1;
    i32 isize = nsmpl * ncp1;

    i32 ic = (n > 0 && withx0) ? (2 * n * n + n) : 0;
    i32 minwls = ncol * ncp1;
    if (withd) minwls += l * m * ncp1;

    i32 ia;
    if (m > 0 && withd) {
        ia = m + (2 * ncol > m ? 2 * ncol : m);
    } else {
        ia = 2 * ncol;
    }

    i32 itau = n * n * m + (ic > ia ? ic : ia);
    if (withx0) itau += l * n;

    i32 ldw2 = isize + (n + (ic > ia ? ic : ia));
    i32 t = 6 * ncol;
    if (n + (ic > ia ? ic : ia) < t) ldw2 = isize + t;

    i32 ldw3 = minwls + (iq * ncp1 + itau);
    if (iq * ncp1 + itau < 6 * ncol) ldw3 = minwls + 6 * ncol;

    if (m > 0 && withd) {
        i32 t2 = isize + 2 * m * m + 6 * m;
        if (t2 > ldw2) ldw2 = t2;
        t2 = minwls + 2 * m * m + 6 * m;
        if (t2 > ldw3) ldw3 = t2;
    }

    i32 ldwork = (ldw2 < ldw3) ? ldw2 : ldw3;
    if (ldwork < 2) ldwork = 2;
    if (m > 0 && withd && ldwork < 3) ldwork = 3;
    ldwork = ldw2;

    i32 liwork = n * m + (withx0 ? n : 0);
    if (withd && m > liwork) liwork = m;
    if (liwork < 1) liwork = 1;

    f64 *dwork = (f64*)malloc(ldwork * sizeof(f64));
    i32 *iwork = (i32*)malloc(liwork * sizeof(i32));
    if (!dwork || !iwork) {
        free(dwork);
        free(iwork);
        Py_DECREF(a_array);
        Py_DECREF(c_array);
        Py_DECREF(u_array);
        Py_DECREF(y_array);
        Py_DECREF(x0_array);
        Py_DECREF(b_array);
        Py_DECREF(d_array);
        return PyErr_NoMemory();
    }

    slicot_ib01qd(jobx0_str, job_str, n, m, l, nsmp,
                  a_data, lda, c_data, ldc,
                  u_data, ldu, y_data, ldy,
                  x0_data, b_data, ldb, d_data, ldd,
                  tol, iwork, dwork, ldwork, &iwarn, &info);

    f64 rcond_w2 = dwork[1];
    f64 rcond_u = (m > 0 && withd) ? dwork[2] : 1.0;

    free(dwork);
    free(iwork);

    PyArray_ResolveWritebackIfCopy(u_array);
    Py_DECREF(a_array);
    Py_DECREF(c_array);
    Py_DECREF(u_array);
    Py_DECREF(y_array);

    if (info < 0) {
        Py_DECREF(x0_array);
        Py_DECREF(b_array);
        Py_DECREF(d_array);
        PyErr_Format(PyExc_ValueError, "Parameter %d had an illegal value", -info);
        return NULL;
    }

    PyObject *result = Py_BuildValue("OOOddii", x0_array, b_array, d_array,
                                      rcond_w2, rcond_u, iwarn, info);
    Py_DECREF(x0_array);
    Py_DECREF(b_array);
    Py_DECREF(d_array);
    return result;
}

/* Python wrapper for ib01rd */
static PyObject* py_ib01rd(PyObject* self, PyObject* args) {
    const char *job_str;
    i32 n, m, l, nsmp;
    f64 tol;
    PyObject *a_obj, *b_obj, *c_obj, *d_obj, *u_obj, *y_obj;
    PyArrayObject *a_array, *b_array, *c_array, *d_array, *u_array, *y_array;
    i32 iwarn, info;

    if (!PyArg_ParseTuple(args, "siiiiOOOOOOd",
                          &job_str, &n, &m, &l, &nsmp,
                          &a_obj, &b_obj, &c_obj, &d_obj,
                          &u_obj, &y_obj, &tol)) {
        return NULL;
    }

    char job = toupper((unsigned char)job_str[0]);

    if (job != 'Z' && job != 'N') {
        PyErr_SetString(PyExc_ValueError, "JOB must be 'Z' or 'N'");
        return NULL;
    }
    if (n < 0) {
        PyErr_SetString(PyExc_ValueError, "N must be non-negative");
        return NULL;
    }
    if (m < 0) {
        PyErr_SetString(PyExc_ValueError, "M must be non-negative");
        return NULL;
    }
    if (l <= 0) {
        PyErr_SetString(PyExc_ValueError, "L must be positive");
        return NULL;
    }
    if (nsmp < n) {
        PyErr_SetString(PyExc_ValueError, "NSMP must be >= N");
        return NULL;
    }

    a_array = (PyArrayObject*)PyArray_FROM_OTF(a_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
    b_array = (PyArrayObject*)PyArray_FROM_OTF(b_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
    c_array = (PyArrayObject*)PyArray_FROM_OTF(c_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
    d_array = (PyArrayObject*)PyArray_FROM_OTF(d_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
    u_array = (PyArrayObject*)PyArray_FROM_OTF(u_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
    y_array = (PyArrayObject*)PyArray_FROM_OTF(y_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);

    if (!a_array || !b_array || !c_array || !d_array || !u_array || !y_array) {
        Py_XDECREF(a_array);
        Py_XDECREF(b_array);
        Py_XDECREF(c_array);
        Py_XDECREF(d_array);
        Py_XDECREF(u_array);
        Py_XDECREF(y_array);
        return NULL;
    }

    i32 lda = n > 0 ? n : 1;
    i32 ldb = (n > 0 && m > 0) ? n : 1;
    i32 ldc = l;
    bool withd = (job == 'N');
    i32 ldd = (withd && m > 0) ? l : 1;
    i32 ldu = (m > 0) ? nsmp : 1;
    i32 ldy = nsmp > 0 ? nsmp : 1;

    const f64 *a_data = (const f64*)PyArray_DATA(a_array);
    const f64 *b_data = (const f64*)PyArray_DATA(b_array);
    const f64 *c_data = (const f64*)PyArray_DATA(c_array);
    const f64 *d_data = (const f64*)PyArray_DATA(d_array);
    const f64 *u_data = (const f64*)PyArray_DATA(u_array);
    const f64 *y_data = (const f64*)PyArray_DATA(y_array);

    npy_intp x0_dims[1] = {n};
    PyArrayObject *x0_array = (PyArrayObject*)PyArray_SimpleNew(1, x0_dims, NPY_DOUBLE);
    if (!x0_array) {
        Py_DECREF(a_array);
        Py_DECREF(b_array);
        Py_DECREF(c_array);
        Py_DECREF(d_array);
        Py_DECREF(u_array);
        Py_DECREF(y_array);
        return NULL;
    }
    f64 *x0_data = (f64*)PyArray_DATA(x0_array);

    i32 nn = n * n;
    i32 nsmpl = nsmp * l;
    i32 iq = n * l;
    i32 ncp1 = n + 1;
    i32 isize = nsmpl * ncp1;
    i32 ic = 2 * nn;
    i32 minwls = n * ncp1;
    i32 itau_calc = ic + l * n;
    i32 ldw1 = isize + 2 * n + ((ic > 4 * n) ? ic : 4 * n);
    i32 ldw2 = minwls + 2 * n + ((iq * ncp1 + itau_calc > 4 * n) ? (iq * ncp1 + itau_calc) : 4 * n);
    i32 ldwork = (ldw1 > ldw2) ? ldw1 : ldw2;
    if (ldwork < 2) ldwork = 2;

    i32 liwork = (n > 0) ? n : 1;

    f64 *dwork = (f64*)malloc(ldwork * sizeof(f64));
    i32 *iwork = (i32*)malloc(liwork * sizeof(i32));
    if (!dwork || !iwork) {
        free(dwork);
        free(iwork);
        Py_DECREF(a_array);
        Py_DECREF(b_array);
        Py_DECREF(c_array);
        Py_DECREF(d_array);
        Py_DECREF(u_array);
        Py_DECREF(y_array);
        Py_DECREF(x0_array);
        return PyErr_NoMemory();
    }

    char job_c[2] = {job, '\0'};
    slicot_ib01rd(job_c, n, m, l, nsmp,
                  a_data, lda, b_data, ldb, c_data, ldc, d_data, ldd,
                  u_data, ldu, y_data, ldy,
                  x0_data, tol, iwork, dwork, ldwork, &iwarn, &info);

    f64 rcond = (n > 0) ? dwork[1] : 1.0;

    free(dwork);
    free(iwork);

    Py_DECREF(a_array);
    Py_DECREF(b_array);
    Py_DECREF(c_array);
    Py_DECREF(d_array);
    Py_DECREF(u_array);
    Py_DECREF(y_array);

    if (info < 0) {
        Py_DECREF(x0_array);
        PyErr_Format(PyExc_ValueError, "Parameter %d had an illegal value", -info);
        return NULL;
    }

    PyObject *result = Py_BuildValue("Odii", x0_array, rcond, iwarn, info);
    Py_DECREF(x0_array);
    return result;
}

/* Python wrapper for sb02mt */
static PyObject* py_sb02mt(PyObject* self, PyObject* args) {
    const char *jobg_str, *jobl_str, *fact_str, *uplo_str;
    i32 n, m;
    PyObject *a_obj, *b_obj, *q_obj, *r_obj, *l_obj, *g_obj;
    PyArrayObject *a_array = NULL, *b_array = NULL, *q_array = NULL;
    PyArrayObject *r_array = NULL, *l_array = NULL, *g_array = NULL;
    i32 oufact, info;

    if (!PyArg_ParseTuple(args, "ssssiiOOOOOO",
                          &jobg_str, &jobl_str, &fact_str, &uplo_str,
                          &n, &m, &a_obj, &b_obj, &q_obj, &r_obj, &l_obj, &g_obj)) {
        return NULL;
    }

    char jobg = toupper((unsigned char)jobg_str[0]);
    char jobl = toupper((unsigned char)jobl_str[0]);
    char fact = toupper((unsigned char)fact_str[0]);
    char uplo = toupper((unsigned char)uplo_str[0]);

    bool ljobg = (jobg == 'G');
    bool ljobl = (jobl == 'N');
    bool lfactc = (fact == 'C');
    bool lfactu = (fact == 'U');
    bool lnfact = (!lfactc && !lfactu);

    if (!ljobg && jobg != 'N') {
        PyErr_SetString(PyExc_ValueError, "Parameter 1 (JOBG) must be 'G' or 'N'");
        return NULL;
    }
    if (!ljobl && jobl != 'Z') {
        PyErr_SetString(PyExc_ValueError, "Parameter 2 (JOBL) must be 'Z' or 'N'");
        return NULL;
    }
    if (lnfact && fact != 'N') {
        PyErr_SetString(PyExc_ValueError, "Parameter 3 (FACT) must be 'N', 'C', or 'U'");
        return NULL;
    }
    if (uplo != 'U' && uplo != 'L') {
        PyErr_SetString(PyExc_ValueError, "Parameter 4 (UPLO) must be 'U' or 'L'");
        return NULL;
    }
    if (n < 0) {
        PyErr_SetString(PyExc_ValueError, "n must be non-negative");
        return NULL;
    }
    if (m < 0) {
        PyErr_SetString(PyExc_ValueError, "m must be non-negative");
        return NULL;
    }

    i32 lda = 1, ldb = 1, ldq = 1, ldr = 1, ldl = 1, ldg = 1;
    f64 *a_data = NULL, *b_data = NULL, *q_data = NULL;
    f64 *r_data = NULL, *l_data = NULL, *g_data = NULL;

    if (ljobl && a_obj != Py_None) {
        a_array = (PyArrayObject*)PyArray_FROM_OTF(a_obj, NPY_DOUBLE,
                                                   NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
        if (!a_array) goto cleanup;
        lda = (i32)PyArray_DIM(a_array, 0);
        a_data = (f64*)PyArray_DATA(a_array);
    }

    b_array = (PyArrayObject*)PyArray_FROM_OTF(b_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (!b_array) goto cleanup;
    ldb = (i32)PyArray_DIM(b_array, 0);
    b_data = (f64*)PyArray_DATA(b_array);

    if (ljobl && q_obj != Py_None) {
        q_array = (PyArrayObject*)PyArray_FROM_OTF(q_obj, NPY_DOUBLE,
                                                   NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
        if (!q_array) goto cleanup;
        ldq = (i32)PyArray_DIM(q_array, 0);
        q_data = (f64*)PyArray_DATA(q_array);
    }

    r_array = (PyArrayObject*)PyArray_FROM_OTF(r_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (!r_array) goto cleanup;
    ldr = (i32)PyArray_DIM(r_array, 0);
    r_data = (f64*)PyArray_DATA(r_array);

    if (ljobl && l_obj != Py_None) {
        l_array = (PyArrayObject*)PyArray_FROM_OTF(l_obj, NPY_DOUBLE,
                                                   NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
        if (!l_array) goto cleanup;
        ldl = (i32)PyArray_DIM(l_array, 0);
        l_data = (f64*)PyArray_DATA(l_array);
    }

    if (ljobg && g_obj != Py_None) {
        g_array = (PyArrayObject*)PyArray_FROM_OTF(g_obj, NPY_DOUBLE,
                                                   NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
        if (!g_array) goto cleanup;
        ldg = (i32)PyArray_DIM(g_array, 0);
        g_data = (f64*)PyArray_DATA(g_array);
    }

    i32 ldwork;
    if (lfactc) {
        ldwork = 1;
    } else if (lfactu) {
        ldwork = (ljobg || ljobl) ? (n * m > 1 ? n * m : 1) : 1;
    } else {
        if (ljobg || ljobl) {
            i32 nm = n * m;
            i32 tmp = 3 * m > nm ? 3 * m : nm;
            ldwork = tmp > 2 ? tmp : 2;
        } else {
            ldwork = 3 * m > 2 ? 3 * m : 2;
        }
    }

    f64 *dwork = (f64*)malloc(ldwork * sizeof(f64));
    i32 *iwork = (i32*)malloc((m > 0 ? m : 1) * sizeof(i32));
    i32 *ipiv = (i32*)malloc((m > 0 ? m : 1) * sizeof(i32));

    if (!dwork || !iwork || !ipiv) {
        free(dwork); free(iwork); free(ipiv);
        PyErr_NoMemory();
        goto cleanup;
    }

    sb02mt(jobg_str, jobl_str, fact_str, uplo_str,
           n, m, a_data, lda, b_data, ldb, q_data, ldq,
           r_data, ldr, l_data, ldl, ipiv, &oufact, g_data, ldg,
           iwork, dwork, ldwork, &info);

    free(dwork);
    free(iwork);
    free(ipiv);

    if (a_array) PyArray_ResolveWritebackIfCopy(a_array);
    if (b_array) PyArray_ResolveWritebackIfCopy(b_array);
    if (q_array) PyArray_ResolveWritebackIfCopy(q_array);
    if (r_array) PyArray_ResolveWritebackIfCopy(r_array);
    if (l_array) PyArray_ResolveWritebackIfCopy(l_array);
    if (g_array) PyArray_ResolveWritebackIfCopy(g_array);

    if (info < 0) {
        Py_XDECREF(a_array);
        Py_XDECREF(b_array);
        Py_XDECREF(q_array);
        Py_XDECREF(r_array);
        Py_XDECREF(l_array);
        Py_XDECREF(g_array);
        PyErr_Format(PyExc_ValueError, "Parameter %d had an illegal value", -info);
        return NULL;
    }

    PyObject *result;
    if (ljobg && !ljobl) {
        result = Py_BuildValue("Oii", g_array, oufact, info);
    } else if (!ljobg && ljobl) {
        result = Py_BuildValue("OOOOii", a_array, b_array, q_array, l_array, oufact, info);
    } else if (ljobg && ljobl) {
        result = Py_BuildValue("OOOOOii", a_array, b_array, q_array, l_array, g_array, oufact, info);
    } else {
        result = Py_BuildValue("ii", oufact, info);
    }

    Py_XDECREF(a_array);
    Py_XDECREF(b_array);
    Py_XDECREF(q_array);
    Py_XDECREF(r_array);
    Py_XDECREF(l_array);
    Py_XDECREF(g_array);

    return result;

cleanup:
    Py_XDECREF(a_array);
    Py_XDECREF(b_array);
    Py_XDECREF(q_array);
    Py_XDECREF(r_array);
    Py_XDECREF(l_array);
    Py_XDECREF(g_array);
    return NULL;
}

/* Python wrapper for sb02nd */
static PyObject* py_sb02nd(PyObject* self, PyObject* args) {
    const char *dico_str, *fact_str, *uplo_str, *jobl_str;
    i32 n, m, p;
    f64 rnorm;
    PyObject *a_obj, *b_obj, *r_obj, *ipiv_obj, *l_obj, *x_obj;
    PyArrayObject *a_array = NULL, *b_array = NULL, *r_array = NULL;
    PyArrayObject *ipiv_array = NULL, *l_array = NULL, *x_array = NULL;
    PyArrayObject *f_array = NULL;
    i32 oufact[2] = {0, 0};
    i32 info;

    if (!PyArg_ParseTuple(args, "ssssiiiOOOOOOd",
                          &dico_str, &fact_str, &uplo_str, &jobl_str,
                          &n, &m, &p, &a_obj, &b_obj, &r_obj,
                          &ipiv_obj, &l_obj, &x_obj, &rnorm)) {
        return NULL;
    }

    char dico = toupper((unsigned char)dico_str[0]);
    char fact = toupper((unsigned char)fact_str[0]);
    char uplo = toupper((unsigned char)uplo_str[0]);
    char jobl = toupper((unsigned char)jobl_str[0]);

    bool discr = (dico == 'D');
    bool lfactc = (fact == 'C');
    bool lfactd = (fact == 'D');
    bool lfactu = (fact == 'U');
    bool withl = (jobl == 'N');
    bool lfacta = lfactc || lfactd || lfactu;
    bool lnfact = !lfacta;

    if (!discr && dico != 'C') {
        PyErr_SetString(PyExc_ValueError, "DICO must be 'D' or 'C'");
        return NULL;
    }
    if ((lnfact && fact != 'N') || (discr && lfactu)) {
        PyErr_SetString(PyExc_ValueError, "FACT must be 'N', 'D', 'C', or 'U' (U not for discrete)");
        return NULL;
    }
    if (uplo != 'U' && uplo != 'L') {
        PyErr_SetString(PyExc_ValueError, "UPLO must be 'U' or 'L'");
        return NULL;
    }
    if (!withl && jobl != 'Z') {
        PyErr_SetString(PyExc_ValueError, "JOBL must be 'Z' or 'N'");
        return NULL;
    }
    if (n < 0) {
        PyErr_SetString(PyExc_ValueError, "n must be non-negative");
        return NULL;
    }
    if (m < 0) {
        PyErr_SetString(PyExc_ValueError, "m must be non-negative");
        return NULL;
    }

    i32 lda = 1, ldb = 1, ldr = 1, ldl = 1, ldx = 1, ldf = 1;
    f64 *a_data = NULL, *b_data = NULL, *r_data = NULL;
    f64 *l_data = NULL, *x_data = NULL, *f_data = NULL;
    i32 *ipiv_data = NULL;

    if (discr && a_obj != Py_None) {
        a_array = (PyArrayObject*)PyArray_FROM_OTF(a_obj, NPY_DOUBLE,
                                                   NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
        if (!a_array) goto cleanup;
        lda = (i32)PyArray_DIM(a_array, 0);
        a_data = (f64*)PyArray_DATA(a_array);
    }

    b_array = (PyArrayObject*)PyArray_FROM_OTF(b_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (!b_array) goto cleanup;
    ldb = (i32)PyArray_DIM(b_array, 0);
    b_data = (f64*)PyArray_DATA(b_array);

    r_array = (PyArrayObject*)PyArray_FROM_OTF(r_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (!r_array) goto cleanup;
    ldr = (i32)PyArray_DIM(r_array, 0);
    r_data = (f64*)PyArray_DATA(r_array);

    ipiv_array = (PyArrayObject*)PyArray_FROM_OTF(ipiv_obj, NPY_INT32,
                                                  NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (!ipiv_array) goto cleanup;
    ipiv_data = (i32*)PyArray_DATA(ipiv_array);

    if (withl && l_obj != Py_None) {
        l_array = (PyArrayObject*)PyArray_FROM_OTF(l_obj, NPY_DOUBLE,
                                                   NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
        if (!l_array) goto cleanup;
        ldl = (i32)PyArray_DIM(l_array, 0);
        l_data = (f64*)PyArray_DATA(l_array);
    }

    x_array = (PyArrayObject*)PyArray_FROM_OTF(x_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (!x_array) goto cleanup;
    ldx = (i32)PyArray_DIM(x_array, 0);
    x_data = (f64*)PyArray_DATA(x_array);

    ldf = m > 1 ? m : 1;
    npy_intp f_dims[2] = {ldf, n};
    npy_intp f_strides[2] = {sizeof(f64), ldf * sizeof(f64)};
    f_data = (f64*)malloc(ldf * n * sizeof(f64));
    if (!f_data) {
        PyErr_NoMemory();
        goto cleanup;
    }
    f_array = (PyArrayObject*)PyArray_New(&PyArray_Type, 2, f_dims, NPY_DOUBLE,
                                          f_strides, f_data, 0, NPY_ARRAY_FARRAY, NULL);
    if (!f_array) {
        free(f_data);
        goto cleanup;
    }
    PyArray_ENABLEFLAGS(f_array, NPY_ARRAY_OWNDATA);

    i32 ldwork;
    if (discr) {
        if (lnfact) {
            i32 tmp = 3 * m > n ? 3 * m : n;
            ldwork = tmp > 2 ? tmp : 2;
        } else {
            i32 tmp1 = n + 3 * m + 2;
            i32 tmp2 = 4 * n + 1;
            ldwork = tmp1 > tmp2 ? tmp1 : tmp2;
        }
    } else {
        if (lfactu) {
            ldwork = 2 * m > 2 ? 2 * m : 2;
        } else {
            ldwork = 3 * m > 2 ? 3 * m : 2;
        }
    }
    i32 nm = n * m;
    ldwork = ldwork > nm ? ldwork : nm;

    f64 *dwork = (f64*)malloc(ldwork * sizeof(f64));
    if (!dwork) {
        PyErr_NoMemory();
        goto cleanup;
    }

    sb02nd(dico_str, fact_str, uplo_str, jobl_str,
           n, m, p, a_data, lda, b_data, ldb, r_data, ldr,
           ipiv_data, l_data, ldl, x_data, ldx, rnorm,
           f_data, ldf, oufact, dwork, ldwork, &info);

    f64 rcond = dwork[1];
    free(dwork);

    if (a_array) PyArray_ResolveWritebackIfCopy(a_array);
    if (b_array) PyArray_ResolveWritebackIfCopy(b_array);
    if (r_array) PyArray_ResolveWritebackIfCopy(r_array);
    if (ipiv_array) PyArray_ResolveWritebackIfCopy(ipiv_array);
    if (l_array) PyArray_ResolveWritebackIfCopy(l_array);
    if (x_array) PyArray_ResolveWritebackIfCopy(x_array);

    if (info < 0) {
        Py_XDECREF(a_array);
        Py_XDECREF(b_array);
        Py_XDECREF(r_array);
        Py_XDECREF(ipiv_array);
        Py_XDECREF(l_array);
        Py_XDECREF(x_array);
        Py_XDECREF(f_array);
        PyErr_Format(PyExc_ValueError, "Parameter %d had an illegal value", -info);
        return NULL;
    }

    npy_intp oufact_dims[1] = {2};
    PyArrayObject *oufact_array = (PyArrayObject*)PyArray_SimpleNew(1, oufact_dims, NPY_INT32);
    if (!oufact_array) {
        Py_XDECREF(a_array);
        Py_XDECREF(b_array);
        Py_XDECREF(r_array);
        Py_XDECREF(ipiv_array);
        Py_XDECREF(l_array);
        Py_XDECREF(x_array);
        Py_XDECREF(f_array);
        return NULL;
    }
    ((i32*)PyArray_DATA(oufact_array))[0] = oufact[0];
    ((i32*)PyArray_DATA(oufact_array))[1] = oufact[1];

    PyObject *result = Py_BuildValue("OOOOdi", f_array, r_array, x_array,
                                     oufact_array, rcond, info);

    Py_XDECREF(a_array);
    Py_XDECREF(b_array);
    Py_XDECREF(r_array);
    Py_XDECREF(ipiv_array);
    Py_XDECREF(l_array);
    Py_XDECREF(x_array);
    Py_XDECREF(f_array);
    Py_XDECREF(oufact_array);

    return result;

cleanup:
    Py_XDECREF(a_array);
    Py_XDECREF(b_array);
    Py_XDECREF(r_array);
    Py_XDECREF(ipiv_array);
    Py_XDECREF(l_array);
    Py_XDECREF(x_array);
    Py_XDECREF(f_array);
    return NULL;
}

/* Python wrapper for ib01oy */
static PyObject* py_ib01oy(PyObject* self, PyObject* args) {
    i32 ns, nmax, n;
    PyObject *sv_obj;
    PyArrayObject *sv_array;
    i32 info;

    if (!PyArg_ParseTuple(args, "iiiO", &ns, &nmax, &n, &sv_obj)) {
        return NULL;
    }

    /* Validate parameters before array conversion */
    if (ns <= 0) {
        PyErr_SetString(PyExc_ValueError, "NS must be positive");
        return NULL;
    }

    if (nmax < 0 || nmax > ns) {
        PyErr_SetString(PyExc_ValueError, "NMAX must be in range [0, NS]");
        return NULL;
    }

    if (n < 0 || n > ns) {
        PyErr_SetString(PyExc_ValueError, "N must be in range [0, NS]");
        return NULL;
    }

    /* Convert SV array */
    sv_array = (PyArrayObject*)PyArray_FROM_OTF(sv_obj, NPY_DOUBLE,
                                                NPY_ARRAY_IN_FARRAY);
    if (sv_array == NULL) {
        return NULL;
    }

    /* Validate SV size */
    npy_intp sv_size = PyArray_SIZE(sv_array);
    if (sv_size < ns) {
        Py_DECREF(sv_array);
        PyErr_SetString(PyExc_ValueError, "SV must have at least NS elements");
        return NULL;
    }

    f64 *sv = (f64*)PyArray_DATA(sv_array);

    /* Call C function */
    SLC_IB01OY(ns, nmax, &n, sv, &info);

    Py_DECREF(sv_array);

    if (info < 0) {
        PyErr_Format(PyExc_ValueError, "Parameter %d had an illegal value", -info);
        return NULL;
    }

    return Py_BuildValue("ii", n, info);
}

/* Python wrapper for mb01sd */
static PyObject* py_mb01sd(PyObject* self, PyObject* args) {
    const char *jobs_str;
    PyObject *a_obj, *r_obj, *c_obj;
    PyArrayObject *a_array, *r_array, *c_array;

    if (!PyArg_ParseTuple(args, "sOOO", &jobs_str, &a_obj, &r_obj, &c_obj)) {
        return NULL;
    }

    char jobs = jobs_str[0];

    a_array = (PyArrayObject*)PyArray_FROM_OTF(a_obj, NPY_DOUBLE,
                                               NPY_ARRAY_FARRAY | NPY_ARRAY_WRITEBACKIFCOPY);
    if (a_array == NULL) {
        return NULL;
    }

    r_array = (PyArrayObject*)PyArray_FROM_OTF(r_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (r_array == NULL) {
        Py_DECREF(a_array);
        return NULL;
    }

    c_array = (PyArrayObject*)PyArray_FROM_OTF(c_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (c_array == NULL) {
        Py_DECREF(a_array);
        Py_DECREF(r_array);
        return NULL;
    }

    npy_intp *a_dims = PyArray_DIMS(a_array);
    i32 m = (i32)a_dims[0];
    i32 n = (i32)a_dims[1];
    i32 lda = m;

    f64 *a_data = (f64*)PyArray_DATA(a_array);
    const f64 *r_data = (const f64*)PyArray_DATA(r_array);
    const f64 *c_data = (const f64*)PyArray_DATA(c_array);

    mb01sd(jobs, m, n, a_data, lda, r_data, c_data);

    Py_DECREF(r_array);
    Py_DECREF(c_array);

    PyArray_ResolveWritebackIfCopy(a_array);

    PyObject *result = (PyObject*)a_array;
    return result;
}

/* Module method definitions */
static PyMethodDef SlicotMethods[] = {
    {"mb04od", py_mb04od, METH_VARARGS,
     "QR factorization of structured block matrix.\n\n"
     "Computes QR factorization of first block column and applies\n"
     "transformations to second block column:\n"
     "  Q' * [R  B] = [R_  B_]\n"
     "       [A  C]   [0   C_]\n\n"
     "Parameters:\n"
     "  uplo (str): 'U' for A upper trapezoidal, 'F' for A full\n"
     "  n (int): Order of R\n"
     "  m (int): Columns in B, C\n"
     "  p (int): Rows in A, C\n"
     "  r (ndarray): n-by-n upper triangular matrix R (F-order), modified\n"
     "  a (ndarray): p-by-n matrix A (F-order), overwritten with Householder vectors\n"
     "  b (ndarray): n-by-m matrix B (F-order), modified\n"
     "  c (ndarray): p-by-m matrix C (F-order), modified\n\n"
     "Returns:\n"
     "  (r, a, b, c, tau): Updated matrices and Householder scalars\n"},

    {"mb04ow", (PyCFunction)py_mb04ow, METH_VARARGS | METH_KEYWORDS,
     "Performs a QR factorization update.\n\n"
     "Parameters:\n"
     "  m (int): Rows of U1\n"
     "  n (int): Rows of T\n"
     "  p (int): Columns of B, C\n"
     "  a (ndarray): Matrix U (lda x m+n), modified in place\n"
     "  t (ndarray): Matrix T (ldt x n), modified in place\n"
     "  x (ndarray): Vector x (m+n), modified in place\n"
     "  b (ndarray): Matrix B (ldb x p), modified in place\n"
     "  c (ndarray): Matrix C (ldc x p), modified in place\n"
     "  d (ndarray): Vector d (p), modified in place\n"
     "  incx (int, optional): Increment for x\n"
     "  incd (int, optional): Increment for d\n\n"
     "Returns:\n"
     "  (a, t, x, b, c, d): Modified arrays\n"},

    {"nf01ay", (PyCFunction)py_nf01ay, METH_VARARGS,
     "Calculate the output of a set of neural networks.\n\n"
     "Parameters:\n"
     "  nsmp (int): Number of training samples\n"
     "  nz (int): Length of each input sample\n"
     "  l (int): Length of each output sample\n"
     "  ipar (ndarray): Integer parameters (nn, ...)\n"
     "  wb (ndarray): Weights and biases\n"
     "  z (ndarray): Input samples (nsmp x nz, F-order)\n\n"
     "Returns:\n"
     "  (y, info): Output samples, exit code\n"},

    {"nf01by", (PyCFunction)py_nf01by, METH_VARARGS,
     "Compute the Jacobian of the error function for a neural network.\n\n"
     "Parameters:\n"
     "  cjte (str): 'C' to compute J'*e, 'N' to skip\n"
     "  nsmp (int): Number of training samples\n"
     "  nz (int): Length of each input sample\n"
     "  l (int): Length of each output sample (must be 1)\n"
     "  ipar (ndarray): Integer parameters (nn, ...)\n"
     "  wb (ndarray): Weights and biases\n"
     "  z (ndarray): Input samples (nsmp x nz, F-order)\n"
     "  e (ndarray): Error vector (nsmp)\n\n"
     "Returns:\n"
     "  (j, jte, info): Jacobian, J'*e, exit code\n"},

    {"nf01bs", (PyCFunction)py_nf01bs, METH_VARARGS,
     "QR factorization of Jacobian in compressed form.\n\n"
     "Parameters:\n"
     "  n (int): Number of columns of J\n"
     "  ipar (ndarray): Integer parameters (st, bn, bsm, bsn)\n"
     "  fnorm (float): Norm of error vector\n"
     "  j (ndarray): Jacobian matrix (ldj x nc, F-order)\n"
     "  e (ndarray): Error vector\n\n"
     "Returns:\n"
     "  (r, e, jnorms, gnorm, ipvt, info): R factor, Q'*e, norms, gradient norm, permutation, exit code\n"},

    {"nf01br", (PyCFunction)py_nf01br, METH_VARARGS,
     "Solve system of linear equations R*x = b or R'*x = b in least squares sense.\n\n"
     "Parameters:\n"
     "  cond (str): 'E' = estimate condition, 'N' = check zeros, 'U' = use rank\n"
     "  uplo (str): 'U' or 'L' - storage scheme\n"
     "  trans (str): 'N'/'T'/'C' - transpose option\n"
     "  n (int): Order of matrix R\n"
     "  ipar (ndarray): Integer parameters (st, bn, bsm, bsn)\n"
     "  r (ndarray): Matrix R (ldr x nc, F-order)\n"
     "  sdiag (ndarray): Diagonal elements (if uplo='L')\n"
     "  s (ndarray): Transpose of last block column (if uplo='L')\n"
     "  b (ndarray): RHS vector b\n"
     "  ranks (ndarray): Ranks\n"
     "  tol (float): Tolerance for rank determination\n\n"
     "Returns:\n"
     "  (b, ranks, info): Solution x, ranks, exit code\n"},

    {"tb01vy", (PyCFunction)py_tb01vy, METH_VARARGS | METH_KEYWORDS,
     "Convert output normal form to state-space representation.\n\n"
     "Parameters:\n"
     "  n (int): System order\n"
     "  m (int): Number of inputs\n"
     "  l (int): Number of outputs\n"
     "  theta (ndarray): Parameter vector (N*(L+M+1)+L*M, F-order)\n"
     "  apply (str, optional): 'A' = apply bijective mapping, 'N' = no mapping (default)\n\n"
     "Returns:\n"
     "  (a, b, c, d, x0, info): State-space matrices, initial state, exit code\n"},

    {"tb01wd", py_tb01wd, METH_VARARGS,
     "Reduce state matrix to real Schur form via orthogonal transformation.\n\n"
     "Parameters:\n"
     "  n (int): Order of state matrix A\n"
     "  m (int): Number of inputs (columns of B)\n"
     "  p (int): Number of outputs (rows of C)\n"
     "  a (ndarray): State dynamics matrix (n x n, F-order)\n"
     "  b (ndarray): Input matrix (n x m, F-order)\n"
     "  c (ndarray): Output matrix (p x n, F-order)\n\n"
     "Returns:\n"
     "  (a, b, c, u, wr, wi, info): Transformed system, Schur vectors, eigenvalues, exit code\n"},

    {"tf01mx", py_tf01mx, METH_VARARGS,
     "Output sequence of linear time-invariant open-loop system.\n\n"
     "Parameters:\n"
     "  n (int): Order of matrix A\n"
     "  m (int): Number of system inputs\n"
     "  p (int): Number of system outputs\n"
     "  ny (int): Number of output vectors to compute\n"
     "  s (ndarray): System matrix [A B; C D] (n+p x n+m, F-order)\n"
     "  u (ndarray): Input sequence (ny x m, F-order)\n"
     "  x (ndarray): Initial state vector (n, F-order)\n\n"
     "Returns:\n"
     "  (y, x, info): Output sequence, final state, exit code\n"},

    {"ma02ad", py_ma02ad, METH_VARARGS,
     "Transpose all or part of a matrix.\n\n"
     "Parameters:\n"
     "  job (str): 'U' for upper triangular, 'L' for lower, else full matrix\n"
     "  a (ndarray): Input matrix (m x n, F-order)\n\n"
     "Returns:\n"
     "  b (ndarray): Transposed matrix (n x m, F-order)\n"},

    {"ma02ed", py_ma02ed, METH_VARARGS,
     "Store by symmetry the upper or lower triangle of a symmetric matrix.\n\n"
     "Parameters:\n"
     "  uplo (str): 'U' if upper triangle given, 'L' if lower triangle given\n"
     "  a (ndarray): Symmetric matrix (n x n, F-order), modified in place\n\n"
     "Returns:\n"
     "  a (ndarray): Completed symmetric matrix with both triangles\n"},

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

    {"mb03ud", (PyCFunction)py_mb03ud, METH_VARARGS | METH_KEYWORDS,
     "Singular value decomposition of upper triangular matrix.\n\n"
     "Computes SVD: A = Q*S*P' where Q, P are orthogonal and S is diagonal\n"
     "with non-negative singular values in descending order.\n\n"
     "Parameters:\n"
     "  n (int): Order of matrix A\n"
     "  a (ndarray): Upper triangular matrix A (n x n, F-order)\n"
     "               If jobp='V', returns P' on exit\n"
     "  jobq (str, optional): 'V' to compute Q, 'N' otherwise (default 'N')\n"
     "  jobp (str, optional): 'V' to compute P', 'N' otherwise (default 'N')\n"
     "  ldwork (int, optional): Workspace size (-1 for query, default auto)\n\n"
     "Returns:\n"
     "  (sv, p, q, info): Singular values, right vectors P', left vectors Q, exit code\n"
     "  - sv: array of singular values (descending order)\n"
     "  - p: P' matrix if jobp='V', else None\n"
     "  - q: Q matrix if jobq='V', else None\n"
     "  - info: 0=success, <0=invalid param, >0=convergence failure\n"},

    {"mb04id", (PyCFunction)py_mb04id, METH_VARARGS | METH_KEYWORDS,
     "QR factorization of matrix with lower-left zero triangle.\n\n"
     "Computes A = Q*R where A has p-by-min(p,m) zero triangle in lower-left.\n"
     "Optionally applies Q' to matrix B.\n\n"
     "Parameters:\n"
     "  n (int): Number of rows of A\n"
     "  m (int): Number of columns of A\n"
     "  p (int): Order of zero triangle\n"
     "  a (ndarray): Matrix A (n x m, F-order), modified in place\n"
     "  b (ndarray, optional): Matrix B (n x l, F-order), modified in place\n"
     "  l (int, optional): Number of columns of B (default 0)\n"
     "  ldwork (int, optional): Workspace size (-1 for query, default auto)\n\n"
     "Returns:\n"
     "  If l>0: (b, tau, info) - transformed B, Householder factors, exit code\n"
     "  If l=0: (tau, info) - Householder factors, exit code\n"
     "  If ldwork=-1: adds optimal workspace size to tuple\n"},

    {"mb04iy", (PyCFunction)py_mb04iy, METH_VARARGS | METH_KEYWORDS,
     "Apply orthogonal transformations from MB04ID to matrix C.\n\n"
     "Applies Q or Q' to matrix C, where Q is product of elementary reflectors\n"
     "as returned by MB04ID (special structure for lower-left zero triangle).\n\n"
     "Parameters:\n"
     "  side (str): 'L' for left (Q*C or Q'*C), 'R' for right (C*Q or C*Q')\n"
     "  trans (str): 'N' for Q, 'T' for Q'\n"
     "  a (ndarray): Reflector storage (lda x k, F-order), modified but restored\n"
     "  tau (ndarray): Reflector scalar factors (k,)\n"
     "  c (ndarray): Matrix C (n x m, F-order), modified in place\n"
     "  p (int, optional): Order of zero triangle (default 0)\n\n"
     "Returns:\n"
     "  (c, info): Transformed matrix and exit code\n"},

    {"mb04kd", py_mb04kd, METH_VARARGS,
     "QR factorization of special structured block matrix.\n\n"
     "Computes Q' * [[R],[A B]] = [[R_bar C],[0 D]]\n\n"
     "Parameters:\n"
     "  uplo (str): 'U' if A is upper trapezoidal, 'F' if A is full\n"
     "  n (int): Order of R and R_bar\n"
     "  m (int): Number of columns of B, C, D\n"
     "  p (int): Number of rows of A, B, D\n"
     "  r (ndarray): Upper triangular matrix R (n x n, F-order)\n"
     "  a (ndarray): Matrix A (p x n, F-order)\n"
     "  b (ndarray): Matrix B (p x m, F-order)\n\n"
     "Returns:\n"
     "  (r_bar, a_out, d, c, tau): Transformed matrices and Householder factors\n"},

    {"mb04oy", py_mb04oy, METH_VARARGS,
     "Apply Householder reflector to matrix [A; B].\n\n"
     "Applies H = I - tau*[1;v]*[1;v]' to (m+1)-by-n matrix [A; B],\n"
     "where A has one row.\n\n"
     "Parameters:\n"
     "  m (int): Number of rows of B\n"
     "  n (int): Number of columns\n"
     "  v (ndarray): Householder vector (m, F-order)\n"
     "  tau (float): Householder scalar\n"
     "  a (ndarray): Matrix A (1 x n, F-order, modified in place)\n"
     "  b (ndarray): Matrix B (m x n, F-order, modified in place)\n\n"
     "Returns:\n"
     "  (a, b): Updated matrices\n"},

    {"mb04ny", py_mb04ny, METH_VARARGS,
     "Apply Householder reflector to matrix [A B] from the right.\n\n"
     "Applies H = I - tau*[1;v]*[1;v]' to m-by-(n+1) matrix [A B],\n"
     "where A has one column.\n\n"
     "Parameters:\n"
     "  m (int): Number of rows of A and B\n"
     "  n (int): Number of columns of B\n"
     "  v (ndarray): Householder vector (1+(n-1)*abs(incv), F-order)\n"
     "  incv (int): Increment between elements of v\n"
     "  tau (float): Householder scalar\n"
     "  a (ndarray): Matrix A (m x 1, F-order, modified in place)\n"
     "  b (ndarray): Matrix B (m x n, F-order, modified in place)\n\n"
     "Returns:\n"
     "  (a, b): Updated matrices\n"},

    {"mb01rb", py_mb01rb, METH_VARARGS,
     "Block triangular symmetric rank-k update (BLAS 3 version).\n\n"
     "Computes R = alpha*R + beta*op(A)*B or R = alpha*R + beta*B*op(A)\n"
     "where only the specified triangle is computed using block algorithm.\n\n"
     "Parameters:\n"
     "  side (str): 'L' for left (R = alpha*R + beta*op(A)*B)\n"
     "              'R' for right (R = alpha*R + beta*B*op(A))\n"
     "  uplo (str): 'U' for upper triangle, 'L' for lower triangle\n"
     "  trans (str): 'N' for op(A)=A, 'T'/'C' for op(A)=A'\n"
     "  m (int): Order of matrix R\n"
     "  n (int): Dimension for product\n"
     "  alpha (float): Scalar multiplier for R\n"
     "  beta (float): Scalar multiplier for product\n"
     "  r (ndarray): m-by-m matrix R (F-order), modified in place\n"
     "  a (ndarray): Matrix A (F-order)\n"
     "  b (ndarray): Matrix B (F-order)\n\n"
     "Returns:\n"
     "  (r, info): Updated R (triangle only) and exit code\n"},

    {"mb01td", py_mb01td, METH_VARARGS,
     "Product of upper quasi-triangular matrices B := A * B.\n\n"
     "Computes A * B where A and B are upper quasi-triangular matrices\n"
     "(block upper triangular with 1x1 or 2x2 diagonal blocks) with the\n"
     "same structure. Result is returned in B.\n\n"
     "Parameters:\n"
     "  a (ndarray): N-by-N upper quasi-triangular matrix A (F-order)\n"
     "  b (ndarray): N-by-N upper quasi-triangular matrix B (F-order)\n\n"
     "Returns:\n"
     "  (b, info): Product A*B in B, exit code (0=success, 1=structure mismatch)\n"},

    {"mb01rx", py_mb01rx, METH_VARARGS,
     "Triangular symmetric rank-k update.\n\n"
     "Computes R = alpha*R + beta*op(A)*B or R = alpha*R + beta*B*op(A)\n"
     "where only the specified triangle is computed.\n\n"
     "Parameters:\n"
     "  side (str): 'L' for left (R = alpha*R + beta*op(A)*B)\n"
     "              'R' for right (R = alpha*R + beta*B*op(A))\n"
     "  uplo (str): 'U' for upper triangle, 'L' for lower triangle\n"
     "  trans (str): 'N' for op(A)=A, 'T'/'C' for op(A)=A'\n"
     "  m (int): Order of matrix R\n"
     "  n (int): Dimension for product\n"
     "  alpha (float): Scalar multiplier for R\n"
     "  beta (float): Scalar multiplier for product\n"
     "  r (ndarray): m-by-m matrix R (F-order), modified in place\n"
     "  a (ndarray): Matrix A (F-order)\n"
     "  b (ndarray): Matrix B (F-order)\n\n"
     "Returns:\n"
     "  (r, info): Updated R (triangle only) and exit code\n"},

    {"mb03od", (PyCFunction)py_mb03od, METH_VARARGS | METH_KEYWORDS,
     "Incremental rank estimation for QR factorization.\n\n"
     "Parameters:\n"
     "  m (int): Number of rows\n"
     "  n (int): Number of columns\n"
     "  a (ndarray): Matrix array (column-major, shape (m,n))\n"
     "  rcond (float): Rank threshold\n"
     "  svlmax (float): Parent matrix singular value estimate\n"
     "  jobqr (str, optional): 'Q' = perform QR (default), 'N' = use existing\n\n"
     "Returns:\n"
     "  (jpvt, rank, sval, info): Pivot indices, rank, singular values, exit code\n"},

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

    {"mb01uy", py_mb01uy, METH_VARARGS,
     "Compute matrix product T := alpha*op(T)*A or T := alpha*A*op(T).\n\n"
     "Parameters:\n"
     "  side (str): 'L' or 'R' - triangular matrix position\n"
     "  uplo (str): 'U' or 'L' - upper/lower triangular\n"
     "  trans (str): 'N'/'T'/'C' - transpose option\n"
     "  m (int): Number of rows of A\n"
     "  n (int): Number of columns of A\n"
     "  alpha (float): Scalar multiplier\n"
     "  t (ndarray): Triangular matrix (F-order), overwritten with result\n"
     "  a (ndarray): Matrix A (F-order)\n\n"
     "Returns:\n"
     "  (t, info): Result matrix and exit code\n"},

    {"mb02yd", py_mb02yd, METH_VARARGS,
     "Solve augmented system A*x = b, D*x = 0 in least squares sense.\n\n"
     "Parameters:\n"
     "  cond (str): 'E' = estimate condition, 'N' = check zeros, 'U' = use rank\n"
     "  n (int): Order of matrix R\n"
     "  r (ndarray): Upper triangular matrix R (n x n, F-order)\n"
     "  ipvt (ndarray): Permutation vector (1-based indices)\n"
     "  diag (ndarray): Diagonal elements of D\n"
     "  qtb (ndarray): First n elements of Q'*b\n"
     "  rank (int): Input rank (COND='U') or 0 otherwise\n"
     "  tol (float): Tolerance for rank determination (COND='E')\n\n"
     "Returns:\n"
     "  (x, rank, info): Solution vector, estimated rank, exit code\n"},

    {"mb02ud", (PyCFunction)py_mb02ud, METH_VARARGS | METH_KEYWORDS,
     "Minimum norm least squares solution using SVD.\n\n"
     "Solves op(R)*X = alpha*B (side='L') or X*op(R) = alpha*B (side='R')\n"
     "where R is upper triangular, using singular value decomposition.\n\n"
     "Parameters:\n"
     "  fact (str): 'N' to compute SVD, 'F' if SVD already available\n"
     "  side (str): 'L' for left, 'R' for right\n"
     "  trans (str): 'N' for op(R)=R, 'T'/'C' for op(R)=R'\n"
     "  jobp (str): 'P' to compute/use pseudoinverse, 'N' otherwise\n"
     "  m (int): Number of rows of B\n"
     "  n (int): Number of columns of B\n"
     "  alpha (float): Scalar multiplier\n"
     "  rcond (float): Rank threshold (not used if fact='F')\n"
     "  r (ndarray): L-by-L upper triangular matrix R (F-order)\n"
     "  b (ndarray): M-by-N matrix B (F-order)\n"
     "  q (ndarray, optional): Q matrix (required if fact='F')\n"
     "  sv (ndarray, optional): Singular values (required if fact='F')\n"
     "  rank (int, optional): Rank (required if fact='F')\n"
     "  rp (ndarray, optional): Pseudoinverse (if fact='F' and jobp='P')\n"
     "  ldwork (int, optional): Workspace size\n\n"
     "Returns:\n"
     "  (x, q, sv, rank, rp, info): Solution, Q matrix, singular values, rank, pseudoinverse, exit code\n"},

    {"md03by", (PyCFunction)py_md03by, METH_VARARGS | METH_KEYWORDS,
     "Compute Levenberg-Marquardt parameter for trust region subproblem.\n\n"
     "Parameters:\n"
     "  cond (str): 'E' = estimate condition, 'N' = check zeros, 'U' = use rank\n"
     "  n (int): Order of matrix R\n"
     "  r (ndarray): Upper triangular matrix R (n x n, F-order)\n"
     "  ipvt (ndarray): Permutation vector (1-based indices)\n"
     "  diag (ndarray): Diagonal scaling D (all nonzero)\n"
     "  qtb (ndarray): First n elements of Q'*b\n"
     "  delta (float): Trust region radius (> 0)\n"
     "  par (float): Initial LM parameter estimate (>= 0)\n"
     "  rank (int): Input rank (COND='U') or 0 otherwise\n"
     "  tol (float): Tolerance for rank determination (COND='E')\n\n"
     "Returns:\n"
     "  (r, par, rank, x, rx, info): Modified R, LM parameter, rank, solution, residual, exit code\n"},

    {"md03bb", (PyCFunction)py_md03bb, METH_VARARGS | METH_KEYWORDS,
     "Compute Levenberg-Marquardt parameter for compressed Jacobian (Wrapper for MD03BY).\n\n"
     "Parameters:\n"
     "  cond (str): 'E' = estimate condition, 'N' = check zeros, 'U' = use rank\n"
     "  n (int): Order of matrix R\n"
     "  r (ndarray): Upper triangular matrix R (n x n, F-order)\n"
     "  ipvt (ndarray): Permutation vector (1-based indices)\n"
     "  diag (ndarray): Diagonal scaling D (all nonzero)\n"
     "  qtb (ndarray): First n elements of Q'*b\n"
     "  delta (float): Trust region radius (> 0)\n"
     "  par (float): Initial LM parameter estimate (>= 0)\n"
     "  rank (int): Input rank (COND='U') or 0 otherwise\n"
     "  tol (float): Tolerance for rank determination (COND='E')\n\n"
     "Returns:\n"
     "  (r, par, rank, x, rx, info): Modified R, LM parameter, rank, solution, residual, exit code\n"},

    {"md03ba", (PyCFunction)py_md03ba, METH_VARARGS,
     "QR factorization with column pivoting for Levenberg-Marquardt.\n\n"
     "Parameters:\n"
     "  n (int): Number of columns of J\n"
     "  ipar (ndarray): Integer parameters (M=ipar[0])\n"
     "  fnorm (float): Norm of error vector\n"
     "  j (ndarray): Jacobian matrix (M x N, F-order)\n"
     "  e (ndarray): Error vector (M)\n\n"
     "Returns:\n"
     "  (r, e, jnorms, gnorm, ipvt, info): R factor, Q'*e, norms, gradient norm, permutation, exit code\n"},

    {"sg03br", py_sg03br, METH_VARARGS,
     "Compute complex Givens rotation in real arithmetic.\n\n"
     "Parameters:\n"
     "  xr (float): Real part of X\n"
     "  xi (float): Imaginary part of X\n"
     "  yr (float): Real part of Y\n"
     "  yi (float): Imaginary part of Y\n\n"
     "Returns:\n"
     "  (c, sr, si, zr, zi): Givens rotation parameters and result\n"},

    {"sg03bw", py_sg03bw, METH_VARARGS,
     "Solve generalized Sylvester equation for small systems.\n\n"
     "Parameters:\n"
     "  trans (str): 'N' or 'T' - equation form\n"
     "  a (ndarray): Upper quasitriangular matrix (m x m, F-order)\n"
     "  e (ndarray): Upper triangular matrix (m x m, F-order)\n"
     "  c (ndarray): Matrix C (n x n, F-order)\n"
     "  d (ndarray): Matrix D (n x n, F-order)\n"
     "  x (ndarray): Right-hand side Y, overwritten with solution (m x n, F-order)\n\n"
     "Returns:\n"
     "  (x, scale, info): Solution matrix, scale factor, and exit code\n"},

    {"sg03bu", py_sg03bu, METH_VARARGS,
     "Solve generalized discrete-time Lyapunov equation for Cholesky factor.\n\n"
     "Parameters:\n"
     "  trans (str): 'N' for TRANS='N', 'T' for TRANS='T'\n"
     "  a (ndarray): Quasitriangular matrix A (n x n, F-order)\n"
     "  e (ndarray): Upper triangular matrix E (n x n, F-order)\n"
     "  b (ndarray): Upper triangular matrix B (n x n, F-order)\n\n"
     "Returns:\n"
     "  (u, scale, info): Cholesky factor U, scale factor, exit code\n"},

    {"sg03bv", py_sg03bv, METH_VARARGS,
     "Solve generalized continuous-time Lyapunov equation for Cholesky factor.\n\n"
     "Parameters:\n"
     "  trans (str): 'N' for TRANS='N', 'T' for TRANS='T'\n"
     "  a (ndarray): Quasitriangular matrix A (n x n, F-order)\n"
     "  e (ndarray): Upper triangular matrix E (n x n, F-order)\n"
     "  b (ndarray): Upper triangular matrix B (n x n, F-order)\n\n"
     "Returns:\n"
     "  (u, scale, info): Cholesky factor U, scale factor, exit code\n"},

    {"sg03bd", py_sg03bd, METH_VARARGS,
     "Solve generalized Lyapunov equation for Cholesky factor.\n\n"
     "Parameters:\n"
     "  dico (str): 'C' for continuous-time, 'D' for discrete-time\n"
     "  fact (str): 'N' to compute factorization, 'F' if factorization supplied\n"
     "  trans (str): 'N' for op(K)=K, 'T' for op(K)=K^T\n"
     "  n (int): Order of matrices A and E\n"
     "  m (int): Number of rows in op(B)\n"
     "  a (ndarray): Matrix A (n x n, F-order)\n"
     "  e (ndarray): Matrix E (n x n, F-order)\n"
     "  b (ndarray): Matrix B (size depends on trans, F-order)\n\n"
     "Returns:\n"
     "  (u, scale, alphar, alphai, beta, info): Cholesky factor, scale, eigenvalues, exit code\n"},

    {"sg03bx", py_sg03bx, METH_VARARGS,
     "Solve 2x2 generalized Lyapunov equation.\n\n"
     "Parameters:\n"
     "  dico (str): 'C' for continuous-time, 'D' for discrete-time\n"
     "  trans (str): 'N' for op(K)=K, 'T' for op(K)=K^T\n"
     "  a (ndarray): Matrix A (2 x 2, F-order)\n"
     "  e (ndarray): Upper triangular matrix E (2 x 2, F-order)\n"
     "  b (ndarray): Upper triangular matrix B (2 x 2, F-order)\n\n"
     "Returns:\n"
     "  (u, scale, m1, m2, info): Cholesky factor, scale, auxiliary matrices, exit code\n"},

    {"tg01fd", py_tg01fd, METH_VARARGS,
     "Orthogonal reduction of descriptor system to SVD-like form.\n\n"
     "Parameters:\n"
     "  compq (str): 'N', 'I', or 'U' - Q computation mode\n"
     "  compz (str): 'N', 'I', or 'U' - Z computation mode\n"
     "  joba (str): 'N', 'R', or 'T' - A22 reduction mode\n"
     "  l (int): Number of rows of A, B, E\n"
     "  n (int): Number of columns of A, E, C\n"
     "  m (int): Number of columns of B\n"
     "  p (int): Number of rows of C\n"
     "  a (ndarray): State dynamics matrix (l x n, F-order)\n"
     "  e (ndarray): Descriptor matrix (l x n, F-order)\n"
     "  b (ndarray): Input/state matrix (l x m, F-order)\n"
     "  c (ndarray): State/output matrix (p x n, F-order)\n"
     "  tol (float): Tolerance for rank determination\n\n"
     "Returns:\n"
     "  (a, e, b, c, q, z, ranke, rnka22, info): Transformed system and ranks\n"},

    {"md03bd", (PyCFunction)(void(*)(void))py_md03bd, METH_VARARGS | METH_KEYWORDS,
     "Levenberg-Marquardt nonlinear least squares optimizer.\n\n"
     "Parameters:\n"
     "  m (int): Number of functions\n"
     "  n (int): Number of variables\n"
     "  x (ndarray): Initial guess (n,) or random if not provided\n"
     "  fcn (callable): Function f(x) returning residual vector (m,)\n"
     "  jac (callable): Function f(x) returning Jacobian matrix (m, n)\n"
     "  itmax (int): Maximum iterations (default 100)\n"
     "  ftol (float): Function tolerance (default sqrt(eps))\n"
     "  xtol (float): Solution tolerance (default sqrt(eps))\n"
     "  gtol (float): Gradient tolerance (default eps)\n\n"
     "Returns:\n"
     "  (x, nfev, njev, fnorm, iwarn, info): Solution, evaluations, norm, status\n"},

    {"ib01nd", (PyCFunction)py_ib01nd, METH_VARARGS,
     "SVD system order via block Hankel.\n\n"
     "Computes SVD of triangular factor R from QR factorization of\n"
     "concatenated block Hankel matrices to determine system order.\n\n"
     "Parameters:\n"
     "  meth (str): 'M' for MOESP, 'N' for N4SID\n"
     "  jobd (str): 'M' or 'N' for MOESP BD computation mode\n"
     "  nobr (int): Number of block rows (nobr > 0)\n"
     "  m (int): Number of system inputs (m >= 0)\n"
     "  l (int): Number of system outputs (l > 0)\n"
     "  r (ndarray): Upper triangular R matrix, dimension 2*(m+l)*nobr x 2*(m+l)*nobr\n"
     "  tol (float): Tolerance for rank estimation (N4SID only)\n\n"
     "Returns:\n"
     "  (r, sv, rcond1, rcond2, iwarn, info): Processed R, singular values, rconds, status\n"},

    {"ib01od", (PyCFunction)py_ib01od, METH_VARARGS,
     "Estimate system order from Hankel singular values.\n\n"
     "Estimates system order based on singular values of triangular factor\n"
     "from QR factorization of concatenated block Hankel matrices.\n\n"
     "Parameters:\n"
     "  ctrl (str): 'C' for user confirmation, 'N' for no confirmation\n"
     "  nobr (int): Number of block rows (nobr > 0)\n"
     "  l (int): Number of system outputs (l > 0)\n"
     "  sv (ndarray): Singular values, dimension (l*nobr), descending order\n"
     "  tol (float): Tolerance (>=0: threshold, 0: default, <0: gap-based)\n\n"
     "Returns:\n"
     "  (n, iwarn, info): Estimated order, warning, exit code\n"},

    {"ib01oy", (PyCFunction)py_ib01oy, METH_VARARGS,
     "User's confirmation of the system order.\n\n"
     "Non-interactive version for library use. Validates parameters and ensures N <= NMAX.\n\n"
     "Parameters:\n"
     "  ns (int): Number of singular values (ns > 0)\n"
     "  nmax (int): Maximum value of system order (0 <= nmax <= ns)\n"
     "  n (int): Estimated system order (0 <= n <= ns)\n"
     "  sv (ndarray): Singular values, dimension (ns), descending order\n\n"
     "Returns:\n"
     "  (n, info): Validated order, exit code\n"},

    {"ib01qd", (PyCFunction)py_ib01qd, METH_VARARGS,
     "Estimate initial state and system matrices B, D.\n\n"
     "Given (A, C) and input/output trajectories, estimates B, D, and x(0)\n"
     "for: x(k+1) = A*x(k) + B*u(k), y(k) = C*x(k) + D*u(k).\n"
     "Matrix A must be in real Schur form.\n\n"
     "Parameters:\n"
     "  jobx0 (str): 'X' to compute x(0), 'N' if x(0) known to be zero\n"
     "  job (str): 'B' for B only (D=0), 'D' for B and D\n"
     "  n (int): System order (n >= 0)\n"
     "  m (int): Number of inputs (m >= 0)\n"
     "  l (int): Number of outputs (l > 0)\n"
     "  a (ndarray): N-by-N state matrix A in Schur form (F-order)\n"
     "  c (ndarray): L-by-N output matrix C (F-order)\n"
     "  u (ndarray): NSMP-by-M input data (F-order)\n"
     "  y (ndarray): NSMP-by-L output data (F-order)\n"
     "  tol (float): Tolerance for rank estimation (<= 0 uses machine eps)\n\n"
     "Returns:\n"
     "  (x0, b, d, rcond_w2, rcond_u, iwarn, info):\n"
     "  - x0: Estimated initial state (n,)\n"
     "  - b: Estimated input matrix B (n x m)\n"
     "  - d: Estimated direct transmission D (l x m)\n"
     "  - rcond_w2: Reciprocal condition of W2\n"
     "  - rcond_u: Reciprocal condition of U (if JOB='D')\n"
     "  - iwarn: Warning (4 = rank-deficient)\n"
     "  - info: Exit code (0=success, 2=SVD failed)\n"},

    {"ib01rd", (PyCFunction)py_ib01rd, METH_VARARGS,
     "Estimate initial state for discrete-time LTI system.\n\n"
     "Given (A,B,C,D) and input/output trajectories, estimates initial state x(0)\n"
     "for: x(k+1) = A*x(k) + B*u(k), y(k) = C*x(k) + D*u(k).\n"
     "Matrix A must be in real Schur form.\n\n"
     "Parameters:\n"
     "  job (str): 'Z' if D is zero, 'N' if D is not zero\n"
     "  n (int): System order (n >= 0)\n"
     "  m (int): Number of inputs (m >= 0)\n"
     "  l (int): Number of outputs (l > 0)\n"
     "  nsmp (int): Number of samples (nsmp >= n)\n"
     "  a (ndarray): N-by-N state matrix A in Schur form (F-order)\n"
     "  b (ndarray): N-by-M input matrix B (F-order)\n"
     "  c (ndarray): L-by-N output matrix C (F-order)\n"
     "  d (ndarray): L-by-M direct transmission D (F-order)\n"
     "  u (ndarray): NSMP-by-M input data (F-order)\n"
     "  y (ndarray): NSMP-by-L output data (F-order)\n"
     "  tol (float): Tolerance for rank estimation (<= 0 uses machine eps)\n\n"
     "Returns:\n"
     "  (x0, rcond, iwarn, info):\n"
     "  - x0: Estimated initial state (n,)\n"
     "  - rcond: Reciprocal condition of triangular factor\n"
     "  - iwarn: Warning (4 = rank-deficient)\n"
     "  - info: Exit code (0=success, 2=SVD failed)\n"},

    {"sb02mt", (PyCFunction)py_sb02mt, METH_VARARGS,
     "Riccati preprocessing - convert coupling weight problems to standard form.\n\n"
     "Computes:\n"
     "  G = B*R^(-1)*B'\n"
     "  A_bar = A - B*R^(-1)*L'\n"
     "  Q_bar = Q - L*R^(-1)*L'\n\n"
     "Parameters:\n"
     "  jobg (str): 'G' to compute G, 'N' to skip\n"
     "  jobl (str): 'Z' if L is zero, 'N' if L is nonzero\n"
     "  fact (str): 'N' R unfactored, 'C' Cholesky, 'U' UdU'/LdL'\n"
     "  uplo (str): 'U' upper triangle, 'L' lower triangle\n"
     "  n (int): Order of A, Q, G (n >= 0)\n"
     "  m (int): Order of R (m >= 0)\n"
     "  a (ndarray or None): Matrix A (n x n, F-order) if jobl='N'\n"
     "  b (ndarray): Matrix B (n x m, F-order)\n"
     "  q (ndarray or None): Matrix Q (n x n, F-order) if jobl='N'\n"
     "  r (ndarray): Matrix R (m x m, F-order)\n"
     "  l (ndarray or None): Matrix L (n x m, F-order) if jobl='N'\n"
     "  g (ndarray or None): Output G (n x n, F-order) if jobg='G'\n\n"
     "Returns:\n"
     "  Variable outputs based on jobg and jobl:\n"
     "  - If jobg='G', jobl='Z': (g, oufact, info)\n"
     "  - If jobg='N', jobl='N': (a, b, q, l, oufact, info)\n"
     "  - If jobg='G', jobl='N': (a, b, q, l, g, oufact, info)\n"},

    {"sb02nd", (PyCFunction)py_sb02nd, METH_VARARGS,
     "Optimal state feedback matrix for optimal control problem.\n\n"
     "Computes:\n"
     "  F = (R + B'XB)^(-1) (B'XA + L')  [discrete-time]\n"
     "  F = R^(-1) (B'X + L')            [continuous-time]\n\n"
     "Parameters:\n"
     "  dico (str): 'D' for discrete-time, 'C' for continuous-time\n"
     "  fact (str): 'N' R unfactored, 'D' R=D'D, 'C' Cholesky, 'U' UdU'/LdL' (continuous only)\n"
     "  uplo (str): 'U' upper triangle, 'L' lower triangle\n"
     "  jobl (str): 'Z' if L is zero, 'N' if L is nonzero\n"
     "  n (int): Order of A, X (n >= 0)\n"
     "  m (int): Number of inputs (m >= 0)\n"
     "  p (int): Rows of D if fact='D' (p >= m for continuous)\n"
     "  rnorm (float): 1-norm of R (required for fact='U' only)\n"
     "  a (ndarray): State matrix A (n x n, F-order), discrete only\n"
     "  b (ndarray): Input matrix B (n x m, F-order)\n"
     "  r (ndarray): Weighting matrix R (m x m or p x m, F-order)\n"
     "  ipiv (ndarray): Pivot indices (m,), input for fact='U'\n"
     "  l (ndarray): Cross weighting L (n x m, F-order) if jobl='N'\n"
     "  x (ndarray): Riccati solution X (n x n, F-order)\n\n"
     "Returns:\n"
     "  (f, r_out, x_out, oufact, rcond, info):\n"
     "  - f: Optimal feedback matrix F (m x n)\n"
     "  - r_out: Factored R or R+B'XB\n"
     "  - x_out: Possibly modified X (for discrete with factored R)\n"
     "  - oufact: Array [fact_type, x_fact_type]\n"
     "  - rcond: Reciprocal condition number\n"
     "  - info: Exit code (0=success, m+1=singular)\n"},

    {"mb01sd", (PyCFunction)py_mb01sd, METH_VARARGS,
     "Scale rows or columns of a matrix by a diagonal matrix.\n\n"
     "Computes one of:\n"
     "  A := diag(R) * A        (jobs='R', row scaling)\n"
     "  A := A * diag(C)        (jobs='C', column scaling)\n"
     "  A := diag(R) * A * diag(C)  (jobs='B', both)\n\n"
     "Parameters:\n"
     "  jobs (str): 'R' for row, 'C' for column, 'B' for both\n"
     "  a (ndarray): M-by-N matrix A (F-order), modified in-place\n"
     "  r (ndarray): Row scale factors, dimension M (not used if jobs='C')\n"
     "  c (ndarray): Column scale factors, dimension N (not used if jobs='R')\n\n"
     "Returns:\n"
     "  a: The scaled matrix\n"},

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
