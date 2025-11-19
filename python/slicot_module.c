/*
 * SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 1996-2025, The SLICOT Team (original Fortran77 code)
 * Copyright (c) 2025, slicot.c contributors (C11 translation)
 */

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdbool.h>
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

    // Create output array with properly allocated data
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

/* Python wrapper for sg03br */
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

/* Module method definitions */
static PyMethodDef SlicotMethods[] = {
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
