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

static PyObject* py_tg01fd(PyObject* self, PyObject* args) {
    const char *compq, *compz, *joba;
    i32 l, n, m, p;
    PyObject *a_obj, *e_obj, *b_obj, *c_obj;
    f64 tol;
    PyArrayObject *a_array, *e_array, *b_array, *c_array;
    i32 ranke, rnka22, info;
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
