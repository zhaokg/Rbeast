#include <math.h>
#include <string.h>

#include "assert.h"
#include "abc_000_warning.h"

#include "abc_ide_util.h"
#include "abc_common.h"
#include "abc_date.h"


#if P_INTERFACE==1  

/////////////////////////////////////////////////////////////////////////////////////////////////
// /  EXPICILTY TO DEINTE ALL THE NEEDED NumPy functions and constants// 
/////////////////////////////////////////////////////////////////////////////////////////////////

#include "limits.h"

void** PyArray_API = NULL;


// Version numbers copied from numpyconfig.h
#define NPY_1_7_API_VERSION 0x00000007
#define NPY_1_8_API_VERSION 0x00000008
#define NPY_1_9_API_VERSION 0x00000009
#define NPY_1_10_API_VERSION 0x0000000a
#define NPY_1_11_API_VERSION 0x0000000a
#define NPY_1_12_API_VERSION 0x0000000a
#define NPY_1_13_API_VERSION 0x0000000b
#define NPY_1_14_API_VERSION 0x0000000c
#define NPY_1_15_API_VERSION 0x0000000c
#define NPY_1_16_API_VERSION 0x0000000d
#define NPY_1_17_API_VERSION 0x0000000d
#define NPY_1_18_API_VERSION 0x0000000d
#define NPY_1_19_API_VERSION 0x0000000d
#define NPY_1_20_API_VERSION 0x0000000e
#define NPY_1_21_API_VERSION 0x0000000e
#define NPY_1_22_API_VERSION 0x0000000f
#define NPY_1_23_API_VERSION 0x00000010
#define NPY_1_24_API_VERSION 0x00000010
#define NPY_1_25_API_VERSION 0x00000011
#define NPY_2_0_API_VERSION 0x00000012
#define NPY_2_1_API_VERSION 0x00000013
#define NPY_2_2_API_VERSION 0x00000013
#define NPY_2_3_API_VERSION 0x00000014
#define NPY_2_4_API_VERSION 0x00000015

typedef unsigned char npy_bool;

#define NPY_SIZEOF_SHORT SIZEOF_SHORT
#define NPY_SIZEOF_INT   SIZEOF_INT
#define NPY_SIZEOF_LONG SIZEOF_LONG
#define NPY_SIZEOF_FLOAT 4
#define NPY_SIZEOF_COMPLEX_FLOAT 8
#define NPY_SIZEOF_DOUBLE 8
#define NPY_SIZEOF_COMPLEX_DOUBLE 16
#define NPY_SIZEOF_LONGDOUBLE 8
#define NPY_SIZEOF_COMPLEX_LONGDOUBLE 16
#define NPY_SIZEOF_PY_INTPTR_T 8
#define NPY_SIZEOF_OFF_T 4
#define NPY_SIZEOF_PY_LONG_LONG 8
#define NPY_SIZEOF_LONGLONG 8

#define NPY_SIZEOF_CHAR 1
#define NPY_SIZEOF_BYTE 1
#define NPY_SIZEOF_DATETIME 8
#define NPY_SIZEOF_TIMEDELTA 8
#define NPY_SIZEOF_INTP  NPY_SIZEOF_PY_INTPTR_T
#define NPY_SIZEOF_UINTP NPY_SIZEOF_PY_INTPTR_T
#define NPY_SIZEOF_HALF 2


#define NPY_BITSOF_BOOL (sizeof(npy_bool) * CHAR_BIT)
#define NPY_BITSOF_CHAR CHAR_BIT
#define NPY_BITSOF_BYTE (NPY_SIZEOF_BYTE * CHAR_BIT)
#define NPY_BITSOF_SHORT (NPY_SIZEOF_SHORT * CHAR_BIT)
#define NPY_BITSOF_INT (NPY_SIZEOF_INT * CHAR_BIT)
#define NPY_BITSOF_LONG (NPY_SIZEOF_LONG * CHAR_BIT)
#define NPY_BITSOF_LONGLONG (NPY_SIZEOF_LONGLONG * CHAR_BIT)
#define NPY_BITSOF_INTP (NPY_SIZEOF_INTP * CHAR_BIT)
#define NPY_BITSOF_HALF (NPY_SIZEOF_HALF * CHAR_BIT)
#define NPY_BITSOF_FLOAT (NPY_SIZEOF_FLOAT * CHAR_BIT)
#define NPY_BITSOF_DOUBLE (NPY_SIZEOF_DOUBLE * CHAR_BIT)
#define NPY_BITSOF_LONGDOUBLE (NPY_SIZEOF_LONGDOUBLE * CHAR_BIT)
 
  
#if NPY_BITSOF_LONG == 8
#define NPY_INT8 NPY_LONG
#define NPY_UINT8 NPY_ULONG
#elif NPY_BITSOF_LONG == 16
#define NPY_INT16 NPY_LONG
#define NPY_UINT16 NPY_ULONG
#elif NPY_BITSOF_LONG == 32
#define NPY_INT32 NPY_LONG
#define NPY_UINT32 NPY_ULONG 
#elif NPY_BITSOF_LONG == 64
#define NPY_INT64 NPY_LONG
#define NPY_UINT64 NPY_ULONG
#elif NPY_BITSOF_LONG == 128
#define NPY_INT128 NPY_LONG
#define NPY_UINT128 NPY_ULONG
#endif

#if NPY_BITSOF_LONGLONG == 8
#  ifndef NPY_INT8
#    define NPY_INT8 NPY_LONGLONG
#    define NPY_UINT8 NPY_ULONGLONG
#  endif
#elif NPY_BITSOF_LONGLONG == 16
#  ifndef NPY_INT16
#    define NPY_INT16 NPY_LONGLONG
#    define NPY_UINT16 NPY_ULONGLONG
#  endif
#elif NPY_BITSOF_LONGLONG == 32
#  ifndef NPY_INT32
#    define NPY_INT32 NPY_LONGLONG
#    define NPY_UINT32 NPY_ULONGLONG
#  endif
#elif NPY_BITSOF_LONGLONG == 64
#  ifndef NPY_INT64
#    define NPY_INT64 NPY_LONGLONG
#    define NPY_UINT64 NPY_ULONGLONG
#  endif
#elif NPY_BITSOF_LONGLONG == 128
#  ifndef NPY_INT128
#    define NPY_INT128 NPY_LONGLONG
#    define NPY_UINT128 NPY_ULONGLONG
#  endif
#elif NPY_BITSOF_LONGLONG == 256
#  define NPY_INT256 NPY_LONGLONG
#  define NPY_UINT256 NPY_ULONGLONG
#endif

#if NPY_BITSOF_INT == 8
#ifndef NPY_INT8
#define NPY_INT8 NPY_INT
#define NPY_UINT8 NPY_UINT
#endif
#elif NPY_BITSOF_INT == 16
#ifndef NPY_INT16
#define NPY_INT16 NPY_INT
#define NPY_UINT16 NPY_UINT
#endif
#elif NPY_BITSOF_INT == 32
#ifndef NPY_INT32
#define NPY_INT32 NPY_INT
#define NPY_UINT32 NPY_UINT 
#endif
#elif NPY_BITSOF_INT == 64
#ifndef NPY_INT64
#define NPY_INT64 NPY_INT
#define NPY_UINT64 NPY_UINT
#endif
#elif NPY_BITSOF_INT == 128
#ifndef NPY_INT128
#define NPY_INT128 NPY_INT
#define NPY_UINT128 NPY_UINT
#endif
#endif

#if NPY_BITSOF_SHORT == 8
#ifndef NPY_INT8
#define NPY_INT8 NPY_SHORT
#define NPY_UINT8 NPY_USHORT
#endif
#elif NPY_BITSOF_SHORT == 16
#ifndef NPY_INT16
#define NPY_INT16 NPY_SHORT
#define NPY_UINT16 NPY_USHORT
#endif
#elif NPY_BITSOF_SHORT == 32
#ifndef NPY_INT32
#define NPY_INT32 NPY_SHORT
#define NPY_UINT32 NPY_USHORT
#endif
#elif NPY_BITSOF_SHORT == 64
#ifndef NPY_INT64
#define NPY_INT64 NPY_SHORT
#define NPY_UINT64 NPY_USHORT
#endif
#elif NPY_BITSOF_SHORT == 128
#ifndef NPY_INT128
#define NPY_INT128 NPY_SHORT
#define NPY_UINT128 NPY_USHORT
#endif
#endif


#if NPY_BITSOF_CHAR == 8
#ifndef NPY_INT8
#define NPY_INT8 NPY_BYTE
#define NPY_UINT8 NPY_UBYTE 
#endif
#elif NPY_BITSOF_CHAR == 16
#ifndef NPY_INT16
#define NPY_INT16 NPY_BYTE
#define NPY_UINT16 NPY_UBYTE
#endif
#elif NPY_BITSOF_CHAR == 32
#ifndef NPY_INT32
#define NPY_INT32 NPY_BYTE
#define NPY_UINT32 NPY_UBYTE
#endif
#elif NPY_BITSOF_CHAR == 64
#ifndef NPY_INT64
#define NPY_INT64 NPY_BYTE
#define NPY_UINT64 NPY_UBYTE
#endif
#elif NPY_BITSOF_CHAR == 128
#ifndef NPY_INT128
#define NPY_INT128 NPY_BYTE
#define NPY_UINT128 NPY_UBYTE
#endif
#endif



#if NPY_BITSOF_DOUBLE == 32
#ifndef NPY_FLOAT32
#define NPY_FLOAT32 NPY_DOUBLE
#define NPY_COMPLEX64 NPY_CDOUBLE
#endif
#elif NPY_BITSOF_DOUBLE == 64
#ifndef NPY_FLOAT64
#define NPY_FLOAT64 NPY_DOUBLE
#define NPY_COMPLEX128 NPY_CDOUBLE
#endif
#elif NPY_BITSOF_DOUBLE == 80
#ifndef NPY_FLOAT80
#define NPY_FLOAT80 NPY_DOUBLE
#define NPY_COMPLEX160 NPY_CDOUBLE
#endif
#elif NPY_BITSOF_DOUBLE == 96
#ifndef NPY_FLOAT96
#define NPY_FLOAT96 NPY_DOUBLE
#define NPY_COMPLEX192 NPY_CDOUBLE
#endif
#elif NPY_BITSOF_DOUBLE == 128
#ifndef NPY_FLOAT128
#define NPY_FLOAT128 NPY_DOUBLE
#define NPY_COMPLEX256 NPY_CDOUBLE
#endif
#endif



#if NPY_BITSOF_FLOAT == 32
#ifndef NPY_FLOAT32
#define NPY_FLOAT32 NPY_FLOAT
#define NPY_COMPLEX64 NPY_CFLOAT
#endif
#elif NPY_BITSOF_FLOAT == 64
#ifndef NPY_FLOAT64
#define NPY_FLOAT64 NPY_FLOAT
#define NPY_COMPLEX128 NPY_CFLOAT
#endif
#elif NPY_BITSOF_FLOAT == 80
#ifndef NPY_FLOAT80
#define NPY_FLOAT80 NPY_FLOAT
#define NPY_COMPLEX160 NPY_CFLOAT
#endif
#elif NPY_BITSOF_FLOAT == 96
#ifndef NPY_FLOAT96
#define NPY_FLOAT96 NPY_FLOAT
#define NPY_COMPLEX192 NPY_CFLOAT
#endif
#elif NPY_BITSOF_FLOAT == 128
#ifndef NPY_FLOAT128
#define NPY_FLOAT128 NPY_FLOAT
#define NPY_COMPLEX256 NPY_CFLOAT
#endif
#endif


typedef Py_intptr_t  npy_intp;
typedef Py_uintptr_t npy_uintp;

/* For specifying array memory layout or iteration order */
typedef enum {
    NPY_ANYORDER = -1,  /* Fortran order if inputs are all Fortran, C otherwise */
    NPY_CORDER = 0,     /* C order */
    NPY_FORTRANORDER = 1,        /* Fortran order */
    NPY_KEEPORDER = 2            /* An order as close to the inputs as possible */
} NPY_ORDER;

typedef struct {
    npy_intp* ptr;
    int len;
} PyArray_Dims;

#define NPY_ARRAY_F_CONTIGUOUS    0x0002
#define NPY_ARRAY_ALIGNED         0x0100
#define NPY_ARRAY_WRITEABLE       0x0400
#define NPY_ARRAY_BEHAVED         (NPY_ARRAY_ALIGNED      |  NPY_ARRAY_WRITEABLE)
#define NPY_ARRAY_FARRAY          (NPY_ARRAY_F_CONTIGUOUS |  NPY_ARRAY_BEHAVED)


#define PyArray_GetNDArrayCVersion     (*(unsigned int (*)(void))    PyArray_API[0])  // Numpy's C API Version
#define PyArray_Type           (*(PyTypeObject *)             PyArray_API[2])
#define PyArray_DescrFromType  (*(PyArray_Descr * (*)(int))   PyArray_API[45])
#define PyArray_Size           (*(npy_intp (*)(PyObject *))   PyArray_API[59])
#define PyArray_NewFromDescr   (*(PyObject * (*)(PyTypeObject *, PyArray_Descr *, int, npy_intp const *, npy_intp const *, void *, int, PyObject *)) \
                              PyArray_API[94])		 
#define PyArray_FromArray      (*(PyObject * (*)(PyArrayObject *, PyArray_Descr *, int))                   PyArray_API[109])		 
#define PyArray_CheckFromAny   (*(PyObject * (*)(PyObject *, PyArray_Descr *, int, int, int, PyObject *))  PyArray_API[108])	
#define PyArray_Newshape       (*(PyObject * (*)(PyArrayObject *, PyArray_Dims *, NPY_ORDER))              PyArray_API[135])		 
#define PyArray_GetPtr         (*(void * (*)(PyArrayObject *, npy_intp const*))                            PyArray_API[160])
#define PyArray_MultiplyList   (*(npy_intp (*)(npy_intp const *, int))                                     PyArray_API[158])


#define PyArray_DIMS(obj)       (((PyArrayObject_fields *)(obj))->dimensions)
#define PyArray_FROM_OF(m,flags) PyArray_CheckFromAny(m, NULL, 0, 0, flags, NULL)		 
#define PyArray_Check(op)        PyObject_TypeCheck(op, &PyArray_Type)		 
#define PyArray_SIZE(m)          PyArray_MultiplyList(PyArray_DIMS(m), PyArray_NDIM(m))
#define NPY_ATTR_DEPRECATE(text)

enum NPY_TYPES {
    NPY_BOOL = 0, NPY_BYTE, NPY_UBYTE, NPY_SHORT, NPY_USHORT,
    NPY_INT, NPY_UINT, NPY_LONG, NPY_ULONG, NPY_LONGLONG, NPY_ULONGLONG,
    NPY_FLOAT, NPY_DOUBLE, NPY_LONGDOUBLE, NPY_CFLOAT, NPY_CDOUBLE, NPY_CLONGDOUBLE,
    NPY_OBJECT = 17, NPY_STRING, NPY_UNICODE, NPY_VOID,
    /* New 1.6 types appended, may be integrated into the above in 2.0.*/
    NPY_DATETIME, NPY_TIMEDELTA, NPY_HALF, NPY_NTYPES, NPY_NOTYPE,
    NPY_CHAR NPY_ATTR_DEPRECATE("Use NPY_STRING"),
    NPY_USERDEF = 256,  /* leave room for characters */
    NPY_NTYPES_ABI_COMPATIBLE = 21 /* The number of types not including the new 1.6 types */
};

typedef Py_hash_t npy_hash_t;
typedef struct  PyArrayObject PyArrayObject;
typedef struct  PyArray_ArrFuncs PyArray_ArrFuncs;
typedef struct  NpyAuxData  NpyAuxData;
typedef struct  NpyAuxData  NpyAuxData;

/* Migrate from Num[y 1.x to Numpy 2.x

 If old code directly accessed fields like descr->elsize, NumPy recommends replacing that with accessors such as PyArray_
ITEMSIZE(arr) or PyDataType_ELSIZE(descr). Also, elsize is now conceptually npy_intp, not int.

 The defintion of _PyArray_Descr below is superfical becauwe we do not directlya access fields rahter we use accesor 
 functions.
*/

typedef struct _PyArray_Descr {
    PyObject_HEAD
        /*
         * the type object representing an
         * instance of this type -- should not
         * be two type_numbers with the same type
         * object.
         */
        PyTypeObject* typeobj;
    char kind;  /* kind for this type */
    char type; /* unique-character representing this type */
    char byteorder; /*  '>' (big), '<' (little), '|'  (not-applicable), or '=' (native).      */
    char flags;    /* flags describing data type */
    int type_num;  /* number representing this type */
    int elsize;    /* element size (itemsize) for this type */
    int alignment; /* alignment needed for this type */
    struct _arr_descr* subarray;   // Non-NULL if this type is is an array (C-contiguous)of some other type
    PyObject* fields; // The fields dictionary for this type For statically defined descr this is always Py_None
    PyObject* names; //An ordered tuple of field names or NULL if no fields are defined
    PyArray_ArrFuncs* f;  // a table of functions specific for each basic data descriptor
    PyObject* metadata;          /* Metadata about this dtype */
    NpyAuxData* c_metadata; // Metadata specific to the C implementation  of the particular dtype. This was added for NumPy 1.7.0.
    npy_hash_t hash; // Cached hash value (-1 if not yet computed).This was added for NumPy 2.0.0.
} PyArray_Descr;



typedef struct tagPyArrayObject_fields {
    PyObject_HEAD
        char* data;  // Pointer to the raw data buffer    
    int nd;      // The number of dimensions, also called 'ndim'    
    npy_intp* dimensions; /* The size in each dimension, also called 'shape' */
    npy_intp* strides;     //Number of bytes to jump to get to the next element in each dimension
    /*
     * This object is decref'd upon
     * deletion of array. Except in the
     * case of WRITEBACKIFCOPY which has
     * special handling.
     *
     * For views it points to the original
     * array, collapsed so no chains of
     * views occur.
     *
     * For creation from buffer object it
     * points to an object that should be
     * decref'd on deletion
     *
     * For WRITEBACKIFCOPY flag this is an
     * array to-be-updated upon calling
     * PyArray_ResolveWritebackIfCopy
     */
    PyObject* base;
    PyArray_Descr* descr;     /* Pointer to type structure */
    int flags;                /* Flags describing array -- see below */
    PyObject* weakreflist;    /* For weak references */
    void* _buffer_info;       /* private buffer info, tagged to allow warning */
} PyArrayObject_fields;




static int       PyArray_TYPE(const PyArrayObject* arr) { return ((PyArrayObject_fields*)arr)->descr->type_num; }
static void*     PyArray_DATA(PyArrayObject* arr) { return ((PyArrayObject_fields*)arr)->data; }
static npy_intp  PyArray_ITEMSIZE(const PyArrayObject* arr) { return ((PyArrayObject_fields*)arr)->descr->elsize; }
static int       PyArray_NDIM(const PyArrayObject* arr) { return ((PyArrayObject_fields*)arr)->nd; }


/////////////////////////////////////////////////////////////////////////////////////////////////
// / DONE WITH  "  EXPICILTY TO DEINTE ALL THE NEEDED NumPy functions and constants// "
/////////////////////////////////////////////////////////////////////////////////////////////////




int  JDN_to_DateNum(int jdn) {
    return jdn - JulianDayNum_from_civil_ag1(1, 1, 1);
}

void StdouFlush(void) {
 // https://stackoverflow.com/questions/69247396/python-c-api-how-to-flush-stdout-and-stderr

    // Two ways to get the sys.stdout.flush

    /*
    PyObject* sysModule = PyImport_AddModule("sys");  //borrowedf
    if (sysModule == NULL) {
        PySys_WriteStdout("Not loaded yet\n ");
        sysModule = PyImport_ImportModule("sys");         // new ref: We don't decref it because we want to re-use it if beast is called again
    }
    PyObject* sout  = PyObject_GetAttrString(sysModule, "stdout");  //new ref
    PyObject* flush = PyObject_GetAttrString(sout, "flush");  //new ref
    */

    PyObject* sysout   = PySys_GetObject("__stdout__");           //borrowedf
    PyObject* pfflush  = PyObject_GetAttrString(sysout, "flush");  //new ref

    int       ok  = PyCallable_Check(pfflush); 
    PyObject* out = PyObject_CallObject(pfflush, NULL);		 //new ref
     
    Py_XDECREF(pfflush);
    Py_XDECREF(out);
}


PyObject* currentModule;
PyObject* classOutput  = NULL;

PyObject* setClassObjects(PyObject* self, PyObject* args){

    if (PyTuple_Size(args) == 1) {
        classOutput= PyTuple_GetItem(args,0);
    }
    return Py_None;
}

#define  Pob PyObject *

void* CvtToPyArray_NewRef(VOIDPTR Y) {

    if (PyArray_Check(Y)) {
    // Accessing view of a NumPy array using the C API
    // https://stackoverflow.com/questions/35000565/accessing-view-of-a-numpy-array-using-the-c-api
        PyArray_Descr* reqDescr = PyArray_DescrFromType(PyArray_TYPE(Y));
        Pob            out      = PyArray_FromArray(Y, reqDescr, NPY_ARRAY_FARRAY);
        return out;
    }   else {
        Pob            out = PyArray_FROM_OF(Y, NPY_ARRAY_FARRAY);
        return out;
    }

    return NULL;
}

typedef PyObject* (*PyGetItem)(PyObject*, Py_ssize_t);
 
static void* PyGetDict(void *ptr) {

    if (PyDict_Check(ptr)) {
        return ptr;
    }

    if (PyList_Check(ptr) || PyTuple_Check(ptr)) {
        return NULL;
    }
        
    if (PyLong_Check(ptr) || PyFloat_Check(ptr)) {
        return NULL;
    }

    if (PyObject_IsInstance(ptr, (PyObject *) & PyBaseObject_Type)) {
     
        // instance's own dict
        PyObject** dictpr = _PyObject_GetDictPtr(ptr);
        if (dictpr != NULL) {
            PyObject* Dict = *dictpr;
            return Dict;
        }

        // class' dictonary
        // not a NumpyArray
        if (!PyArray_Check(ptr) && Py_TYPE(ptr)->tp_dict != NULL) {
            PyObject* Dict = Py_TYPE(ptr)->tp_dict;
            return Dict;
        }
    }

    return NULL;
}

static void* PyGetDictItem(void *ptr, int idx) {
           
    PyObject* values = PyDict_Values(ptr);              // new ref
    Pob       item   = PyList_GetItem(values,idx);
    Py_DECREF(values);
    return item;
}

static void* PyGetDictItemString(void *ptr, char *fldname) {

    PyObject* item = PyDict_GetItemString(ptr, fldname); //borrowed ref
    if (item) { 
        return item; 
    }
   
    PyObject* keys = PyDict_Keys(ptr);              // new ref
    int       n    = (int) PyList_Size(keys);
    item           = NULL;
    for (int i = 0; i < n; i++) {
        char      tmpName[100 + 1];
        PyObject* tmpkey = PyList_GetItem(keys, i);  // borrowed ref
        int       len    = GetCharArray(tmpkey, tmpName, 100);
        if (len > 0 && strcicmp(tmpName, fldname) == 0) {
            item = PyDict_GetItem(ptr, tmpkey); //borrwed ref
            break;
        }
    }
    Py_DECREF(keys);
    return item;

}

static void* PyGetDictItemString123(void *ptr, char *fldname, int nPartial) {

    PyObject* item = PyDict_GetItemString(ptr, fldname); //borrowed ref
    if (item) { 
        return item; 
    }

    PyObject* keys = PyDict_Keys(ptr);              // new ref
    int       n    = (int) PyList_Size(keys);
    item = NULL;
    for (int i = 0; i < n; i++) {
        char      tmpName[100 + 1];
        PyObject* tmpkey = PyList_GetItem(keys, i);  // borrowed ref
        int       len    = GetCharArray(tmpkey, tmpName, 100);
        if (len > 0 && strcicmp_nfirst(tmpName, fldname,nPartial) == 0) {
            item = PyDict_GetItem(ptr, tmpkey); //borrwed ref
            break;
        }
    }
    Py_DECREF(keys);
    return item;

}

 
int IsClass(void* ptr, char* class) { return 0; }
int IsCell(void* ptr) { return 0; }
int IsChar(void* ptr) {   

    if (ptr == NULL)                  return 0;
    if (PyUnicode_Check(ptr))         return 1;    
    if (PyArray_Check(ptr) && PyArray_TYPE(ptr) == NPY_STRING)      return 1;    
    if (PyArray_Check(ptr) && PyArray_TYPE(ptr) == NPY_UNICODE)     return 1;
    

    PyGetItem pyGetItem = NULL;
    if      (PyList_Check(ptr))        pyGetItem = PyList_GetItem;
    else if (PyTuple_Check(ptr))       pyGetItem = PyTuple_GetItem; //PyDict_GetItem;
    
    
    if(pyGetItem) {      
        int sz = (int) PyObject_Size(ptr);
        int ok = 1;
        for (int i = 0; i < sz; ++i) {
            if ( !PyUnicode_Check(pyGetItem(ptr, i) ) ) {
                ok = 0; 
                break;
            }
        }
        if (ok) { return 1; }
    }

    if (PyArray_Check(ptr) && (PyArray_TYPE(ptr) == NPY_STRING || PyArray_TYPE(ptr) == NPY_UNICODE)  ){
        return 1;
    }

    if (PyArray_Check(ptr)  && PyArray_TYPE(ptr)==NPY_OBJECT) {
        void** base = (void**) PyArray_DATA(ptr);
        int sz      = (int)    PyArray_Size(ptr);
        int ok      = 1;
        for (int i = 0; i < sz; ++i) {
            if (!PyUnicode_Check(base[i])) {
                ok = 0;
                break;
            }
        }
        if (ok) { return 1; }
    }

    return 0;
}
int IsStruct(void* ptr) {

    if (ptr == NULL)  return 0;
    if (PyList_Check(ptr) || PyDict_Check(ptr) || PyTuple_Check(ptr) )     return 1L;
        
    Pob dict = PyGetDict(ptr);
    if (dict) {
        if (PyUnicode_Check(ptr)) {
        // a string also has an dict but is not treated as an object here
            return 0;
        }  else {
            return 1;
        }        
    }

    return 0L;
}
int IsEmpty( void* ptr) {  return (ptr == Py_None || GetNumberOfElements(ptr) == 0);}

//https://numpy.org/devdocs/reference/c-api/dtype.html#c.NPY_COMPLEX64
//https://docs.python.org/3/c-api/concrete.html
int IsSingle(void* ptr) {
    //PyLong_Check(ptr),
    //PyObject_IsInstance(ptr, (PyObject*)&PyLong_Type) 
    
    /**********************************/
    //Python's float is a double type
    /**********************************/
    if (PyArray_Check(ptr) && PyArray_TYPE(ptr) == NPY_FLOAT32) {
        return 1L;        
    }
    return 0L;
}

int IsDouble(void* ptr) {
    //PyLong_Check(ptr),
    //PyObject_IsInstance(ptr, (PyObject*)&PyLong_Type) 
 
    if (PyFloat_Check(ptr))                                         return 1L; 
    if (PyArray_Check(ptr) && PyArray_TYPE(ptr) == NPY_FLOAT64)     return 1L;

    PyGetItem pyGetItem = NULL;
    if      (PyList_Check(ptr))       pyGetItem = PyList_GetItem;
    else if (PyTuple_Check(ptr))      pyGetItem = PyTuple_GetItem; //PyDict_GetItem;

    if (pyGetItem) {
        int sz = (int) PyObject_Size(ptr);
        int ok = 1;
        for (int i = 0; i < sz; ++i) {
            if (!PyFloat_Check(pyGetItem(ptr, i))) {
                ok = 0;
                break;
            }
        }
        if (ok) { return 1; }
    }
 

    if (PyArray_Check(ptr) && PyArray_TYPE(ptr) == NPY_OBJECT) {
        void** base = (void**)PyArray_DATA(ptr);
        int sz = (int)PyArray_Size(ptr);
        int ok = 1;
        for (int i = 0; i < sz; ++i) {
            if (!PyFloat_Check(base[i])) {
                ok = 0;
                break;
            }
        }
        if (ok) { return 1; }
    }

    return 0;
}
int IsInt16(void* ptr) {   
    if (PyArray_Check(ptr) && PyArray_TYPE(ptr) == NPY_INT16) {
        return 1L;
    }
    return 0;
}
int IsInt32(void* ptr) {
    //PyLong_Check(ptr),
    //PyObject_IsInstance(ptr, (PyObject*)&PyLong_Type) 
    if (PyLong_Check(ptr)) {
        return 1L;
    }
    if (PyArray_Check(ptr) && PyArray_TYPE(ptr) == NPY_INT32) {
        return 1L;
    }

    PyGetItem pyGetItem = NULL;
    if      (PyList_Check(ptr))       pyGetItem = PyList_GetItem;
    else if (PyTuple_Check(ptr))      pyGetItem = PyTuple_GetItem; //PyDict_GetItem;
    if (pyGetItem) {
        int sz = (int) PyObject_Size(ptr);
        int ok = 1;
        for (int i = 0; i < sz; ++i) {
            if ( !PyLong_Check( pyGetItem(ptr, i) ) ) {
                ok = 0;
                break;
            }
        }
        if (ok) { return 1; }
    }
   
    if (PyArray_Check(ptr) && PyArray_TYPE(ptr) == NPY_OBJECT) {
        void** base = (void**)PyArray_DATA(ptr);
        int sz = (int)PyArray_Size(ptr);
        int ok = 1;
        for (int i = 0; i < sz; ++i) {
            if (!PyLong_Check(base[i])) {
                ok = 0;
                break;
            }
        }
        if (ok) { return 1; }
    }
    return 0;
 
}

int IsInt64(void* ptr) {
     if (PyArray_Check(ptr) && PyArray_TYPE(ptr) == NPY_INT64) {
         return 1L;
    }
    return 0;
}

int IsLogical(void* ptr) {
     // a subclass of integer
    if (PyBool_Check(ptr)) {
        return 1L;
    }  
    if (PyArray_Check(ptr) && PyArray_TYPE(ptr) == NPY_BOOL) {
        return 1L;
    }

    PyGetItem pyGetItem = NULL;
    if      (PyList_Check(ptr))       pyGetItem = PyList_GetItem;
    else if (PyTuple_Check(ptr))      pyGetItem = PyTuple_GetItem; //PyDict_GetItem;

    if (pyGetItem) {
        int sz = (int) PyObject_Size(ptr);
        int ok = 1;
        for (int i = 0; i < sz; ++i) {
            if (!PyBool_Check(pyGetItem(ptr, i))) {
                ok = 0;
                break;
            }
        }
        if (ok) { return 1; }
    }
     
    if (PyArray_Check(ptr) && PyArray_TYPE(ptr) == NPY_OBJECT) {
        void** base = (void**)PyArray_DATA(ptr);
        int sz = (int)PyArray_Size(ptr);
        int ok = 1;
        for (int i = 0; i < sz; ++i) {
            if (!PyBool_Check(base[i])) {
                ok = 0;
                break;
            }
        }
        if (ok) { return 1; }
    }

    return 0;
}

int IsNumeric(void* ptr) {
    return IsDouble(ptr) || IsSingle(ptr) || IsInt32(ptr)|| IsInt16(ptr) || IsInt64(ptr) || IsLogical(ptr);
}


VOID_PTR GetFieldByIdx(VOID_PTR strucVar, I32 ind) {
 
    if (PyList_Check(strucVar)) {
        PyObject* out = PyList_GetItem(strucVar, ind); //borrowed ref
        return out;
    }  

    if (PyTuple_Check(strucVar)) {
        PyObject* out = PyTuple_GetItem(strucVar, ind); //borrowed ref
        return out;
    }

    if (PyDict_Check(strucVar)) {
        Pob list = PyDict_Values(strucVar);//  new ref
        Pob item = PyList_GetItem(list, ind); // borrowed ref
        Py_DECREF(list);
        return item;
    }

    Pob dict = PyGetDict(strucVar);
    if (dict) { 
            Pob list = PyDict_Values(strucVar);//  new ref
            Pob item = PyList_GetItem(list, ind); // borrowed ref
            Py_DECREF(list);
            return item;   
    }

    return NULL;
}


void  GetFieldNameByIdx(VOID_PTR strucVar, I32 ind0, char *str, int buflen) {

    Pob dict   = PyGetDict(strucVar);
    
    if (dict) {
        PyObject*  keys = PyDict_Keys(dict);              // new ref
        Py_ssize_t n    = PyList_Size(keys);

        PyObject* tmpkey = PyList_GetItem(keys, ind0);  // borrowed ref
        if (IsChar(tmpkey)) {
            int       len = GetCharArray(tmpkey, str, buflen);
        }   else {
            str[0] = 0;
        }        

        Py_DECREF(keys);
    }    else {
        str[0] = 0;
    }
    
    
}

void *CreateStructVar(FIELD_ITEM* fieldList, int nfields); //Done

void  DestoryStructVar(VOID_PTR strutVar) {  Py_XDECREF(strutVar);}

void  RemoveField(FIELD_ITEM* fieldList, int nfields, char* fieldName); // Done. IDE-independent

void AddStringAttribute(VOID_PTR listVar, const char* field, const char* value) {

    Pob dict = PyGetDict(listVar);
    if (dict) {
        Pob pyValue = PyUnicode_FromString(value);
        PyObject_SetAttrString(listVar, field, pyValue); // do not steal ref and value is increfed.
        Py_XDECREF(pyValue);  // decref it so
    }

}
void AddIntegerAttribute(VOID_PTR listVar, const char* field, I32 value) {

    Pob dict = PyGetDict(listVar);
    if (dict) {
        Pob pyValue = PyLong_FromLong(value);
        PyObject_SetAttrString(listVar, field, pyValue); // do not steal ref and value is increfed.
        Py_XDECREF(pyValue);  // decref it so
    }

}
void RemoveAttribute(VOID_PTR listVar, const char* field) {

}

I32   GetConsoleWidth() {
    return 85;

    PyObject* osModule         = PyImport_AddModule("os");  //borrowedf
    if (osModule == NULL) {
        osModule = PyImport_ImportModule("os");  // new ref
        // We don't decref it because we want to re-use it if beast is called again
    }

    if (osModule == NULL) {
        return 80;
    }

    PyObject* get_console_size = PyObject_GetAttrString(osModule, "get_terminal_size"); //new ref

 
    // Pob result     = PyObject_CallObject(get_console_size, 0);                               // new ref
    // Pob result     = PyObject_CallFunction(get_console_size, NULL);
    Pob result     = PyObject_CallMethod(osModule, "get_terminal_size",NULL);                    // new ref
    Pob col        = PyObject_GetAttrString(result,"columns");  // new ref
    
    int width = PyLong_AsLong(col);

    Py_XDECREF(get_console_size);
    Py_XDECREF(result);
    Py_XDECREF(col);
    return width;
}


I32 GetCharArray(void* ptr, char* dst, int n) {
    dst[0] = 0;
    if (ptr ==NULL || !IsChar(ptr)) return 0;

    if (PyUnicode_Check(ptr)) {
        //int len;
        /* STUPID ME. it took me one day to figure out this bug:
        * if using int, the upper exta 32-bits will be reset to zero,
        * this will corrupt the pushed 'si' regiser, the si regiseter
        * seems to hold the prhs passed in the mainFunction
        */
        Py_ssize_t  len;
        const   char* str = PyUnicode_AsUTF8AndSize(ptr, &len);
        strncpy(dst, str, n);
        return (I32) len;
    }


    int idx = 0;
    return  GetCharVecElem(ptr, idx, dst, n);
 
    
}

I32 GetCharVecElem(void* ptr, int idx, char* dst, int n) {

    Py_ssize_t len = 0; // it has to be of Py_ssize_t type,, as needed for PyUnicode_AsUTF8AndSize
     
    PyObject* tmpItem = NULL;
    if (PyUnicode_Check(ptr) && idx == 0) {
        tmpItem = ptr;
    } else   if (PyList_Check(ptr)) {
        tmpItem = PyList_GetItem(ptr, idx);
    } else  if (PyTuple_Check(ptr)) {
        tmpItem = PyTuple_GetItem(ptr, idx);
    }

    // It is a list or dict item
    if (tmpItem && PyUnicode_Check(tmpItem)) {   
        const char* str = PyUnicode_AsUTF8AndSize(tmpItem, &len);
        len = min(len,n - 1);
        memcpy(dst, str, len);
        dst[len] = 0;
        return (I32) len;
    }

    //If a Numpy string array
    if (PyArray_Check(ptr) && PyArray_TYPE(ptr)==NPY_STRING) {
        const char* base   = PyArray_DATA(ptr);
        int         elsize = PyArray_ITEMSIZE(ptr);
        char       *str    = base + elsize * idx;
        len     = min(elsize, n-1);
        memcpy(dst, str, len);
        dst[len] = 0;
        return (I32) len;
    }

    //If a Numpy Unicode array
    if (PyArray_Check(ptr) && PyArray_TYPE(ptr) == NPY_UNICODE) {
        char*   base   = PyArray_DATA(ptr);
        size_t  elsize = PyArray_ITEMSIZE(ptr);
        char* str = base + elsize * idx;
        len = min(elsize/4, n - 1);  // Assuming a Unicode char has 4 bytse
        for (int i = 0; i < len; i++) {
            dst[i] = str[i * 4];
        }
        dst[len] = 0;
        return (I32) len;
    }

    //If a Numpy object array
    if (PyArray_Check(ptr) && PyArray_TYPE(ptr) == NPY_OBJECT) {
        void**   base   = PyArray_DATA(ptr);
        size_t   elsize = PyArray_ITEMSIZE(ptr);

        void* tmpItem = base[idx];

        // It is a list or dict item
        if (tmpItem && PyUnicode_Check(tmpItem)) {
            const char* str = PyUnicode_AsUTF8AndSize(tmpItem, &len);
            len = min(len, n - 1);
            memcpy(dst, str, len);
            dst[len] = 0;
            return (I32) len;
        }
 
    }

    dst[len] = 0;
    return 0;
}


void* GetField123(const void* structVar, char* fname, int nPartial) {

    Pob dict = NULL;

    if (PyDict_Check((void*)structVar)) {
        dict = (Pob) structVar;
    }   else {
        dict = PyGetDict( (void *) structVar);
    }
 
    if (dict) {
        return PyGetDictItemString123(dict, fname, nPartial);
    } 

    return NULL;
}
void* GetField(const void* structVar, char* fname) {
 
    if (PyDict_Check(structVar)) {
        PyObject* item = PyGetDictItemString((void*)structVar, fname);
        return item;
    
    }  else {
 
        if (PyObject_HasAttrString((void*)structVar, fname)) {
            PyObject *attr= PyObject_GetAttrString((void*)structVar, fname); // new ref
            Py_DECREF(attr);  // so the returned 
            return  attr;
        }        
    }

    return NULL;

}

static F64 NumpyGetNumericElemt(void* ptr, int ind) {

    int    ndim = PyArray_NDIM(ptr);
    if (ndim != 1) {
        return getNaN();
    }

    npy_intp indices[] = { ind };    
    void     *dptr       = PyArray_GetPtr(ptr, indices);
    int      dtype      = PyArray_TYPE(ptr);
    if (dtype == NPY_INT16) { 
        return *(int16_t*)dptr;
    }
    else if (dtype == NPY_INT32) {
        return *(int32_t*)dptr;
    }
    else if (dtype == NPY_INT64) {
        return (F64)  * (int64_t*)dptr;
    }
    else if (dtype == NPY_FLOAT) {
        return *(float*)dptr;
    }
    else if (dtype == NPY_DOUBLE) {
        return *(double*)dptr;
    }

    return getNaN();
}

F64   GetScalar(const void* ptr) {
    
    Pob item = (Pob) ptr;
    if (PyList_Check(ptr) ) {
        item = PyList_GetItem( (Pob) ptr,0);
    }
    if (PyTuple_Check(ptr)) {
        item = PyTuple_GetItem( (Pob) ptr, 0);
    }    
    if (PyDict_Check(ptr)) {
        item = PyGetDictItem( (Pob) ptr, 0);
    }

    if (PyLong_Check(ptr)) {
        return PyLong_AsLong(item);
    }
    if (PyFloat_Check(ptr)) {
        return PyFloat_AsDouble(item);
    }

    // For NumPy
    if (PyArray_Check(ptr)) {
        return NumpyGetNumericElemt( (Pob) ptr, 0);
    }
    return getNaN();

}

F64   GetNumericElement(const void* Y, I32 idx0) {
    //
    I32 idx = idx0; // zero-based
    Pob item = NULL;
    if (PyList_Check(Y)) {
        item = PyList_GetItem( (Pob) Y, idx);
    }
    if (PyTuple_Check(Y)) {
        item = PyTuple_GetItem( (Pob) Y, idx);
    }
    if (PyDict_Check(Y)) {
        item = PyGetDictItem( (Pob) Y, idx);
    }

    if (item) {
        if (PyLong_Check(item)) {
            return PyLong_AsLong(item);
        }
        if (PyFloat_Check(item)) {
            return PyFloat_AsDouble(item);
        }
    }
    
    // For NumPy
    if (PyArray_Check(Y)) {
        return NumpyGetNumericElemt( (Pob) Y, idx);
    }
    return getNaN();
}


void* GetData(const void* ptr) { 
    return PyArray_DATA( (Pob) ptr);
}

/*
int   GetDataType(VOID_PTR Y) {

    // For NUMPY ONLY
    if (!PyArray_Check(Y)) {
        return DATA_UNKNOWN;
    }

    int dtype = PyArray_TYPE(Y);
    if (dtype == NPY_INT16) {  return DATA_INT16; }
    if (dtype == NPY_INT32) {  return DATA_INT32; }
    if (dtype == NPY_INT64) { return DATA_INT64; }
    if (dtype == NPY_FLOAT) { return DATA_FLOAT; }
    if (dtype == NPY_DOUBLE) { return DATA_DOUBLE; }

    return DATA_UNKNOWN;
}
*/

int  GetDim1(const void* ptr) {
    if (PyArray_Check(ptr)) {
        npy_intp* dims = PyArray_DIMS(ptr);
        return (int) dims[0];
    }
    if (PyList_Check(ptr) || PyTuple_Check(ptr) ){
        return (int) PyObject_Size(ptr);
    }
    return -9999L;
}

int  GetDim2(const void* ptr) {

    if (!PyArray_Check(ptr)) {
        return -9999L;
    }
    npy_intp* dims = PyArray_DIMS(ptr);
    return dims[1];
}

int  GetNumOfDim(const void* ptr) {
    if (PyArray_Check(ptr)) {
        return  PyArray_NDIM(ptr);        
    }

    if (PyList_Check(ptr) || PyTuple_Check(ptr)) {
        return 1L;
    }       
    return -9999L;
}

void GetDimensions(const void* ptr, int dims[], int ndims) {
    if (PyArray_Check(ptr)) {
        int       N      = min(ndims, PyArray_NDIM(ptr));
        npy_intp* npdims = PyArray_DIMS(ptr);
        for (int i = 0; i < N; i++) {
            dims[i] = npdims[i];
        }
    }

    if (PyList_Check(ptr) || PyTuple_Check(ptr) ) {
        dims[0] = PyObject_Size(ptr);
    }   
}
 
void * SetDimensions(const void* ptr, int dims[], int ndims) {
    if (PyArray_Check(ptr)) {
        //PyObject* PyArray_Newshape(PyArrayObject * self, PyArray_Dims * newshape, NPY_ORDER order);
        //NPY_ANYORDER: keep the orignal order
        
        PyArray_Dims newshape;
        npy_intp  newdims[100];
        newshape.len = ndims;
        newshape.ptr = newdims;
         for (int i = 0; i < ndims; i++) {
             newdims[i] = dims[i];
        }

         PyObject *newptr=PyArray_Newshape(ptr, &newshape, NPY_ANYORDER);

         return newptr;

    }
    return NULL;
    if (PyList_Check(ptr) || PyTuple_Check(ptr)) {
      //  dims[0] = PyObject_Size(ptr);
    }
}

int  GetNumberOfElements(const void* ptr) {
    
    if (PyArray_Check(ptr)) {
        return PyArray_SIZE(ptr);
    }

    if (PyList_Check(ptr) || PyTuple_Check(ptr)) {
        return PyObject_Size(ptr);
    }

    if (PyUnicode_Check(ptr) ) {
    // A string has a dict but  is not treated as an object and instead as a 1 scalar
    // This must come before the "next" PyGetDict(ptr);"
        return 1;
    }

    PyObject* dict = PyGetDict(ptr);
    if (dict) {
        return PyObject_Size(dict);
    }

    if (PyLong_Check(ptr) || PyFloat_Check(ptr)) {
        return 1;
    }

    return 0;
}

I32  GetNumberOfFields(const void* structVar) {

    Pob dict = PyGetDict(structVar);

    if (dict) {
        return PyDict_Size(dict);
    }
    return -999;
}


int HaveEqualDimesions(const void* p1, const void* p2);
int CopyNumericObjToF32Arr(F32PTR outmem, VOID_PTR infield, int N);
 
void* CreateNumVar(DATA_TYPE dtype, int* dims, int ndims, VOIDPTR* data_ptr) {

    npy_intp dimsnp[10];
    for (int i = 0; i < ndims; i++) {
        dimsnp[i] = dims[i];
    }

    int dtypenp;
    if   (dtype == DATA_INT16)
        dtypenp = NPY_INT16;
    else if (dtype == DATA_INT32)
        dtypenp = NPY_INT32;
    else if (dtype == DATA_FLOAT)
        dtypenp = NPY_FLOAT;
    else if (dtype == DATA_DOUBLE)
        dtypenp = NPY_DOUBLE;
    else if (dtype == DATA_INT64)
        dtypenp = NPY_INT64;
    else {
        *data_ptr = NULL;
        return NULL;
    }

    PyObject* array = PyArray_NewFromDescr(
        &PyArray_Type, 
        PyArray_DescrFromType(dtypenp) /*reference stolen*/,
        ndims, dimsnp, NULL, NULL, NPY_ARRAY_F_CONTIGUOUS|NPY_ARRAY_WRITEABLE|NPY_ARRAY_ALIGNED, NULL); // n ew ref

    if (array && data_ptr) {
        *data_ptr = PyArray_DATA(array);
    }
 
    return array;
}

/*
Create a new ojbet

typeobjet.c: object_new->(_PyObject_InitializeDict(dictobject.c)->_PyObject_ComputedDictPointer

_PyObject_MakeTpCall->tp_call[instancemethod_call]->

type_call-.(type_new, type_init)

inst_dict is created during _PyObjectDict_SetItem for the type with Py_FLAGS_Heap_type set.

*/

static void* pyCreateStructObject() {
    PyObject* obj;
    /*    if (PyErr_Occurred()) {
        PyErr_PrintEx(0);
        PyErr_Clear(); // this will reset the error indicator so you can run Python code again
    }    */
   // PyObject* empTuple = PyTuple_New(0);
   // obj = PyObject_CallObject((PyObject*)&BarObject_Type, NULL/*no args*/); //new ref 
   // obj = PyObject_CallFunction((PyObject*)&BarObject_Type, NULL); //new ref
   // obj = PyObject_CallMethod(currentModule, "pyobject", NULL);
    obj = PyObject_CallObject(classOutput, NULL);
    //r_printf("pyCreateStructObject  %#x   \n", obj);
    return obj;
}
static void pySetAddField(void *obj, char *fname, void *value) {
    if (value == NULL) return ;

    PyObject_SetAttrString(obj, fname, value); // do not steal ref and value is increfed.
    Py_XDECREF(value);  // decref it so    
}

void* CreateStructVar(FIELD_ITEM* fieldList, int nfields) {

    // Find the sentinal element to update nfields
	int nfields_actual = 0;
	for (int i = 0; i < nfields; ++i) {
        nfields_actual++;
        if (fieldList[i].name[0] == 0) {
            nfields_actual--;
            break;
        }
	}
	nfields = nfields_actual;

	PyObject * _restrict out;
	{
		//char * fldNames[100];		
		//for (int i = 0; i < nfields; i++) { fldNames[i] = fieldList[i].name; }					
		//mwSize dims_2d[2] = { 1,1 };
		//out = mxCreateStructArray(2, dims_2d, nfields, fldNames);
        out = pyCreateStructObject(); 
	}

	for (int i = 0; i < nfields; i++){	 

		if (fieldList[i].ptr == NULL) continue;

        Pob optr = NULL;
		if (fieldList[i].type == DATA_STRUCT){
            optr = (void*)fieldList[i].ptr;            
        }  else {
            optr = CreateNumVar(fieldList[i].type, fieldList[i].dims, fieldList[i].ndim, fieldList[i].ptr);
        }
		
        if (optr) {
            pySetAddField(out, fieldList[i].name, optr);  //pyptr de-refed inside.
        }        
	}

	return (void*)out; 
}

/////////////////////////////////////////////////////
// Define a class type with a dict
/////////////////////////////////////////////////////

typedef struct {
    PyObject_HEAD
    PyObject* inst_dict;
    PyObject* weakreflist;
} BarObject;

static PyMemberDef BarObject_members[] = {
    {"__dict__", T_OBJECT_EX, offsetof(BarObject, inst_dict), 0,    "instance's dictonary"},
    {NULL}  /* Sentinel */
};

//https://stackoverflow.com/questions/1954494/python-instance-method-in-c
//https://nedbatchelder.com/text/whirlext.html

static PyObject * BarObject_new(PyTypeObject* type, PyObject* args, PyObject* kwds) {

    // !PyArg_ParseTupleAndKeywords(args, kwds, "|OOi", kwlist, &first, &last,   &self->number))

    BarObject* self = type->tp_alloc(type, 0);

    if (self) {
        self->inst_dict   = PyDict_New();
        self->weakreflist = Py_None;
       Py_INCREF(Py_None);
    }   

    // abc_ide_util_python.c:996 : 14 : warning : format ‘ % x’ expects argument of type ‘unsigned int’, but argument 3 has type ‘PyObject *
    //r_printf("New called...%#x  inst_dict %#x \n", self, self->inst_dict);
    r_printf("New called...%p  inst_dict %p \n", self, self->inst_dict);
    
    return self;
}
static int BarObject_init(BarObject* self, PyObject* args, PyObject* kwds) {

    static char* kwlist[] = { "first", "last", "number", NULL };
    PyObject* first = NULL, * last = NULL, * tmp;

    // !PyArg_ParseTupleAndKeywords(args, kwds, "|OOi", kwlist, &first, &last,   &self->number))

    //r_printf("Init called... %#x\n", self);
    r_printf("Init called... %p\n", self);
    return 0;
}

static void  BarObject_dealloc(BarObject* self) {
    self->weakreflist = NULL;
    Py_DECREF(Py_None);

    Py_XDECREF(self->inst_dict);
    Py_TYPE(self)->tp_free((PyObject*)self);        
}



// Using a variable name (like x) in the standard Python REPL is equivalent to print(repr(x))
// 
// https://github.com/python/mypy/issues/1325
// https://stackoverflow.com/questions/5061251/create-a-python-type-from-c-that-implements-a-dict
// https://docs.python.org/3/c-api/typeobj.html#typedef-examples
 PyTypeObject BarObject_Type = {
  PyObject_HEAD_INIT(NULL)
  .tp_new       = BarObject_new, //PyType_GenericNew,
  .tp_name      = "Rbeast.pyobject",
  .tp_base      = NULL,
  .tp_basicsize = sizeof(BarObject),
  .tp_getattro  = PyObject_GenericGetAttr,
  .tp_setattro  = PyObject_GenericSetAttr,
   //BarObject_Type.tp_getattr = PyObject_GenericGetAttr; //This field is deprecated
   //BarObject_Type.tp_setattr = PyObject_GenericSetAttr; //This field is deprecated  
  .tp_flags          = Py_TPFLAGS_DEFAULT, //Py_TPFLAGS_BASETYPE|Py_TPFLAGS_TYPE_SUBCLASS, //| Py_TPFLAGS_HAVE_GC | Py_TPFLAGS_DICT_SUBCLASS;
    //BarObject_Type.tp_itemsize = 0;
  .tp_dictoffset     = offsetof(BarObject, inst_dict),
  .tp_weaklistoffset = 0, //offsetof(BarObject, weakreflist),
  .tp_members        = BarObject_members,
  .tp_doc      = "Doc string for class Bar in module Foo.",
  .tp_init     = (initproc)BarObject_init,        
  .tp_dealloc  = (destructor)BarObject_dealloc,
};

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


I32  CheckInterrupt() {
// https:// stackoverflow.com/questions/14707049/allowing-ctrl-c-to-interrupt-a-python-c-extension/33652496#33652496
// https:// pybind11.readthedocs.io/en/stable/faq.html#how-can-i-properly-handle-ctrl-c-in-long-running-functions
// https:// docs.python.org/3/c-api/exceptions.html#c.PyErr_CheckSignals
    return PyErr_CheckSignals() != 0; 
    
}
void ConsumeInterruptSignal() { return; }

#endif

#include "abc_000_warning.h"