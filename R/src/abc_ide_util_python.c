#include <math.h>
#include <string.h>

#include "assert.h"
#include "abc_000_warning.h"

#include "abc_ide_util.h"
#include "abc_common.h"
#include "abc_date.h"

#if P_INTERFACE==1  

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
       
    char tmpName[100 + 1];
    PyObject* keys = PyDict_Keys(ptr);              // new ref
    int       n    = PyList_Size(keys);
    item = NULL;
    for (int i = 0; i < n; i++) {
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

   
    char tmpName[100 + 1];
    PyObject* keys = PyDict_Keys(ptr);              // new ref
    int       n    = PyList_Size(keys);
    item = NULL;
    for (int i = 0; i < n; i++) {
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

    if (ptr == NULL) return 0;

    if (PyUnicode_Check(ptr)) {
        return 1;
    } 
    if (PyArray_Check(ptr) && PyArray_TYPE(ptr) == NPY_STRING) {
        return 1;
    }
    if (PyArray_Check(ptr) && PyArray_TYPE(ptr) == NPY_UNICODE) {
        return 1;
    }

    PyGetItem pyGetItem = NULL;
    if      (PyList_Check(ptr))        pyGetItem = PyList_GetItem;
    else if (PyTuple_Check(ptr))       pyGetItem = PyDict_GetItem;
 
    
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
    if (ptr == NULL) return 0;

    if (PyList_Check(ptr) || PyDict_Check(ptr) || PyTuple_Check(ptr) ) {
        return 1L;
    }
    
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
 
    if (PyFloat_Check(ptr)) {
        return 1L;
    } 
    if (PyArray_Check(ptr) && PyArray_TYPE(ptr) == NPY_FLOAT64) {
        return 1L;
    }

    PyGetItem pyGetItem = NULL;
    if (PyList_Check(ptr))           pyGetItem = PyList_GetItem;
    else if (PyTuple_Check(ptr))     pyGetItem = PyDict_GetItem;
    if (pyGetItem) {
        int sz = PyObject_Size(ptr);
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
    if      (PyList_Check(ptr))      pyGetItem = PyList_GetItem;
    else if (PyTuple_Check(ptr))     pyGetItem = PyDict_GetItem;
    if (pyGetItem) {
        int sz = PyObject_Size(ptr);
        int ok = 1;
        for (int i = 0; i < sz; ++i) {
            if (!PyLong_Check(pyGetItem(ptr, i))) {
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
    if      (PyList_Check(ptr))      pyGetItem = PyList_GetItem;
    else if (PyTuple_Check(ptr))     pyGetItem = PyDict_GetItem;
    if (pyGetItem) {
        int sz = PyObject_Size(ptr);
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
        PyObject* keys = PyDict_Keys(dict);              // new ref
        int       n    = PyList_Size(keys);

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

void* CreateStructVar(FIELD_ITEM* fieldList, int nfields); //Done

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
        char* str = PyUnicode_AsUTF8AndSize(ptr, &len);
        strncpy(dst, str, n);
        return len;
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
        char* str = PyUnicode_AsUTF8AndSize(tmpItem, &len);
        len = min(len,n - 1);
        memcpy(dst, str, len);
        dst[len] = 0;
        return len;
    }

    //If a Numpy string array
    if (PyArray_Check(ptr) && PyArray_TYPE(ptr)==NPY_STRING) {
        char* base   = PyArray_DATA(ptr);
        int   elsize = PyArray_ITEMSIZE(ptr);
        char* str    = base + elsize * idx;
        len     = min(elsize, n-1);
        memcpy(dst, str, len);
        dst[len] = 0;
        return len;
    }

    //If a Numpy Unicode array
    if (PyArray_Check(ptr) && PyArray_TYPE(ptr) == NPY_UNICODE) {
        char* base   = PyArray_DATA(ptr);
        int   elsize = PyArray_ITEMSIZE(ptr);
        char* str = base + elsize * idx;
        len = min(elsize/4, n - 1);  // Assuming a Unicode char has 4 bytse
        for (int i = 0; i < len; i++) {
            dst[i] = str[i * 4];
        }
        dst[len] = 0;
        return len;
    }

    //If a Numpy object array
    if (PyArray_Check(ptr) && PyArray_TYPE(ptr) == NPY_OBJECT) {
        void** base  = PyArray_DATA(ptr);
        int    elsize = PyArray_ITEMSIZE(ptr);

        void* tmpItem = base[idx];

        // It is a list or dict item
        if (tmpItem && PyUnicode_Check(tmpItem)) {
            char* str = PyUnicode_AsUTF8AndSize(tmpItem, &len);
            len = min(len, n - 1);
            memcpy(dst, str, len);
            dst[len] = 0;
            return len;
        }
 
    }

    dst[len] = 0;
    return 0;
}


void* GetField123(const void* structVar, char* fname, int nPartial) {

    Pob dict = NULL;

    if (PyDict_Check((void*)structVar)) {
        dict = structVar;
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
    
    Pob item = ptr;
    if (PyList_Check(ptr) ) {
        item = PyList_GetItem(ptr,0);
    }
    if (PyTuple_Check(ptr)) {
        item = PyTuple_GetItem(ptr, 0);
    }    
    if (PyDict_Check(ptr)) {
        item = PyGetDictItem(ptr, 0);
    }

    if (PyLong_Check(ptr)) {
        return PyLong_AsLong(item);
    }
    if (PyFloat_Check(ptr)) {
        return PyFloat_AsDouble(item);
    }

    // For NumPy
    if (PyArray_Check(ptr)) {
        return NumpyGetNumericElemt(ptr, 0);
    }
    return getNaN();

}
F64   GetNumericElement(const void* Y, I32 idx0) {
    //
    I32 idx = idx0; // zero-based
    Pob item = NULL;
    if (PyList_Check(Y)) {
        item = PyList_GetItem(Y, idx);
    }
    if (PyTuple_Check(Y)) {
        item = PyTuple_GetItem(Y, idx);
    }
    if (PyDict_Check(Y)) {
        item = PyGetDictItem(Y, idx);
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
        return NumpyGetNumericElemt(Y, idx);
    }
    return getNaN();
}

void* GetData(const void* ptr) { 
    return PyArray_DATA(ptr);
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
        return dims[0];
    }
    if (PyList_Check(ptr) || PyTuple_Check(ptr) ){
        return PyObject_Size(ptr);
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
// https://stackoverflow.com/questions/14707049/allowing-ctrl-c-to-interrupt-a-python-c-extension/33652496#33652496
// https://pybind11.readthedocs.io/en/stable/faq.html#how-can-i-properly-handle-ctrl-c-in-long-running-functions
// https://docs.python.org/3/c-api/exceptions.html#c.PyErr_CheckSignals
    return PyErr_CheckSignals() != 0; 
    
}
void ConsumeInterruptSignal() { return; }

#endif

#include "abc_000_warning.h"