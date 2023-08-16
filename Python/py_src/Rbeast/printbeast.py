from . import Rbeast as cb

def print(o, index = 0):
    """
   USAGE: printbeast(o, index)

   o :  the time series analysis output from  beast123; o
   should contain results for mulltiple time series

   index: if o contains results for more than 1 time
   series, index specifies for which time series the result is printed.
   If o is the result for a 3D stacked cube, index will be a vector of 2
   integer to specify the row and col of the desired pixel. If o contains
   only one time series, index will be ignored

  Contact info : To report bug or get help, do not hesitate to contact Kaiguang Zhao
   at <strong>zhao.1423@osu.edu</strong>.
    """
    if hasattr(o, 'marg_lik'):
        if isinstance(index, int):
           index=index+1
        else:
           index = [x+1 for x in index]   
        cb.Rbeast('print',o,index);

def obj_repr(o):    
    N    = len( o.__dict__.keys())
    s    = "Object of %d field(s):\n\n" % (N)
    s    = my_own_repr(o, s ,0+1)
    return s
   
import numpy as np   
def my_own_repr(o, s, left):

    dict           = list( o.__dict__.keys())
    maxFldLen      = max([len(x) for x in dict])
    N              = len(dict)
    for i in range( N ):
        key = dict[i]
        s = s +'%*s' % (left,'')+ "%-*.*s : " % (maxFldLen,maxFldLen, key)
        e = o.__dict__[key]
        if isinstance(e, np.ndarray):
            etype = e.dtype
            if e.size == 1:
                if etype == np.int16 or etype== np.int32 or etype== np.int64:
                   s = s + "%d \n" % (e[0])
                elif etype == np.float32 or etype== np.float64:
                   s = s + "%g \n" % (e[0])
                else:
                   s = s + "%g \n" % (e[0])
                continue

            if  etype == np.int16:
                typestr ='int16'
            elif etype == np.int32:
                typestr = 'int32'
            elif etype == np.int64:
                typestr = 'int64'
            elif etype == np.float32:
                typestr = 'float32'
            elif etype== np.float64:
                typestr = 'float64'
            else:
                typestr = 'others'

            dims = e.shape
            tmp = "[%d" % ( dims[0])
            for j in range(len(dims)-1):
                tmp = tmp + "x%d" % (dims[j+1])
            s = s +tmp +" " + typestr + "] \n"
        #if isinstance(e, BeastOutputClass):
        if isinstance(e, type(o)):        
            s = s+ "[ 1 object with %d fields] \n" % (len( e.__dict__.keys()))
            s = my_own_repr(e,s, left+ maxFldLen+3)
        if isinstance(e, str):
            s = s+ "'"+e+ "'\n"
    return s


#cbeast.pyobject.__getitem__ = extract

