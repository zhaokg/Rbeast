from numpy import ndarray, squeeze, array as nparray
from array import array as arrtype

IsTpye = isinstance

def force_convert_to_numpy(x): 
      if  IsTpye(x, list) or IsTpye(x, tuple) or IsTpye(x, arrtype):
           if not IsTpye(x[0], str):    # a sloppy way to test if the array is a string array or not
                 # this is not really needed here bcz they can be converted in the C code
                 return squeeze(nparray(x))
           else:
                 return x   # we assume if the first elm is a sting and the whole array is a string array,
      elif IsTpye(x, ndarray):
           return squeeze( x )
      elif hasattr(x,'to_numpy') and callable(getattr(x,'to_numpy')):
           x = getattr(x,'to_numpy')()
           return squeeze(x)
      else:
          print("WARNING: the data has an uknonwn format but will be forced into numpy.ndarray. If it fails, pls explicilty convert your data into numpy arrays.");
          #raise ValueError('Unknown formats for the input data.')
          return x

