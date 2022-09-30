from . import Rbeast as cb
from numpy import ndarray, squeeze

def beast123(Y, metadata=None, prior=None, mcmc=None, extra=None ):
     
      isNumpyInput = False;
      if   isinstance(Y, list) or isinstance(Y, tuple):
            isNumpyInput = False
      elif isinstance(Y, ndarray):
            isNumpyInput = True
      elif hasattr(Y,'to_numpy'):
            Y = getattr(Y,'to_numpy')()
            isNumpyInput = True
      else:
            raise ValueError('Unknown formats for the input Y.')
      
      if isNumpyInput:
            Y=squeeze(Y)
   
      o=cb.Rbeast('beastv4',Y, metadata, prior, mcmc, extra)
      return (o)


