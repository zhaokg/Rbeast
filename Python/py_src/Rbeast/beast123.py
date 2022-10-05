from . import Rbeast as cb
from .cvt_to_numpy import force_convert_to_numpy
#Y = force_convert_to_numpy(Y)

def beast123(Y, metadata=None, prior=None, mcmc=None, extra=None ):
     
      Y = force_convert_to_numpy(Y)
      
      if hasattr(metadata, 'time'):
            time = metadata.time
            if hasattr(time, "year"):
                time.year = force_convert_to_numpy(time.year)
                if hasattr(time, "month"):
                    time.month = force_convert_to_numpy(time.month)
                if hasattr(time, "day"):
                    time.day   = force_convert_to_numpy(time.day)   
                if hasattr(time, "doy"):
                    time.doy   = force_convert_to_numpy(time.doy)   
            elif hasattr(time,'datestr') or hasattr(time,'dateStr'):
                pass
            else: 
                #then, we assume time is a numeric vector
                time = force_convert_to_numpy(time) 
                
            metadata.time = time
            
      o = cb.Rbeast('beastv4',Y, metadata, prior, mcmc, extra)
      return (o)


