import os, importlib
from numpy import genfromtxt as np_genfromtxt,  load as np_load



def load_example(name):
    """
    Possible example datasets are:
      nile       : annual streamflow of the River Nile from 1871
      googletrend: monthly Google search popularity of the keyword 'beach'
      yellowstone: biweekly NDVI (greenness) of a Yellowstone site
      ohio       : satelite reflectanes and NDVI of an Ohio site
      co2        : monthly CO2 concentrations
      covid19    : weekly  new cases and deaths from COVID-19
      simdata    : a 774x3 matrix: three time series of length 774 
      imagestack : stacked images of NDVI over a small area
    """
    name = name.lower()    
    
    if not (name == 'nile' or name == 'googletrend' or name == 'imagestack' or name == 'ohio'  or 
       name == 'yellowstone' or name == 'co2' or name == 'covid19' or name == 'simdata'):
        print("Possible example datasets are:")    
        print("  nile       : annual streamflow of the River Nile from 1871")    
        print("  googletrend: monthly Google search popularity of the keyword 'beach'")    
        print("  yellowstone: biweekly NDVI (greenness) of a Yellowstone site' ")    
        print("  ohio       : satelite reflectanes and NDVI of an Ohio site  ")    
        print("  co2        : monthly CO2 concentrations  ")    
        print("  covid19    : weekly  new cases and deaths from COVID-19  ")   
        print("  simdata    : a 774x3 matrix: three time series of length 774 ")          
        print("  imagestack : stacked images of NDVI over a small area  \n\n")    
        return None
        
    data_path_csv = os.path.join(os.path.dirname(__file__), 'data', name + '.csv')
    data_path_npy = os.path.join(os.path.dirname(__file__), 'data', name + '.npy')
  
    #https://stackoverflow.com/questions/301134/how-can-i-import-a-module-dynamically-given-its-name-as-string    
    #the code below is not used any longer but kept here as a referenece
    
    # spec = importlib.util.find_spec('pandas')
    # if (spec is None   ):
    #    printf("Error; the 'pandas' package is needed to load the example dataset!\n")            
    # pd = importlib.import_module('pandas')

    try:
        import pandas as pd
        if (name == 'nile'):         
            data = pd.read_csv(data_path_csv)            
            return data.flow, data.year
        if (name == 'googletrend'):         
            data = pd.read_csv(data_path_csv)            
            return data.y, data.time 
        if (name == 'yellowstone'):         
            data = pd.read_csv(data_path_csv)            
            return data.ndvi, data.date
        if (name == 'co2'):         
            data = pd.read_csv(data_path_csv)            
            return data.co2, data.time             
        if (name == 'ohio' or name=='covid19'):         
            data = pd.read_csv(data_path_csv)            
            return data              
    except Exception as error:
        if (name == 'nile' or name == 'googletrend' or name == 'yellowstone' or name == 'co2'):
            x   = np_genfromtxt(data_path_csv, delimiter=',', skip_header=1)
            #data = {'year': x[:, 0], 'flow': x[:, 1]}
            return x[:,1], x[:,0]
            
    if (name == 'simdata'): 
    # np.save('y:/ndvi3d.npy',np.array([a,b,e],dtype='object'),allow_pickle=True)
        x   = np_genfromtxt(data_path_csv, delimiter=',', skip_header=1)
        return x  
        
    if (name == 'imagestack'): 
    # np.save('y:/ndvi3d.npy',np.array([a,b,e],dtype='object'),allow_pickle=True)
        data = np_load(data_path_npy,allow_pickle=True) 
        return data[0], data[1], data[2]   
   
    return None
   
   # the code below is not used any more.
    if (name == 'ndvi'):    
        f1   = os.path.join(os.path.dirname(__file__), 'data',  'ndvi.npy')
        ndvi = np_load(f1)
        f1 = os.path.join(os.path.dirname(__file__), 'data', 'ndvi_fyear.npy')
        year = np_load(f1)
        f1 = os.path.join(os.path.dirname(__file__), 'data', 'ndvi_filename.txt')
        with open(f1) as f:
            datestr = [line.rstrip() for line in f]
        return ndvi, year, datestr

