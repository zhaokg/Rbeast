import numpy as np
import os


def load_example(name):
    name = name.lower()    
    if not (name == 'nile' or name == 'beach' or name == 'ndvi'):
        return None

    try:
        import pandas as pd            
        data_path = os.path.join(os.path.dirname(__file__), 'data', name + '.csv')
        if (name == 'nile'):         
            data = pd.read_csv(data_path)            
            return data.flow, data.year
        if (name == 'beach'):         
            data = pd.read_csv(data_path)            
            return data.y, data.time            
    except Exception as error:        
        if (name == 'nile'):
            x   = np.genfromtxt(data_path, delimiter=',', skip_header=1)
            #data = {'year': x[:, 0], 'flow': x[:, 1]}
            return x[:,0], x[:,1]
        if (name == 'beach'):
            x = np.genfromtxt(data_path, delimiter=',', skip_header=1)
            #data = {'time': x[:, 0], 'y': x[:, 1]}
            return x[:,1], x[:,0]
        
    if (name == 'ndvi'):        
        f1   = os.path.join(os.path.dirname(__file__), 'data',  'ndvi.npy')
        ndvi = np.load(f1)
        f1 = os.path.join(os.path.dirname(__file__), 'data', 'ndvi_fyear.npy')
        year = np.load(f1)
        f1 = os.path.join(os.path.dirname(__file__), 'data', 'ndvi_filename.txt')
        with open(f1) as f:
            datestr = [line.rstrip() for line in f]
        return ndvi, year, datestr

