#https://gist.github.com/CMCDragonkai/510ce9456a0429f616baa243d1de3dbf
#https://stackoverflow.com/questions/14295277/not-exporting-functions-from-python-module

#Any directory with __init__.py is considered a package in Python.
#
#Any python files inside a package is considered a module.
#
#Modules contain functions and other bindings that is always exported.
#
#If you are outside the package, and you want to import a module from a package:
#
#from package import module
#module.use_function()

#However this is only because the __init__.py is empty. We can allow the package to export functions and modules as well. Which means this should be possible:

#import package
#package.module.use_function()
#package.use_function()
#To do this, the package's __init__.py must contain something like:

#from . import module
#from .mod import *
#Therefore the __init__.py is kind of like an exporting script for the package. Similar to how I use index.js in my JavaScript #packages. This means you can use the package __init__.py as a sort of staging point for all exports. And anything that isn't# #exported here is hidden. Unless the user really wants to acquire it, in which case they can explicitly use the module.

#To actually have encapsulation at the module level. You need to use _ prefix on your bindings. This is different from __ prefix #used in classes. Note that using the * import will ignore module bindings that have _ prefixed. But they can still be accessed #explicitly if the module is directly accessed. See: https://stackoverflow.com/a/1547160/582917

#I think this form of using __init__.py is the best way. Users shouldn't need to use the from ... import syntax unless they need #to import specific bindings and don't want to use qualified modules.


from .beast        import beast        as beast
from .beast_irreg  import beast_irreg  as beast_irreg
from .beast123     import beast123     as beast123
from .plotbeast    import plot         as plot
from .extractbeast import extract      as extract
from .printbeast   import print        as print, obj_repr
from .load_example import load_example as load_example
from .             import Rbeast       as cbeast
 

# https://stackoverflow.com/questions/2356399/tell-if-python-is-in-interactive-mode
def isInteractiveMode():
    import sys 
    return hasattr(sys, 'ps1') 


# https://stackoverflow.com/questions/67631/how-do-i-import-a-module-given-the-full-path
# https://docs.python.org/3/library/importlib.html#importing-a-source-file-directly
# https://stackoverflow.com/questions/10675054/how-to-import-a-module-in-python-with-importlib-import-module
# https://stackoverflow.com/questions/63156606/aliases-for-imported-modules-with-importlib
def LoadPlot():
    import sys,importlib, os
    mod       = sys.modules[__name__]    
    mod_path  = os.path.join(os.path.dirname(__file__), 'plotbeast.py')
    spec      = importlib.util.spec_from_file_location("plotbeast", mod_path)
    plotbeast = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(plotbeast)
    #sys.modules["module.name"] = foo    
    setattr(mod,'plot',  getattr(plotbeast, 'plot') )

# https://stackoverflow.com/questions/67631/how-do-i-import-a-module-given-the-full-path
# https://docs.python.org/3/library/importlib.html#importing-a-source-file-directly
# https://stackoverflow.com/questions/10675054/how-to-import-a-module-in-python-with-importlib-import-module
# https://stackoverflow.com/questions/63156606/aliases-for-imported-modules-with-importlib
def LoadPlot_v2():
    import sys, importlib
    mod       = sys.modules[__name__]    
    plotbeast = importlib.import_module('.plotbeast',__name__)
    setattr(mod,'plot',  getattr(plotbeast, 'plot') ) 

#if isInteractiveMode:
##     LoadPlot()
#      LoadPlot_v2 # this works 


# https://stackoverflow.com/questions/1957780/how-to-override-the-operator-in-python
# https://stackoverflow.com/questions/57413453/is-it-possible-to-override-getitem-at-instance-level-in-python     

class BeastOutputClass:
    pass
BeastOutputClass.__getitem__ = extract    
BeastOutputClass.__repr__    = obj_repr

class args:
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self,key, value)
    pass


def RemoveUselessMembers():
    import sys    
    dict = sys.modules[__name__].__dict__.copy()    
    mod  = sys.modules[__name__]
    membersKept = ['args','obj_repr','beast','print','beast_irreg', 'beast123', 'plot','cbeast', '__dict__','BeastOutputClass','load_example']
    for x in dict:
        if x not in membersKept:
            delattr(mod, x)
            
#RemoveUselessMembers()

cbeast.setClassObjects(BeastOutputClass)



