import os, sys, glob
from  setuptools import setup, find_packages, Extension, find_namespace_packages

#import numpy as np

# chichen-and-egg problems:  setup should install numpy but numpy is needed to get the numpy headers in np.get_include()
# and run the setup before it is installed
# 
# https://stackoverflow.com/questions/56008251/setup-requires-does-not-seem-to-install-dependencies
# https://stackoverflow.com/questions/35516059/python-is-not-installing-dependencies-listed-in-install-requires-of-setuptools
# https://stackoverflow.com/questions/19919905/how-to-bootstrap-numpy-installation-in-setup-py/21621689

# A few solutions: https://stackoverflow.com/questions/19919905/how-to-bootstrap-numpy-installation-in-setup-py/21621689
# This should now (since 2018-ish) be solved by adding numpy as a buildsystem dependency in pyproject.toml, 
# so that pip install makes numpy available before it runs setup.py.

# We fix the problem by using pyproject.toml
# Alternatively, we also included a solution from 
# https://stackoverflow.com/questions/19919905/how-to-bootstrap-numpy-installation-in-setup-py/21621689
# by using the deferring the import till it has been installed

# https://stackoverflow.com/questions/68282771/python-packaging-build-requirements-in-pyproject-toml-vs-setup-requires
class get_numpy_include(object):
    """Defer numpy.get_include() until after numpy is installed."""
    def __str__(self):
        from numpy import get_include 
        # this is the only function we need here:  not sure if this help reduce 
        # the possibility of "runtimeerror: module compiled against api version 0xe but this version of numpy is 0xd"
        return get_include()


# https://stackoverflow.com/questions/49640513/unicodedecodeerror-charmap-codec-cant-decode-byte-0x9d-in-position-x-charac
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname),encoding='utf8').read()
    #return open(os.path.join(os.path.dirname(__file__), fname)).read()


def is_platform_windows():
    return sys.platform == "win32" or sys.platform == "cygwin"

def is_platform_mac():
    return sys.platform == "darwin"
    

extralibs  = []
if is_platform_windows():
   extralibs = ["kernel32", "user32", "gdi32" ]

filenames    = glob.glob('ext_src/*.c');
filenames.append("ext_src/abc_ioFlushcpp.cpp")
      
modules = Extension(
            "Rbeast.Rbeast",
            sources       = filenames,
            #include_dirs  = [np.get_include(), "ext_src/"],   #Bad option bcz np needs to be installed first
            #include_dirs  = [get_numpy_include(), "ext_src/"],#Bad option bcz np needs to be installed first as build-time dependencies                  
            include_dirs  = ["ext_src/"],                 
            define_macros = [('P_RELEASE','1'),('R_INTERFACE','0')],
            libraries     = extralibs
        )
        
#packages = find_packages( include=['exampleproject','exampleproject.*','data'],  exclude=['figures', 'output', 'notebooks'])
#packages = find_packages( exclude=['figures', 'output', 'notebooks'])           
#packages = find_packages( exclude=['tests'])           
packages = find_packages( include = ['Rbeast'],exclude=['Rbeast.src'] )     
packages = find_namespace_packages( where='./py_src', exclude=['build','tests','extension_src'])
# print(packages)

# - this is the DISTRUBTION NAME not the IMPORT NAME
# PyPI requiers the normalized lowecase name. For some reason,
# cibuildwhell v2.xx always generates a uppercase name if name = "Rbeast"
# As as workaround, we explictly set the distrubtion name to the lowercase for CP2.7
pkgname = 'Rbeast'
if sys.version_info[0] == 3  and sys.version_info[1] == 7:
   pkgname = 'rbeast'

setup(

    #name             = "rbeast",   
    # - this is the DISTRUBTION NAME not the IMPORT NAME
    # - PyPI enforces PEP 503 normalization, which means it treats Rbeast and rbeast as the same project,
    # - but it requires the uploaded files to strictly follow the lowercase naming convention.
    name             = pkgname,    
    version          = '0.1.25',
    description      = "Bayesian changepoint detection and time series decomposition",
    author           = "Kaiguang Zhao",
    author_email     = 'zhao.1423@osu.edu',
    url              = 'https://github.com/zhaokg/Rbeast',    
    keywords         = ['changepoint', 'structural breaks','time series decomposition' 'time series analysis', 'trend analysis'],
    python_requires  = '>=3.5',  
    long_description = read('README.md'),
    long_description_content_type = "text/markdown",        
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",        
        "Programming Language :: Python :: 3.9",                
        "Programming Language :: Python :: 3.10",                        
        "Programming Language :: Python :: 3.11",   
        "Programming Language :: Python :: 3.12",   
        "Programming Language :: Python :: 3.13",   
        "Programming Language :: Python :: 3.14",   
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Hydrology",
        "Topic :: Scientific/Engineering :: Atmospheric Science",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
        "Topic :: Scientific/Engineering",
        "License :: OSI Approved :: MIT License",
    ],    
    #setup_requires      = ['numpy'],                                # Deprecated in favor of pyproject.toml
    #install_requires     =  ['numpy>=1.17.3', 'matplotlib>=2.2.0'], # ['numpy>=1.10', 'matplotlib>=2.2.0'],
    install_requires     =  ['numpy>=1.10.0'],                       # remove the depedence on matlplotlib; give a warning if matplot is missing
    #entry_points        ={  'console_scripts': ['mycommand=exampleproject.data:main1'] },    

    # ~ packages and package_dir arguments in the setup() function are what determine the import name.   
    packages             = packages,    
    package_dir          = {"Rbeast": "py_src/Rbeast","Rbeast.data": "py_src/Rbeast/data", '': '.'},      
    
    include_package_data = True,      
    package_data         = {'Rbeast': ['data/nile.csv'],'Rbeast.data': ['googletrend.csv'] ,'Rbeast': ['data/*.npy','data/*.txt','data/*.csv']},
    
    exclude_package_data = {'': ['*.c','*.cpp']},       
    ext_modules          = [modules]
)
