import os,sys
import numpy as np
from   setuptools import setup, find_packages, Extension,find_namespace_packages
import glob


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

def is_platform_windows():
    return sys.platform == "win32" or sys.platform == "cygwin"

def is_platform_mac():
    return sys.platform == "darwin"
    
filenames    = glob.glob('ext_src/*.c');
filenames.append("ext_src/abc_ioFlushcpp.cpp")

extralibs  = []
if is_platform_windows():
   extralibs=["kernel32", "user32", "gdi32" ]
      
modules = Extension(
            "Rbeast.Rbeast",
            sources       = filenames,
            include_dirs  = [np.get_include(), "ext_src/"],
            define_macros = [('P_RELEASE','1'),('R_INTERFACE','0')],
            libraries     = extralibs
        )
        
#packages = find_packages( include=['exampleproject','exampleproject.*','data'],  exclude=['figures', 'output', 'notebooks'])
#packages = find_packages( exclude=['figures', 'output', 'notebooks'])           
#packages = find_packages( exclude=['tests'])           
packages = find_packages( include =['Rbeast'],exclude=['Rbeast.src'])     
packages = find_namespace_packages(where='./py_src', exclude=['build','tests','extension_src'])
#print(packages)
#print(packages)
setup(
    name             = "Rbeast",   
    version          = '0.1.8',
    description      = "Python package for Bayesian changepoint detection and time series decomposition",
    author           = "Kaiguang Zhao",
    author_email     = 'zhao.1423@osu.edu',
    url              = 'https://github.com/zhaokg/Rbeast',    
    keywords         = ('changepoint', 'structural breaks','time series decomposition' 'time series analysis', 'trend analysis',),
    python_requires  = '>=3.7',  
    long_description = read('README.md'),
    long_description_content_type = "text/markdown",        
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",        
        "Programming Language :: Python :: 3.9",                
        "Programming Language :: Python :: 3.10",                        
        "Programming Language :: Python :: 3.11",   
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Hydrology",
        "Topic :: Scientific/Engineering :: Atmospheric Science",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
        "Topic :: Scientific/Engineering",
        "License :: OSI Approved :: MIT License",
    ],    
    install_requires=['numpy', 'matplotlib>=2.2.0'],
    #entry_points={  'console_scripts': ['mycommand=exampleproject.data:main1'] },    
    packages             = packages,    
    package_dir          = {"Rbeast": "py_src/Rbeast","Rbeast.data": "py_src/Rbeast/data", '': '.'},      
    include_package_data = True,      
    package_data         = {'Rbeast': ['data/nile.csv'],'Rbeast.data': ['beach.csv'] ,'Rbeast': ['data/*.npy','data/*.txt']},
    exclude_package_data = {'': ['*.c','*.cpp']},       
    ext_modules          = [modules]
)
