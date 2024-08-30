##    Installation for Matlab or Octave:
 
*  **Method 1**: Manually download the files in this folder (e.g, beast.m and testdata) to your local folder
             Then, add your local folder to Matlab's search path by running `addpath( genpath('c:\beast') )`.
             Make sure to replace `c:\beast` with your own path.
 
* **Method 2**:  The recommended way to install is to simply run the following line of command in Matlab's console:
 ```
   % (a) Install to a temp folder
   eval(webread('http://b.link/rbeast',weboptions('cert','')))  

   % (b) Or install it to a chosen folder by first setting the beastPath variable
   beastPath = 'c:\beast'                 
   eval(webread('http://b.link/rbeast',weboptions('cert','')))  % install to the folder specified by beastPath    
```
   
##    Uninstall

* **METHOD 1**: Manually delete all the downloaded files 
* **METHOD 2**: Automatically remove the files by running `rbeast_uninstall`
 
 
##      Usage and example

#### Run 'help beast', 'help beast123', or 'help beast_irreg' for usage and examples  
 
#### Major functions available:

* beast:            handle a single regular time series
* beast_irreg:      handle a single irregular time series
* beast123:         handle one or more time seires or stacked images 
* rbeast_uninstall: remove the installed files from the harddisk
* rbeast_update:    check github and update to the latest BEAST version\n");

#### Examples
```
    load('Nile.mat')             % Nile river annual streamflow: trend-only data
    o=beast(Nile, 'start', 1871, 'season','none') 
    printbeast(o))
    plotbeast(o))
    help beast
    help beast123
```

      


##    Compilation of the mex binary from the C/C++ source files

The BEAST program includes a Matlab mex library compiled from the C soure code (e.g., Rbeast.mexw64 for Windows, Rbeast.mexa64 for Linux, Rbeast.mexmaci64 for MacOS) and some Matlab wrapper functions (e.g.,beast.m, and beast123.m) similar to the R interface, plus several test datasets (e.g., Nile.mat, and co2.mat). We generated the Matlab mex binary library for Win10, Ubuntu 22.04, and macOS High Sierra 10.13. If they fail on your machine, the mex library can be compiled from the C source code files under Rbeast\Source. If needed, we are happy to work with you to compile for your specific OS or machines. Additional information on compilations from the C source is given below.

### Step 1. Download the C source files
* **Method 1**: Manually download all the C source files from the `Source` folder at https://github.com/zhaokg/Rbeast to your local folder.
* **Method 2**: After installing Rbeast in Matlab, run `rbeast_src_download` to automatically download the source files to your beast folder

### Step 2. Install a C/C++ compiler (e.g., gcc and clang)
Many choices are possible; given below are the options I used.
* **Window**: CRAN's distribution of [MinGW gcc](https://cran.r-project.org/bin/windows/Rtools/rtools43/files/rtools43-5863-5818.exe), with more information available [here](https://cran.r-project.org/bin/windows/Rtools/rtools43/rtools.html). By default, gcc is installed as `‪C:\rtools43\mingw64\bin\gcc.exe`.
  
* **Linux/Ubuntu**: Linux/Ubuntu should come with gcc pre-installed. If not, use `sudo apt install build-essential`.
* **MacOS**       : Install xcode to get the clang compiler. Alternatively, use `brew` to install, [following](https://stackoverflow.com/questions/8674546/how-to-update-llvm-and-clang-on-mac-os-x)
  ```
  # Install Homebrew
  /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
  
  # Use brew to install clang
  brew install llvm  
  ```
### Step 3a. Use Matlab's mex to compile
If the mex command is apropriately set up in Matlab, here are two ways to build the Rbeast.mex binary file:
* **Method 1**: Run `rbeast_src_download` first to download the C files, and then run `rbeast_src_compile` to compile using `mex`.
* **Method 2**: Call `mex` manually in Matlab using the following command, assuming that Matlab's current directory is the source folder:
    - Linux/MacOS
      ```
      mex -v CFLAGS='-pthread -DM_RELEASE -UUSE_MEX_CMD -fPIC -O2 -Wall -std=gnu99 -mfpmath=sse -msse2 -mstackrealign' -lmwservices -lut  -lpthread *.c -output Rbeast
      ```
    - Windows
      ```
      mex -v CFLAGS='-pthread -DM_RELEASE -UUSE_MEX_CMD -fPIC -O2 -Wall -std=gnu99 -mfpmath=sse -msse2 -mstackrealign' -lmwservices -lut  -lpthread  -lkernel32 -lgdi32 -luser32  *.c -output Rbeast
      ```
### Step 3b. Build the mex binary using gcc or clang
If your matlab's mex doesn't work or Step3a fails for some reason, we can circument mex to explicity build the binary using gcc and clang.  Matlab's mex file is nothing but 
a dynamic libray (e.g., a .dll file on Windows, and shared object file on Linux). 

THe key for successful compilation is to locate the Matlab's paths of include headers and libraries (e.g., import libs for Windows and run-time libraries for non-Windows systems). Here are  typical paths using some machines as examples. 

* macOS-intel：
   - include path = `'/Applications/MATLAB_R2019b.app/extern/include'`
   - library path =`/Applications/MATLAB_R2019b.app/bin/maci64/`.
* macOS-apple silicon：
   - include path = `'/Applications/MATLAB_R2023b.app/extern/include/'`
   - library path =`/Applications/MATLAB_R2023b.app/bin/maca64/`.
* Linux/Ubunbu：
   - include path = `'/MATLAB/extern/include`'
   - library path =`/MATLAB/bin/glnxa6`.
* Windows for Visual Studio：
  - include path = `'C:/Program Files/MATLAB/R2019a/extern/include'`
  - library path =`C:\Program Files\MATLAB\R2019a\extern\lib\win64\microsoft`.
* Windows for gcc：
  - include path = `'C:/Program Files/MATLAB/R2019a/extern/include'`
  - library path =`C:\Program Files\MATLAB\R2019a\extern\lib\win64\mingw64`.
* Windows/Octave：
   - include path = `'C:\Program Files\GNU Octave\Octave-8.2.0\mingw64\include\octave-8.2.0\octave'`
   - library path =`C:\PROGRA~1\GNUOCT~1\OCTAVE~1.0\mingw64\lib\octave\8.2.0`.

Suppose that the current working dirctory in the OS terminal/shell is the source folder where the beast c files are saved. Here are some examples of the commands to build Rbeast.mex：

* MacOS-Apple Silcon:   `clang  -mmacosx-version-min=10.13 -dynamiclib  -fPIC -I/usr/local/include -I/Applications/MATLAB_R2023b.app/extern/include/ -DM_RELEASE -Wall -O2  -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress -L/Library/Frameworks/R.framework/Resources/lib -L/usr/local/lib -L/Applications/MATLAB_R2023b.app/bin/maca64/ -msse2 -mstackrealign -lpthread -lm -lut -lmwservices -lmat -lmex -lmx *.c -o Rbeast.mexmaca64`
  
* MacOS-Intel:   `clang  -mmacosx-version-min=10.13 -dynamiclib  -fPIC -I/usr/local/include -I/Applications/MATLAB_R2019b.app/extern/include/ -DM_RELEASE -Wall -O2  -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress -L/Library/Frameworks/R.framework/Resources/lib -L/usr/local/lib -L/Applications/MATLAB_R2019b.app/bin/maci64/ -msse2 -mstackrealign -lpthread -lm -lut -lmwservices -lmat -lmex -lmx *.c -o Rbeast.mexmaci64`

* Windows for Ocatve `gcc -shared -fPIC -pthread -O2 -DM_RELEASE -DO_INTERFACE -I"C:\Program Files\GNU Octave\Octave-8.2.0\mingw64\include\octave-8.2.0\octave"    *.c -LC:\PROGRA~1\GNUOCT~1\OCTAVE~1.0\mingw64\lib\octave\8.2.0  -Wl,--export-all-symbols   -loctinterp -loctave -lkernel32 -lgdi32 -luser32  -o Rbeast.mex`

* Windows for Matlab `gcc -shared -fPIC -pthread -O2 -DM_RELEASE -DO_INTERFACE -I"C:/Program Files/MATLAB/R2019a/extern/include"  *.c -LC:\Program Files\MATLAB\R2019a\extern\lib\win64\mingw64  -Wl,--export-all-symbols  -lm -lut -lmwservices -lmat -lmex -lmx -lkernel32 -lgdi32 -luser32  -o Rbeast.mexw64`

* Linux/ubuntu `gcc -shared -fPIC -fPIC -O2 -Wall -std=gnu99 -mfpmath=sse -msse2 -mstackrealign -O2 -DM_RELEASE -DO_INTERFACE -I/MATLAB/extern/include  -L/MATLAB/bin/glnxa6  -Wl,--export-all-symbols  -lm -lut -lmwservices -lmat -lmex -lmx  *.c -o Rbeast.mexa64`

Make sure to use your own Matlab include and library paths. If the Rbeast C/C++ source folder is not the current working direcotry, replace `*.c` in th above with `YOUR_BEAST_SOURCE_PATH/*c`.

 
 