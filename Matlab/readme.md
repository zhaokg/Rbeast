## Installation for MATLAB or Octave

### Method 1: Manual installation

Manually download the files in this Github folder (e.g., the m scripts and test data) to a local folder. Then add that folder to the MATLAB search path:

```matlab
% Replace `C:\beast` with the path to your local installation folder

addpath(genpath('C:\beast')) 
```



### Method 2: Automatic installation

The recommended installation method is to run the following command in the MATLAB console:

```matlab
% Install to a temporary folder
eval(webread('http://b.link/rbeast', weboptions('cert', '')))
```

To install to a specific folder, define `beastPath` before running the installer command:

```matlab
beastPath = 'C:\beast';                                        % Set your target folder
eval(webread('http://b.link/rbeast', weboptions('cert', '')))  % Install to beastPath
```

---

## Uninstall

- **Method 1:** Manually delete the downloaded files.
- **Method 2:** Run the following command in MATLAB:

```matlab
rbeast_uninstall
```

---

## Usage and examples

Run one of the following commands in MATLAB or Octave for usage instructions and examples:

```matlab
help beast
help beast123
help beast_irreg
```

### Major functions

- `beast`: Handles a single regular time series.
- `beast_irreg`: Handles a single irregular time series.
- `beast123`: Handles one or more time series or stacked images.
- `rbeast_uninstall`: Removes the installed files from the hard disk.
- `rbeast_update`: Checks GitHub and updates Rbeast to the latest version.

### Example

```matlab
load('Nile.mat')                         % Nile River annual streamflow: trend-only data
o = beast(Nile, 'start', 1871, 'season', 'none');
printbeast(o)
plotbeast(o)
help beast
help beast123
```

---

## Compiling the MEX binary from the C/C++ source files

Rbeast includes a MATLAB MEX library, such as `Rbeast.mexw64` for Windows, `Rbeast.mexa64` for Linux, `Rbeast.mexmaci64` for Intel macOS, and `Rbeast.mexmaca64` for Apple Silicon macOS, which is compiled from the C source code. The package also includes MATLAB wrapper functions such as `beast.m` and `beast123.m`, as well as several test datasets such as `Nile.mat` and `co2.mat`.

Precompiled MEX binaries are provided for selected platforms, including Windows 10, Ubuntu 22.04, macOS High Sierra 10.13 on Intel, and recent macOS versions on Apple Silicon. If the precompiled binary does not work on your machine, you can compile the MEX library from the C source files under the `Rbeast/Source` folder. If needed, we are happy to help compile Rbeast for your specific operating system or machine.

### Step 1: Download the C source files

- **Method 1:** Manually download all C source files from the [`Source`](https://github.com/zhaokg/Rbeast/tree/master/Source) folder.
- **Method 2:** After installing Rbeast in MATLAB, run the following command to automatically download the source files to your Rbeast folder:

```matlab
rbeast_src_download
```

### Step 2: Install a C/C++ compiler

The compiler depends on your operating system. The following are examples of compilers that have been used to build Rbeast.

#### Windows

For Windows, one option is the CRAN distribution of MinGW GCC through Rtools. More information is available from the [Rtools website](https://cran.r-project.org/bin/windows/Rtools/).

Typical GCC locations are:

- Rtools43: `C:\rtools43\mingw64\bin\gcc.exe`
- Rtools45: `C:\rtools45\x86_64-w64-mingw32.static.posix\bin\gcc.exe`

The installer may not configure the system path automatically. If needed, add the GCC directory to the `PATH` so that `gcc` can be accessed from `cmd.exe` or PowerShell.

In `cmd.exe`:

```cmd
set PATH=%PATH%;C:\rtools45\x86_64-w64-mingw32.static.posix\bin
```

In PowerShell:

```powershell
$env:Path = $env:Path + ";C:\rtools45\x86_64-w64-mingw32.static.posix\bin"
```

#### Linux/Ubuntu

GCC is often preinstalled on Linux. If needed, install the standard build tools:

```bash
sudo apt install build-essential
```

#### macOS

Install Xcode or the Xcode command-line tools to get `clang`. Alternatively, install LLVM/Clang with [Homebrew](https://stackoverflow.com/questions/8674546/how-to-update-llvm-and-clang-on-mac-os-x):

```bash
# Install Homebrew
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install LLVM/Clang
brew install llvm
```

---

## Step 3a: Build the MEX binary in MATLAB using `mex`

If the MATLAB `mex` command is properly configured, there are two ways to build the Rbeast MEX binary within MATLAB's console.

### Method 1: Use the Rbeast helper functions

Run the following commands in MATLAB:

```matlab
rbeast_src_download
rbeast_src_compile
```

### Method 2: Call `mex` manually

Assuming the MATLAB current folder is the source folder, use one of the following commands.

#### Linux/macOS

```matlab
mex -v CFLAGS='-pthread -DM_RELEASE -UUSE_MEX_CMD -fPIC -O2 -Wall -std=gnu99 -mfpmath=sse -msse2 -mstackrealign' -lmwservices -lut -lpthread *.c -output Rbeast
```

#### Windows

```matlab
mex -v CFLAGS='-pthread -DM_RELEASE -UUSE_MEX_CMD -fPIC -O2 -Wall -std=gnu99 -mfpmath=sse -msse2 -mstackrealign' -lmwservices -lut -lpthread -lkernel32 -lgdi32 -luser32 *.c -output Rbeast
```

---

## Step 3b: Build the MEX binary directly using GCC or Clang

If MATLAB's `mex` command does not work, you can build the MEX binary directly with GCC or Clang. A MATLAB MEX file is a dynamic library: a `.mexw64` file on Windows, a `.mexa64` file on Linux, and a `.mexmaci64` or `.mexmaca64` file on macOS.

The key requirement is to locate MATLAB's include and library folders. Typical paths are listed below. Replace these paths with the locations on your own machine.

| Platform | Include path | Library path |
|---|---|---|
| macOS Intel | `/Applications/MATLAB_R2019b.app/extern/include` | `/Applications/MATLAB_R2019b.app/bin/maci64` |
| macOS Apple Silicon | `/Applications/MATLAB_R2023b.app/extern/include` | `/Applications/MATLAB_R2023b.app/bin/maca64` |
| Linux/Ubuntu | `/MATLAB/extern/include` | `/MATLAB/bin/glnxa64` |
| Windows with Visual Studio | `C:/Program Files/MATLAB/R2019a/extern/include` | `C:/Program Files/MATLAB/R2019a/extern/lib/win64/microsoft` |
| Windows with MinGW GCC | `C:/Program Files/MATLAB/R2019a/extern/include` | `C:/Program Files/MATLAB/R2019a/extern/lib/win64/mingw64` |
| Windows with Octave | `C:/Program Files/GNU Octave/Octave-8.2.0/mingw64/include/octave-8.2.0/octave` | `C:/PROGRA~1/GNUOCT~1/OCTAVE~1.0/mingw64/lib/octave/8.2.0` |

The commands below assume that the current terminal folder is the source folder containing the Rbeast `.c` files.

---

### macOS Apple Silicon

```bash
clang -arch arm64 -mmacosx-version-min=12.0 -dynamiclib -fPIC \
 ./*.c \
 -I"/Applications/MATLAB_R2025a.app/extern/include" \
 -DM_RELEASE -DNDEBUG -D_REENTRANT -DMATLAB_MEX_FILE \
 -DMATLAB_DEFAULT_RELEASE=R2017b -DUSE_MEX_CMD \
 -fno-common -fwrapv -ffp-contract=off -fexceptions -O2 \
 -L"/Applications/MATLAB_R2025a.app/bin/maca64" \
 -weak-lmx -weak-lmex -weak-lmat -lpthread -lm -lmwservices -lut \
 -o Rbeast.mexmaca64
```

This is a single command split across multiple lines with backslashes (`\`). Masure sure that the backslash should be the final character on each continued line, with no trailing spaces after it. For other OS environments and consoles, the line continutaiton mark may be the blackslash ***\\***, caret ***^***, or backtick ***`***.

---

### macOS Intel

```bash
clang -mmacosx-version-min=10.13 -dynamiclib -fPIC \
 ./*.c \
 -I"/Applications/MATLAB_R2019b.app/extern/include" \
 -DM_RELEASE -DNDEBUG -D_REENTRANT -DMATLAB_MEX_FILE \
 -DMATLAB_DEFAULT_RELEASE=R2017b -DUSE_MEX_CMD \
 -Wall -msse2 -mstackrealign -fno-common -fexceptions -O2 \
 -L"/Applications/MATLAB_R2019b.app/bin/maci64" \
 -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress \
 -lpthread -lm -lut -lmwservices -lmat -lmex -lmx \
 -o Rbeast.mexmaci64
```

---

### Linux/Ubuntu

```bash
gcc -shared -fPIC -pthread -std=gnu99 \
 ./*.c \
 -I/MATLAB/extern/include \
 -DM_RELEASE -DNDEBUG -D_REENTRANT -DMATLAB_MEX_FILE \
 -DMATLAB_DEFAULT_RELEASE=R2017b -DUSE_MEX_CMD \
 -Wall -mfpmath=sse -msse2 -mstackrealign -fno-common -fexceptions -O2 \
 -L/MATLAB/bin/glnxa64 \
 -Wl,--export-all-symbols \
 -lm -lut -lmwservices -lmat -lmex -lmx \
 -o Rbeast.mexa64
```

---

### Windows for MATLAB with MinGW GCC

#### Earlier MinGW/GCC builds

Some earlier Windows MinGW/GCC environments expand wildcards such as `./*.c` automatically. If this works on your system, the following command can be used in `cmd.exe`:

```cmd
gcc -shared -fPIC -pthread -O2 ^
 ./*.c ^
 -I"C:/Program Files/MATLAB/R2019a/extern/include" ^
 -DM_RELEASE -DNDEBUG -D_REENTRANT -DMATLAB_MEX_FILE ^
 -DMATLAB_DEFAULT_RELEASE=R2017b -DUSE_MEX_CMD ^
 -Wall -mfpmath=sse -msse2 -mstackrealign -fno-common -fexceptions -O2 ^
 -L"C:/Program Files/MATLAB/R2019a/extern/lib/win64/mingw64" ^
 -L"C:/Program Files/MATLAB/R2019a/bin/win64" ^
 -Wl,--export-all-symbols ^
 -lmx -lmex -lmat -lut -lm ^
 -l:libmwservices.dll ^
 -lkernel32 -lgdi32 -luser32 ^
 -o Rbeast.mexw64
```

#### Recent MinGW/GCC builds in `cmd.exe`

On Windows, wildcard expansion depends on the shell and runtime environment. In some recent MinGW/GCC distributions, commands such as `gcc ./*.c ...` may pass the literal string `./*.c` to GCC rather than expanding it to the list of source files. If that happens, use a GCC response file:

```cmd
del "%TEMP%\cfilelist.txt" 2>nul
for %f in (*.c) do @echo %f>>"%TEMP%\cfilelist.txt"
type "%TEMP%\cfilelist.txt"

gcc -shared -fPIC -pthread -O2 ^
 @"%TEMP%\cfilelist.txt" ^
 -I"C:/Program Files/MATLAB/R2019a/extern/include" ^
 -DM_RELEASE -DNDEBUG -D_REENTRANT -DMATLAB_MEX_FILE ^
 -DMATLAB_DEFAULT_RELEASE=R2017b -DUSE_MEX_CMD ^
 -Wall -mfpmath=sse -msse2 -mstackrealign -fno-common -fexceptions -O2 ^
 -L"C:/Program Files/MATLAB/R2019a/extern/lib/win64/mingw64" ^
 -L"C:/Program Files/MATLAB/R2019a/bin/win64" ^
 -Wl,--export-all-symbols ^
 -lmx -lmex -lmat -lut -lm ^
 -l:libmwservices.dll ^
 -lkernel32 -lgdi32 -luser32 ^
 -o Rbeast.mexw64

del "%TEMP%\cfilelist.txt"
```

The first three lines create a temporary response file containing the list of `.c` files in the current folder. The response file is then passed to GCC as `@"%TEMP%\cfilelist.txt"`. After compilation, the temporary file is deleted.

In `cmd.exe`, the caret (`^`) must be the final character on a continued line. Do not put spaces after it.

#### Recent MinGW/GCC builds in PowerShell

```powershell
$src = Get-ChildItem .\*.c | ForEach-Object { $_.FullName }

gcc -shared -fPIC -pthread -O2 `
 $src `
 -I"C:/Program Files/MATLAB/R2019a/extern/include" `
 -DM_RELEASE -DNDEBUG -D_REENTRANT -DMATLAB_MEX_FILE `
 -DMATLAB_DEFAULT_RELEASE=R2017b -DUSE_MEX_CMD `
 -Wall -mfpmath=sse -msse2 -mstackrealign -fno-common -fexceptions -O2 `
 -L"C:/Program Files/MATLAB/R2019a/extern/lib/win64/mingw64" `
 -L"C:/Program Files/MATLAB/R2019a/bin/win64" `
 "-Wl,--export-all-symbols" `
 -lmx -lmex -lmat -lut -lm `
 -l:libmwservices.dll `
 -lkernel32 -lgdi32 -luser32 `
 -o Rbeast.mexw64
```

PowerShell uses the backtick (`` ` ``) for line continuation. The backtick must be the final character on the line, with no trailing spaces after it. The linker flag `"-Wl,--export-all-symbols"` is quoted because the comma has special meaning in PowerShell.

#### Note on `libmwservices.dll` and `ioFlush`

Some builds use the undocumented MATLAB internal function `ioFlush`, which is exported by `libmwservices.dll`. MATLAB does not provide a standard MinGW import library for this internal symbol, so the command above links directly against the DLL using:

```cmd
-l:libmwservices.dll
```

This requires the MATLAB runtime folder, such as `C:/Program Files/MATLAB/R2019a/bin/win64`, to be included in the library search path with `-L`. Because `ioFlush` is an undocumented internal MATLAB function, this interface may change across MATLAB releases.

---

### Windows for Octave with MinGW GCC

#### Earlier MinGW/GCC builds

```cmd
gcc -shared -fPIC -pthread -O2 ^
 ./*.c ^
 -I"C:/Program Files/GNU Octave/Octave-8.2.0/mingw64/include/octave-8.2.0/octave" ^
 -DM_RELEASE -DO_INTERFACE -DNDEBUG -D_REENTRANT -DMATLAB_MEX_FILE ^
 -DMATLAB_DEFAULT_RELEASE=R2017b -DUSE_MEX_CMD ^
 -Wall -mfpmath=sse -msse2 -mstackrealign -fno-common -fexceptions -O2 ^
 -L"C:/PROGRA~1/GNUOCT~1/OCTAVE~1.0/mingw64/lib/octave/8.2.0" ^
 -Wl,--export-all-symbols ^
 -loctinterp -loctave ^
 -lkernel32 -lgdi32 -luser32 ^
 -o Rbeast.mex
```

#### Recent MinGW/GCC builds in `cmd.exe`

```cmd
del "%TEMP%\cfilelist.txt" 2>nul
for %f in (*.c) do @echo %f>>"%TEMP%\cfilelist.txt"
type "%TEMP%\cfilelist.txt"

gcc -shared -fPIC -pthread -O2 ^
 @"%TEMP%\cfilelist.txt" ^
 -I"C:/Program Files/GNU Octave/Octave-8.2.0/mingw64/include/octave-8.2.0/octave" ^
 -DM_RELEASE -DO_INTERFACE -DNDEBUG -D_REENTRANT -DMATLAB_MEX_FILE ^
 -DMATLAB_DEFAULT_RELEASE=R2017b -DUSE_MEX_CMD ^
 -Wall -mfpmath=sse -msse2 -mstackrealign -fno-common -fexceptions -O2 ^
 -L"C:/PROGRA~1/GNUOCT~1/OCTAVE~1.0/mingw64/lib/octave/8.2.0" ^
 -Wl,--export-all-symbols ^
 -loctinterp -loctave ^
 -lkernel32 -lgdi32 -luser32 ^
 -o Rbeast.mex

del "%TEMP%\cfilelist.txt"
```

#### Recent MinGW/GCC builds in PowerShell

```powershell

gcc -shared -fPIC -pthread -O2 `
  (Get-ChildItem .\*.c | ForEach-Object Name) `
 -I"C:/Program Files/GNU Octave/Octave-8.2.0/mingw64/include/octave-8.2.0/octave" `
 -DM_RELEASE -DO_INTERFACE -DNDEBUG -D_REENTRANT -DMATLAB_MEX_FILE `
 -DMATLAB_DEFAULT_RELEASE=R2017b -DUSE_MEX_CMD `
 -Wall -mfpmath=sse -msse2 -mstackrealign -fno-common -fexceptions -O2 `
 -L"C:/PROGRA~1/GNUOCT~1/OCTAVE~1.0/mingw64/lib/octave/8.2.0" `
 "-Wl,--export-all-symbols" `
 -loctinterp -loctave `
 -lkernel32 -lgdi32 -luser32 `
 -o Rbeast.mex
```

---

## Final notes

- Replace all MATLAB, Octave, include, and library paths with the paths on your own machine.
- If the Rbeast C/C++ source folder is not the current working directory, replace `*.c` or `./*.c` with the path to your source files.
- In `cmd.exe`, use `%f` for interactive commands. In a batch file, use `%%f` instead.
- In PowerShell, quote `"-Wl,--export-all-symbols"` because the comma has special meaning.
- When linking with GCC/MinGW, place source files or object files (e.g., `./*.c`) before the libraries (e.g., `-lmx -lmex `) that resolve their symbols.
- In GCC, add the `-save-temps` to keep all the intermdiate files generated in the compilation
- In GCC, add `Wl,--out-implib,Rbeast.dll.a` if an import lib is neeed.
