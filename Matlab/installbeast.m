%
% <strong>rbeast_install/installbeast:</strong> Download and install the BEAST matlab scripts
%                                      
%  'rbeast_install.m' and 'installbeast.m' are exaclty the same. installbeast.m is an online version of
%   rbeast_install.m.
%
%  'http://b.link/beast'  points to an online version of installbeast.m.
%
% Download the Matlab scripts and example datasets of BEAST from Github to a local folder
% specified by the variable "beastPath". If beastPath doesn't exist, a temporary folder will be 
% used. This is essentially the same as running:
%
%  <strong>eval( webread('http://b.link/rbeast',  weboptions('Cert','')) ) </strong>
%

%% 
% system("gcc -c -fPIC -pthread -DNDEBUG -DM_RELEASE -DMATLAB_DEFAULT_RELEASE=R2017b -I/MATLAB/extern/include -O2 -Wall  -std=gnu99 -mfpmath=sse -msse2 -mstackrealign  *.c")
% system("g++ -c -fPIC -pthread -DNDEBUG -DM_RELEASE -DMATLAB_DEFAULT_RELEASE=R2017b -I/MATLAB/extern/include  -O2 -Wall    -mfpmath=sse -msse2 -mstackrealign  *.cpp")
% system("gcc -shared -pthread -L/MATLAB/bin/glnxa64 -lmx -lmex -lmat -lm -lut -lmwservices *.o -o Rbeast.mexa64")

% mex -v CFLAGS='$CFLAGS  -DM_RELEASE -Wall -Wl,-v' -lmwservices -lut *.c -output Rbeast
% mex -v CFLAGS='-DM_RELEASE -Wall -Wl,-v' -lmwservices -lut *.c -output Rbeast

% mex -v CFLAGS='-DM_RELEASE -UUSE_MEX_CMD -fpic' -lmwservices -lut *.c -output Rbeast
% mex -v CFLAGS='-DM_RELEASE -UUSE_MEX_CMD -fPIC -O2 -Wall -std=gnu99 -march=native' -lmwservices -lut *.c -output Rbeast

%% Compilation for MEX

%%% Linux (Rbeast.mexa64):  Compiled using gcc on Ubunttu (Matlab online)
%
% mex -v CFLAGS='-pthread -DM_RELEASE -UUSE_MEX_CMD -fPIC -O2 -Wall -std=gnu99 -mfpmath=sse -msse2 -mstackrealign' -lmwservices -lut  -lpthread *.c -output Rbeast


%%% MacOS (Rbeast.mexmaci64): Compiled using clang on MacOS10.13 (VM)
%
%clang  -mmacosx-version-min=10.13 -dynamiclib  -fPIC -I"/Library/Frameworks/R.framework/Resources/include" -I/usr/local/include -I/Applications/MATLAB_R2020a.app/extern/include/ -DM_RELEASE -Wall -O2  -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress -L/Library/Frameworks/R.framework/Resources/lib -L/usr/local/lib -L/Applications/MATLAB_R2020a.app/bin/maci64/ -msse2 -mstackrealign -lpthread -lm -lut -lmwservices -lmat -lmex -lmx *.c -o Rbeast.mexmaci64
%clang  -mmacosx-version-min=10.13 -dynamiclib  -fPIC -I"/Library/Frameworks/R.framework/Resources/include" -I/usr/local/include -I/Applications/MATLAB_R2019b.app/extern/include/ -DM_RELEASE -Wall -O2  -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress -L/Library/Frameworks/R.framework/Resources/lib -L/usr/local/lib -L/Applications/MATLAB_R2019b.app/bin/maci64/ -msse2 -mstackrealign -lpthread -lm -lut -lmwservices -lmat -lmex -lmx *.c -o Rbeast.mexmaci64


%%%  Octave (Rbeast.mex): Compiled using gcc on Windows 10
%
% mex -v : get the compiler option and include and lib paths
%
% the src files must preceded the libraries (otherwise, mexPrintf etc. is
% not found) https://stackoverflow.com/questions/45135/why-does-the-order-in-which-libraries-are-linked-sometimes-cause-errors-in-gcc
%
% gcc -shared -fPIC -pthread -O2 -DM_RELEASE -DO_INTERFACE -I"C:\Program Files\GNU Octave\Octave-8.2.0\mingw64\include\octave-8.2.0\octave"    f:/rpk/src/*.c -LC:\PROGRA~1\GNUOCT~1\OCTAVE~1.0\mingw64\lib\octave\8.2.0  -Wl,--export-all-symbols   -loctinterp -loctave -lkernel32 -lgdi32 -luser32  -o Rbeast.mex
%
% I am using Rtools' mingw-gcc: c:\rtools43\mingw64\bin\
% gcc -shared -fPIC -pthread -O2 -DM_RELEASE -DO_INTERFACE -I"C:\Program Files\GNU Octave\Octave-8.2.0\mingw64\include\octave-8.2.0\octave"    f:/rpk/src/*.c -LC:\PROGRA~1\GNUOCT~1\OCTAVE~1.0\mingw64\lib\octave\8.2.0  -Wl,--export-all-symbols   -loctinterp -loctave -lkernel32 -lgdi32 -luser32  -o f:\rpk\mat\Rbeast.mex
%
% beastPath='a:\xx\'
% eval( webread('https://go.osu.edu/rbeast',weboptions('Cert','')) )
% eval( webread('http://bit.ly/loadbeast',  weboptions('Cert','')) )
% eval( webread('https://github.com/zhaokg/Rbeast/raw/master/Matlab/installbeast.m', weboptions('Cert','')) )

%%
clear Rbeast % just in case that an exisitng version has been loaded 

if ~exist('beastPath','var') 
    warning("The variable 'beastPath' does not exist; a temporaby folder is used instead");
    beastPath =  fullfile(tempdir(), "Rbeast");
    disp( strcat("BEAST installation Path: " , beastPath) );
end

if isempty(beastPath)     
    beastPath =  fullfile(tempdir(), "Rbeast");
    %disp("BEAST installation Path: " + beastPath);
	disp( strcat ("BEAST installation Path: " , beastPath) );
end


if ~exist(beastPath,"dir")
    beast_success=mkdir(beastPath);
    if ~ beast_success
        error("Cannot create or access the beast path specified.");
    end
else
    beast_tmpfile = fullfile(beastPath, num2str(datenum(date())+rand(1),'%10.6f') );
    beast_fid     = fopen(beast_tmpfile,'w+');
    if beast_fid==-1
        error("Cannot wiret or access the beast path specified.");
    else
        fclose(beast_fid);
        delete(beast_tmpfile);
    end
end

datapath = fullfile(beastPath,'testdata');
if ~exist(datapath,"dir")
    mkdir(datapath);
end

%%
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
rpath    = "https://github.com/zhaokg/Rbeast/raw/master/Matlab/";
%%
datalist={  'Nile.mat',  'ohioNDVI.mat',   'simData.mat',   'covid19.mat', 'googletrend.mat', ...
            'imageStack.mat',   'YellowstoneNDVI.mat', 'co2.mat'};

for i=1:numel(datalist)  
	fn    = datalist{i};                      %fn    = string(datalist{i});
    lfile = fullfile(datapath,fn);   
	rfile = strcat(rpath, "testdata/", fn);   %rfile = rpath+"testdata/"+fn;
	if isOctave
		urlwrite(rfile,lfile);
	else
	    websave(lfile, rfile,weboptions('CertificateFilename',[]));
	end
    fprintf('Downloaded: %s\n', lfile);
end

%% on the safe side, get all the mx library for all file systems
codelist={ 'Rbeast.mex', 'Rbeast.mexw64','Rbeast.mexmaci64', 'Rbeast.mexa64', ...
           'beast.m',   'beast123.m',    'beast_irreg.m' , 'extractbeast.m', 'plotbeast.m',   'printbeast.m', ...
		   'rbeast_install.m', 'rbeast_uninstall.m' , 'rbeast_update.m', 'rbeast_version.m','rbeast_path.m', ...
		   'rbeast_src_compile.m','rbeast_src_download.m','readme.md'};
       
% 'installrbeast.m' is only available online and not downloaded to local paths     

for i=1:numel(codelist)
    fn    = codelist{i}; %fn    = string(codelist{i});
    lfile = fullfile(beastPath,fn);
    rfile = strcat(rpath,fn); %rfile = rpath+fn;
	if isOctave
		urlwrite(rfile,lfile);
	else
	    websave(lfile, rfile,weboptions('CertificateFilename',[]));
	end
    fprintf('Downloaded: %s\n', lfile);
end

%%
addpath(beastPath);
addpath(datapath);
% addpath(genpath(beastPath) );
savepath
%%

%clc
fprintf('\n\n');
fprintf('*** <strong>Rbeast</strong> was installed at %s\n', beastPath);
fprintf("*** '%s' and '%s' are added to the search path. \n", beastPath, datapath);
% Not needed bcz savepath will save the searach path permanantly
% fprintf("     Make sure to add these two paths back (e.g., addpath()) after re-starting Matlab\n");
fprintf('\n');
fprintf("*** <strong>Major functions available</strong>:\n");
fprintf("    <strong>beast</strong>: handle a single regular time series\n");
fprintf("    <strong>beast_irreg</strong>: handle a single irregular time series\n");
fprintf("    <strong>beast123</strong>: handle one or more time seires or stacked images \n");
fprintf("    <strong>rbeast_uninstall</strong>: remove the installed files from the harddisk\n");
fprintf("    <strong>rbeast_update</strong>: check github and update to the latest BEAST version, if any\n");
fprintf("\n");
fprintf("*** <strong>Examples</strong>\n");
fprintf("    load('Nile.mat')             %% Nile river annual streamflow: trend-only data\n");
fprintf("    o=beast(Nile, 'start', 1871, 'season','none') \n");
fprintf("    printbeast(o)\n");
fprintf("    plotbeast(o)\n\n");

fprintf("*** Run <strong>'help beast'</strong>, <strong>'help beast123'</strong>, or <strong>'help beast_irreg'</strong> for usage and examples\n");


%% Not used here, but kept as a quick reminder
% https://stackoverflow.com/questions/24923384/how-to-get-matlab-to-determine-if-the-os-is-windows-or-mac-so-to-find-all-seri
if ismac()
   rbeastFile='Rbeast.mexmaci64';
elseif isunix() % true for linux and mac
   rbeastFile='Rbeast.mexa64';
elseif ispc()
   rbeastFile='Rbeast.mexw64';   
end

clearvars datapath codepath fn lfile rfile datalist codelist beast_success beast_fid beast_tmpfile rbeastFile isOctave 