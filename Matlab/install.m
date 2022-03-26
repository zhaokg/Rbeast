%beastPath='a:\xx\'
eval( webread('https://go.osu.edu/rbeast',weboptions('Cert','')) )

if ~exist('beastPath','var')
    warning("The variable 'beastPath' doesnot exist; a temporaby folder is used instead");
    beastPath =  tempdir()+"Rbeast\";
    disp("BEAST installation Path: " +beastPath)
end

if ~exist(beastPath,"dir")
    success=mkdir(beastPath)
    if ~success
        error("Cannot create or access the beast path specified.");
    end
else
    tmpfile=fullfile(beastPath, num2str(datenum(date())+rand(1),'%10.6f') );
    fid=fopen(tmpfile,'w+');
    if fid==-1
        error("Cannot wiret or access the beast path specified.");
    else
        fclose(fid);
        delete(tmpfile);
    end
end

datapath=fullfile(beastPath,'testdata');
if ~exist(datapath,"dir")
    mkdir(datapath);
end
%%

rpath ="https://github.com/zhaokg/Rbeast/raw/master/Matlab/";
fn ="Nile.mat";
lfile=beastPath+fn;
rfile=rpath+fn;

datalist={   'Nile.mat',  'ohioNDVI.mat',   'simData.mat',   'covid19.mat', ...
    'imageStack.mat',   'YellowstoneNDVI.mat'};

for i=1:numel(datalist)
    fn=string(datalist{i});
    lfile=fullfile(datapath,fn);
    rfile=rpath+"testdata/"+fn;
    websave(lfile, rfile,weboptions('Cert',[])); 
    fprintf('Downloaded: %s\n', lfile);
end

codelist={   'beast.m',   'beast123.m',    'beast_irreg.m' , 'extractbeast.m' ...
    'plotbeast.m',   'printbeast.m',   'installrbeast.m'};


for i=1:numel(datalist)
    fn=string(codelist{i});
    lfile=fullfile(beastPath,fn);
    rfile=rpath+fn;
    websave(lfile, rfile,weboptions('Cert',[]));    
    fprintf('Downloaded: %s\n', lfile);
end
%%
addpath(beastPath);
addpath(datapath);
addpath(genpath(beastPath) );
savepath
%%
clc
fprintf('... Rbeast was installed at %s\n', beastPath);
fprintf("... '%s' and '%s' are added to the search path. \n     Make sure to add these two paths back (e.g., addpath()) after re-starting Matlab\n", beastPath, datapath);
fprintf("... Run <strong>'help beast'</strong>, <strong>'help beast123'</strong>, or <strong>'help beast_irreg'</strong> for usage and examples\n");
