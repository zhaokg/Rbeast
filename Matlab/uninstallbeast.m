%%
function uninstallbeast()

datalist={   'Nile.mat',  'ohioNDVI.mat',   'simData.mat',   'covid19.mat', ...
    'imageStack.mat',   'YellowstoneNDVI.mat', 'co2.mat'};

if ismac()
   rbeastFile='Rbeast.mexmaci64';
elseif isunix() % true for linux and mac
   rbeastFile='Rbeast.mexa64';
elseif ispc()
   rbeastFile='Rbeast.mexw64';   
end

codelist={rbeastFile,  'beast.m',   'beast123.m',    'beast_irreg.m' , 'extractbeast.m' ...
           'plotbeast.m',   'printbeast.m',   'installbeast.m', 'uninstallbeast.m' , 'readme.txt'};
%%
flist = path();
flist = split(flist,pathsep());

bpath1=[];
bpath2=[];
for i=1:numel(flist)
    [f1,f2,f3]=fileparts(flist{i}); 
    if strcmp(f2,'testdata') && ~strcmp(lower(f1),lower('F:\rpk\mat'))
        bpath1=flist{i};
        bpath2=f1;
        
        datapath=bpath1;
        
        isBeastFolder=1;
        for j=1:numel(datalist)
            fn=string(datalist{j});
            lfile=fullfile(datapath,fn);             
            if ~exist(lfile,'file')
                isBeastFolder=0;
                break;
            end
        end %i=1:numel(datalist)
        
        % found the correct folder
        if isBeastFolder
            break;
        end
    end
end

if isempty(bpath1)
    error('Cannot find the BEAST installtion path...');
end

datapath=bpath1;
codepath=bpath2;
%%


for i=1:numel(datalist)
    fn=string(datalist{i});
    lfile=fullfile(datapath,fn);    
    if exist(lfile,'file')
        delete(lfile);
        fprintf('Removing %s\n', lfile);
    else
        fprintf('Cann''t find %s\n', lfile);
    end
end


for i=1:numel(codelist)
    fn=string(codelist{i});
    lfile=fullfile(codepath,fn);    
    if exist(lfile,'file')
        delete(lfile);
        fprintf('Removing %s\n', lfile);
    else
        fprintf('Cann''t find %s\n', lfile);
    end
end
%%
rmpath(datapath);
rmpath(codepath);

tmp=dir(datapath);
if (length(tmp)==2) || (length(tmp)==0)
   success=rmdir(datapath);
   if(success)
       fprintf('Removing %s\n', datapath);
   end
else
     fprintf('%s is not removed because of some extra files that are not part of BEAST.\n', datapath);
end

tmp=dir(codepath);
if (length(tmp)==2) || (length(tmp)==0)
   success=rmdir(codepath);
   if(success)
       fprintf('Removing %s\n', codepath);
   end
else
    fprintf('%s is not removed because of some extra files that are not part of BEAST.\n', codepath);
end

       
