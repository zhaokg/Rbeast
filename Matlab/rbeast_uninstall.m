function rbeast_uninstall()
%%
codeitself = webread('https://b.link/delbeast',weboptions('cert',''));
%%
if ismac()
   rbeastFile='Rbeast.mexmaci64';
elseif isunix() % true for linux and mac
   rbeastFile='Rbeast.mexa64';
elseif ispc()
   rbeastFile='Rbeast.mexw64';   
end

datalist={   'Nile.mat',  'ohioNDVI.mat',   'simData.mat',   'covid19.mat', 'googletrend.mat', ...
             'imageStack.mat',   'YellowstoneNDVI.mat', 'co2.mat'};
         
codelist={'Rbeast.mexw64','Rbeast.mexmaci64', 'Rbeast.mexa64', 'beast.m',   'beast123.m',    'beast_irreg.m' , 'extractbeast.m' ...
           'plotbeast.m',   'printbeast.m',   'rbeast_uninstall.m' , 'rbeast_update.m', ...
           'rbeast_version.m','rbeast_path.m', 'readme.txt'};


%'installrbeast.m' is only available online and not downloaded to local paths       
oldcodelist={'installbeast.m', 'uninstallbeast.m'};
%%
beastpath = rbeast_path();
if isempty(beastpath)
    error('Cannot find the BEAST installtion path...');
end
codepath  = beastpath;
datapath  = fullfile(beastpath,'testdata');
 
for i=1:numel(datalist) 
    lfile=fullfile(datapath,  string(datalist{i}) );    
    if exist(lfile,'file')
        delete(lfile);
        fprintf('Removing %s\n', lfile);
    else
        fprintf('Cann''t find %s\n', lfile);
    end
end


for i=1:numel(codelist)
    lfile=fullfile(codepath, string(codelist{i}));    
    if exist(lfile,'file')
        delete(lfile);
        fprintf('Removing %s\n', lfile);
    else
        fprintf('Cann''t find %s\n', lfile);
    end
end

for i=1:numel(oldcodelist)
    lfile=fullfile(codepath, string(oldcodelist{i}));    
    if exist(lfile,'file')
        delete(lfile);
        fprintf('Removing %s\n', lfile);
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

       
