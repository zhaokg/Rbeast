function rbeast_uninstall()
% 
% <strong>rbeast_uninstall:</strong> Remove the beast's program from the
% local machine

%% Useless code
codeitself = webread('https://b.link/delbeast',weboptions('CertificateFilename',''));

if ismac()
   rbeastFile='Rbeast.mexmaci64';
elseif isunix() % true for linux and mac
   rbeastFile='Rbeast.mexa64';
elseif ispc()
   rbeastFile='Rbeast.mexw64';   
end

%%
beastpath = rbeast_path();
if isempty(beastpath)
    error('Cannot find the BEAST installtion path...');
end

codepath  = beastpath;
datapath  = fullfile(beastpath,'testdata');
srcpath   = fullfile(beastpath,'source' ); 
%%
datalist={   'Nile.mat',  'ohioNDVI.mat',   'simData.mat',   'covid19.mat', 'googletrend.mat', ...
             'imageStack.mat',   'YellowstoneNDVI.mat', 'co2.mat'};
         
codelist={ 'Rbeast.mex', 'Rbeast.mexw64','Rbeast.mexmaci64', 'Rbeast.mexa64', ...
           'beast.m',   'beast123.m',    'beast_irreg.m' , 'extractbeast.m', 'plotbeast.m',   'printbeast.m', ...
		   'rbeast_install.m', 'rbeast_uninstall.m' , 'rbeast_update.m', 'rbeast_version.m','rbeast_path.m', ...
		   'rbeast_src_compile.m','rbeast_src_download.m','readme.txt','readme.md'};

% readme.txt has be renamed to readme.md. We keep readme.txt here for the old versions.		    

remove_manyfiles(datapath, datalist); % delete(file)
remove_manyfiles(codepath, codelist); % delete(file)

%'installrbeast.m' is only available online and not downloaded to local paths       
%oldcodelist={'installbeast.m', 'uninstallbeast.m'}; 
%remove_manyfiles(codepath, codelist);
 
%% If there is a src folder, also delete everything
allsrcFile = fullfile(srcpath, '*.*' );
if exist(srcpath,"dir")   
   delete(allsrcFile);  % delete('PathToFolder\*')
   rmdir(srcpath);
end

%%
rmpath(datapath);     % remove the path from the Matlab's search paths
rmpath(codepath);     % remove the path from the Matlab's search paths
 
remove_emptyfolder(datapath);
remove_emptyfolder(codepath);

end


%% Helper functions

%  remove a file and print the status
function remove_file(localFile)        
    if exist(localFile,'file')
        delete(localFile);
        fprintf('Removing %s\n', localFile);
    else
        fprintf('Cann''t find %s\n', localFile);
    end
end
 
% remove many files
function remove_manyfiles(folder, flist)        
 for i=1:numel(flist)
    fn = fullfile(folder, flist{i});
    remove_file(fn);
 end
end
 
% remove a file and print the status
function remove_file_quiet(localFile)        
    if exist(localFile,'file')
        delete(localFile);
        fprintf('Removing %s\n', localFile);    
    end
end
 
% remove a file and print the status
function remove_manyfiles_quite(folder, flist)        
 for i=1:numel(flist)
    fn = fullfile(folder, flist{i});
    remove_file_quiet(fn);
 end
end

% delete a folder
function remove_emptyfolder(fld)  
    flist = dir(fld);
    if (length(flist)==2) || (length(flist)==0)
       success=rmdir(fld);
       if(success)
           fprintf('Removing %s\n', fld);
       end
    else
         fprintf('%s is not removed because of some extra files or subfolders in the folder.\n', fld);
    end
end