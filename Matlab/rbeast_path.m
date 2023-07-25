function beastpath=rbeast_path()
% <strong>rbeast_path:</strong> A silly way to enumerate the current search path to find the beast
% install path, if any
%%
beastpath=[];
%%
datalist={'Nile.mat',  'ohioNDVI.mat',   'simData.mat',   'covid19.mat', 'imageStack.mat',   'YellowstoneNDVI.mat', 'co2.mat'};

%flist = split(path(),pathsep());
flist = strsplit(path(),pathsep());
 
for i=1:numel(flist)
    
    [fpath,fname,fext] = fileparts(flist{i}); 
	fpath_lower        = lower(fpath);
    if strcmp(fname,'testdata') && ~strcmp(fpath_lower,'f:\rpk\mat') && isempty( strfind(fpath_lower,':\rpk\mat') ) ...
	   && isempty( strfind(fpath_lower,'rpk\mat') ) && isempty( strfind(fpath_lower,'rpk/mat') )

        datapath=flist{i};   %...\testdata
       
        isBeastFolder=1;
        for j=1:numel(datalist)         
            lfile=fullfile( datapath, datalist{j} );             
            if ~exist(lfile,'file')
                isBeastFolder=0;
                break;
            end
        end %i=1:numel(datalist)
        
        % Found the correct folder
        if isBeastFolder
            beastpath = fpath;
            break;
        end
    end
end

if isempty(beastpath)
    %warning('Cannot find the BEAST installtion path...');
end
%datapath=bpath1;
%codepath=bpath2;