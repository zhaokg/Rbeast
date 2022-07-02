function beastpath=rbeast_path()
% a silly way to enumerate the current search path to find the beast
% install path
datalist={'Nile.mat',  'ohioNDVI.mat',   'simData.mat',   'covid19.mat', 'imageStack.mat',   'YellowstoneNDVI.mat', 'co2.mat'};
%%
beastpath=[];
%%
flist = split(path(),pathsep());
 
for i=1:numel(flist)
    
    [fpath,fname,fext] = fileparts(flist{i}); 
    if strcmp(fname,'testdata') && ~strcmp(lower(fpath),lower('F:\rpk\mat'))

        datapath=flist{i};   %...\testdata
        
        isBeastFolder=1;
        for j=1:numel(datalist)         
            lfile=fullfile( datapath, string(datalist{j}) );             
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
    error('Cannot find the BEAST installtion path...');
end
%datapath=bpath1;
%codepath=bpath2;