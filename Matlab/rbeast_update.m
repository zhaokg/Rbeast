% <strong>rbeast_update:</strong>  Check if a newer version is available from GitHub or not. If yes,
% the current local verson will be automatically deleted and the new
% version will be donwloaded and installed

beastPath=rbeast_path();

if ~isempty(beastPath)
    fprintf("Local version:\n");
    rbeast_version;
    localVersion =rbeastGitHubVersion;
    isRemoteVersionAvaiable=0;
    
    try
        fprintf("Github version:\n");
        eval(webread('https://github.com/zhaokg/Rbeast/raw/master/Matlab/rbeast_version.m',weboptions('CertificateFilename','')))
        isRemoteVersionAvaiable=1;
    end
    if (isRemoteVersionAvaiable==0)
        error("Can't access Github");
    end
    remoteVersion =rbeastGitHubVersion;

    if(remoteVersion>localVersion)
        strMsgBeast = sprintf('A new version v%f is available from Github. Press y to install it: \n',remoteVersion);
        answerBeast = input(strMsgBeast,'s');
        if strcmp(answerBeast,'y') 
            hasBeastSrcFld = exist(fullfile(beastPath,'source'),'dir');
            rbeast_uninstall();
            eval(webread('http://b.link/rbeast',weboptions('CertificateFilename','')));
            if (hasBeastSrcFld)
                fprintf('\nContinue to re-downoad the C/C++ source filers \n' );
                pause(3);
                rbeast_src_download();
            end
        end
    else
        fprintf('The latest version is being used and no update is available from GitHub. \n' );
    end
end


if isempty(beastPath)
      eval(webread('https://github.com/zhaokg/Rbeast/raw/master/Matlab/installbeast.m',weboptions('CertificateFilename','')));    
end
%%
clearvars  rbeastGitHubVersion localVersion remoteVersion isRemoteVersionAvaiable strMsgBeast hasBeastSrcFld answerBeast