beastPath=rbeast_path();
if ~isempty(beastPath)
    rbeast_version;
    localVersion =rbeastGitHubVersion;

    isDone=0;
    try
        eval(webread('https://github.com/zhaokg/Rbeast/raw/master/Matlab/rbeast_version.m',weboptions('cert','')))
        isDone=1;
    end

    if (isDone==0)
        error("Can't access Github");
    end

    remoteVersion =rbeastGitHubVersion;

    if(remoteVersion>localVersion)
        str=sprintf('A new version v%f is available from Github. Press y to install it: \n',remoteVersion);
        y = input(str,'s');
        if strcmp(y,'y')        
            rbeast_uninstall();
            eval(webread('http://b.link/rbeast',weboptions('cert','')));
        end
    else
        fprintf('The latest version is being used and no update is available from GitHub. \n' );
    end
end


if isempty(beastPath)
      eval(webread('http://b.link/rbeast',weboptions('cert','')));    
end