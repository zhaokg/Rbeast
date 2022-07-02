   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %       Installation:
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Method 1: Manually copy all the files in this folder (e.g, beast.m and testdata) to your local folder
             Then, add your local folder to Matlab's search path by running:

             % Replcae c:\beast with your own drive
             >  addpath(genpath('c:\beast')) 

   Method 2: (Recommendded)
	    The easiest way to install is to simply run the following line of command in Matlab's console:
 
             >  eval(webread('http://b.link/beast',weboptions('cert','')))  % install to a temp folder
        
             Or install it to a chosen folder by first setting the beastPath variable

             >  beastPath='c:\beast'
             >  eval(webread('http://b.link/beast',weboptions('cert',''))) % install to the folder specified by beastPath


   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %       Uninstall
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   METHOD 1: Simply delete all the downloaded files manually
   METHOD 2: Automatically remove the files by running:
           
             > rbeast_uninstall

 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %       Usage and example
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

*** Major functions available:
    beast:            handle a single regular time series
    beast_irreg:      handle a single irregular time series
    beast123:         handle one or more time seires or stacked images 
    rbeast_uninstall: remove the installed files from the harddisk
    rbeast_update:    check github and update to the latest BEAST version\n");
fprintf("\n");
*** Examples
    load('Nile.mat')             % Nile river annual streamflow: trend-only data
    o=beast(Nile, 'start', 1871, 'season','none') 
    printbeast(o))
    plotbeast(o))

*** Run 'help beast', 'help beast123', or 'help beast_irreg' for usage and examples            