% Aug 12, 2024: Remove the extra blank line when pringing. Also added the hasOutlier argument
% Dec 12, 2023: Added the the sseg.leftmargin, sseg.rightmargin, tseg.leftmargin, and
%   tseg.rightmargin argumentst upon the requst of Racine Nassau from University of Washington
% July 12, 2023: Added the rbeast_mex_compile.m: not always working but
%    may be usefulf or some users
% Feb 20, 2024: Add the output of CI for trend slope, at the requst of Dr. Igor Rizaev
% July 10, 2023: Fix a silly typo (toder -> torder):  Thanks to michalk vasnicka.
% July 8, 2023: Add the support for Octave:  Thanks to Mario Rossi
% Dec 8, 2002: Fixed the num_samples_per_period error for beast_irreg
%   Handle missing values, if any, in the time: Thanks to Gary Iacobucci <garyiaco@buffalo.edu> 
% Nov 24., 2022: The API has changed a little bit: the most important
%   change is that 'freq' is deprecated in favor of 'period'
% Oct 24., 2022: Thanks to Robert Martin: Suppressed unnecessary warnings 
%   related to whichDimIsTime=1
% Oct 13., 2022: Add the googletrend.mat file in the dataset list
% Oct 2., 2022: fix a small bug "int ParseInputData( BEAST2_IO_PTR _OUT_ io)"
%   The pixel indices are 1-based not 0-based. This bug has no effects on the 
%   model results but just improve the robutness of the program

% The larger the number, the higher/latest version it is
rbeastGitHubVersion = 0.9481;
disp(strcat('  matlab-beast v', num2str(rbeastGitHubVersion)));