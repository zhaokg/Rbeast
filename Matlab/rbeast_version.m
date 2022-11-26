% Nov 24., 2022: The API has changed a little bit: the most important
% change is that 'freq' is deprecated in favor of 'period'

% Oct 24., 2022: Thanks to Robert Martin: Suppressed unnecessary warnings 
% related to whichDimIsTime=1

% Oct 13., 2022: Add the googletrend.mat file in the dataset list


% Oct 2., 2022: fix a small bug "int ParseInputData( BEAST2_IO_PTR _OUT_ io)"
% The pixel indices are 1-based not 0-based. This bug has no effects on the 
% model results but just improve the robutness of the program

% The larger the number, the higher/latest version it is
rbeastGitHubVersion = 0.945;