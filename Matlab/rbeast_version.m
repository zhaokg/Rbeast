
% Oct 2., 2022: fix a small bug "int ParseInputData( BEAST2_IO_PTR _OUT_ io)"
% The pixel indices are 1-based not 0-based. This bug has no effects on the 
% model results but just improve the robutness of the program

% The larger the number, the higher/latest version it is
rbeastGitHubVersion = 0.9412;