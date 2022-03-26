function x=extractbeast(o, index)
%  
%   USAGE: <strong>extractbeast(o, index) </strong>
%
%   <strong>o</strong>:  the time series analysis output from  beast123; o
%   should contain results for mulltiple time series
%
%   <strong>index </strong>: if o contains results for more than 1 time
%   series, index specifies for which time series the result is extracted.
%   If o is the result for a 3D stacked cube, index will be a vector of 2
%   integer to specify the row and col of the desired pixel. If o contains
%   only one time series, index will be ignored

%      
 if length(o.marg_lik)==1
     x=o;
 else
     x=Rbeast('tsextract',o,index);
 end
    
 