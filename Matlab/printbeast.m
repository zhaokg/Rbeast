function x=printbeast(o, index)
%   USAGE: <strong>printbeast(o, index) </strong>
%
%   <strong>o</strong>:  the time series analysis output from  beast123; o
%   should contain results for mulltiple time series
%
%   <strong>index </strong>: if o contains results for more than 1 time
%   series, index specifies for which time series the result is printed.
%   If o is the result for a 3D stacked cube, index will be a vector of 2
%   integer to specify the row and col of the desired pixel. If o contains
%   only one time series, index will be ignored
%
%   <strong>Contact info</strong>: To report bug or get help, do not hesitate to contact Kaiguang Zhao
%   at <strong>zhao.1423@osu.edu</strong>.
%
 if strcmp(o.class,'beast')
     if nargin==1
         index=1;
     end
     Rbeast('print',o,index);
 else
     error('the input has to be an output from the BEAST functions');
 end
    
 