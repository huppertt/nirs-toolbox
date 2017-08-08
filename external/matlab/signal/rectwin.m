function w = rectwin(n_est)
%RECTWIN Rectangular window.
%   W = RECTWIN(N) returns the N-point rectangular window.
%
%   % Example:
%   %   Create a 64-point rectangular window and display the result in 
%   %   WVTool:
%
%   L=64;
%   wvtool(rectwin(L))
%  
%   See also BARTHANNWIN, BARTLETT, BLACKMANHARRIS, BOHMANWIN, 
%            FLATTOPWIN, NUTTALLWIN, PARZENWIN, TRIANG, WINDOW.

%   Copyright 1988-2002 The MathWorks, Inc.

error(nargchk(1,1,nargin,'struct'));
[n,w,trivialwin] = check_order(n_est);
if trivialwin, return, end;

w = ones(n,1);


% [EOF] 

