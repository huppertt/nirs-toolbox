function w = bartlett(n_est)
%BARTLETT Bartlett window.
%   W = BARTLETT(N) returns the N-point Bartlett window.
%
%   % Example:
%   %   Create a 64-point Bartlett window and display the result using 
%   %   WVTool.
%
%   L=64;
%   wvtool(bartlett(L))
%
%   See also BARTHANNWIN, BLACKMANHARRIS, BOHMANWIN, FLATTOPWIN,
%            NUTTALLWIN, PARZENWIN, RECTWIN, TRIANG, WINDOW.

%   Copyright 1988-2002 The MathWorks, Inc.

narginchk(1,1);
% Cast to enforce Precision Rules
n_est = signal.internal.sigcasttofloat(n_est,'double','bartlett','N',...
  'allownumeric');

[n,w,trivialwin] = check_order(n_est);
if trivialwin, return, end;

w = 2*(0:(n-1)/2)/(n-1);
if rem(n,2)
    % It's an odd length sequence
    w = [w w((n-1)/2:-1:1)]';
else
    % It's even
    w = [w w(n/2:-1:1)]';
end

% [EOF] bartlett.m
