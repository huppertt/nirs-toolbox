function w = kaiser(n_est,bta)
%KAISER Kaiser window.
%   W = KAISER(N) returns an N-point Kaiser window in the column vector W.
% 
%   W = KAISER(N,BTA) returns the BETA-valued N-point Kaiser window.
%       If omitted, BTA is set to 0.500.
%
%   % Example:
%   %   Create a 200-point Kaiser window with a beta of 2.5 and display 
%   %   the result using WVTool.
%
%   w = kaiser(200,2.5);
%   wvtool(w)
%
%   See also CHEBWIN, GAUSSWIN, TUKEYWIN, WINDOW.

%   Copyright 1988-2013 The MathWorks, Inc.

narginchk(1,2);

% Default value for the BETA parameter.
if nargin < 2 || isempty(bta), 
    bta = 0.500;
end
% Cast to enforce Precision Rules
n_est = signal.internal.sigcasttofloat(n_est,'double','kaiser','N',...
  'allownumeric');
bta = signal.internal.sigcasttofloat(bta,'double','kaiser','BTA',...
  'allownumeric');

[nn,w,trivialwin] = check_order(n_est);
if trivialwin, return, end;

nw = round(nn);
bes = abs(besseli(0,bta));
odd = rem(nw,2);
xind = (nw-1)^2;
n = fix((nw+1)/2);
xi = (0:n-1) + .5*(1-odd);
xi = 4*xi.^2;
w = besseli(0,bta*sqrt(1-xi/xind))/bes;
w = abs([w(n:-1:odd+1) w])';

    
% [EOF] kaiser.m
