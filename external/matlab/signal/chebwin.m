function w = chebwin(n_est, r)
%CHEBWIN Chebyshev window.
%   CHEBWIN(N) returns an N-point Chebyshev window in a column vector.
% 
%   CHEBWIN(N,R) returns the N-point Chebyshev window with R decibels of
%   relative sidelobe attenuation. If omitted, R is set to 100 decibels.
%
%   % Example:
%   %   Create a 64-point Chebyshev window with 100 dB of sidelobe 
%   %   attenuation and display the result using WVTool.
%
%   L=64;
%   wvtool(chebwin(L))
%
%   See also TAYLORWIN, GAUSSWIN, KAISER, TUKEYWIN, WINDOW.

%   Author: James Montanaro
%   Reference: E. Brigham, "The Fast Fourier Transform and its Applications" 

%   Copyright 1988-2013 The MathWorks, Inc.

narginchk(1,2);

% Default value for R parameter.
if nargin < 2 || isempty(r), 
    r = 100.0;
end

% Cast to enforce precision rules. 
n_est = signal.internal.sigcasttofloat(n_est,'double','chebwin','N','allownumeric');
r = signal.internal.sigcasttofloat(r,'double','chebwin','R','allownumeric');

[n,w,trivialwin] = check_order(n_est);
if trivialwin, return, end;

if r < 0,
    error(message('signal:chebwin:MustBePositive'));
end

w = chebwinx(n,r);


% [EOF] chebwin.m
