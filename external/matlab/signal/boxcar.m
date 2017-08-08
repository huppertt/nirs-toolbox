function w = boxcar(n_est)
%BOXCAR Boxcar window.
%   BOXCAR still works but maybe removed in the future. Use
%   RECTWIN instead.  Type help RECTWIN for details.
% 
%   See also BARTLETT, BLACKMAN, CHEBWIN, HAMMING, HANN, KAISER
%   and TRIANG.

%   Copyright 1988-2004 The MathWorks, Inc.

[n,w,trivialwin] = check_order(n_est);
if trivialwin, return, end;

w = ones(n,1);

% [EOF] boxcar.m
