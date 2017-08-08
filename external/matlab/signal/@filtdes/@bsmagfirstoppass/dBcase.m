function dBcase(h,d)
%DBCASE Handle the dB case.
%
% This should be a private method.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

convertmag(h,d,...
    {'Dstop', 'Dpass2'},...
    {'Astop', 'Apass2'},...
    {'stop', 'pass'},...
    'todb');

% [EOF]
