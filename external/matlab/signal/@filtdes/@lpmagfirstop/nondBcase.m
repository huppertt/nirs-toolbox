function nondBcase(h,d)
%nondBcase Handle the linear case.
%
% This should be a private method.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

convertmag(h,d,...
    {'Astop'},...
    {'Dstop'},...
    {'stop'},...
    @tolinear);

% [EOF]
