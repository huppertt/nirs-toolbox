function nondBcase(h,d)
%nondBcase Handle the linear case.
%
% This should be a private method.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

convertmag(h,d,...
    {'Apass1', 'Apass2'},...
    {'Dpass1', 'Dpass2'},...
    {'pass', 'pass'},...
    @tolinear);

% [EOF]
