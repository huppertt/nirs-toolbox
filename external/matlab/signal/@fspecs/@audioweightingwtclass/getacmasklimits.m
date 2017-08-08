function [upperClass1 lowerClass1 upperClass2 lowerClass2] = getacmasklimits(~)
%GETACMASKUPPERLOWERLIMITS   Get the mask upper and lower limits for A and C
%weighting.

%   Copyright 2009 The MathWorks, Inc.

% Tolerances as defined in IEC 61672-1 2002-05 standard for A and C weighting

% Upper limit class 1
upperClass1 = [ 3.5 3 2.5 2.5 2.5 2 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.4 1.4 ...
    1.4 1.4 1.4 1.4 1.1 1.4 1.6 1.6 1.6 1.6 1.6 2.1 2.1 2.1 2.6 3 3.5 4 ].';

% Lower limit class 1
lowerClass1 = [ inf inf 4.5 2.5 2 2 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.4 1.4 ...
    1.4 1.4 1.4 1.4 1.1 1.4 1.6 1.6 1.6 1.6 1.6 2.1 2.6 3.1 3.6 6 17 inf ].';

% Upper limit class 2
upperClass2 = [ 5.5 5.5 5.5 3.5 3.5 3.5 2.5 2.5 2.5 2.5 2 2 2 2 1.9 1.9 1.9 ...
    1.9 1.9 1.9 1.4 1.9 2.6 2.6 3.1 3.1 3.6 4.1 5.1 5.6 5.6 6 6 6 ].';

% Lower limit class 2
lowerClass2 = [ inf inf inf 3.5 3.5 3.5 2.5 2.5 2.5 2.5 2 2 2 2 1.9 1.9 1.9 ...
    1.9 1.9 1.9 1.4 1.9 2.6 2.6 3.1 3.1 3.6 4.1 5.1 5.6 inf inf inf inf ].';

% [EOF]