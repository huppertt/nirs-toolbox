function h = iirpeak
%IIRPEAK Construct an IIRPEAK object

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

h = filtdes.iirpeak;

filterType_construct(h);

% [EOF]
