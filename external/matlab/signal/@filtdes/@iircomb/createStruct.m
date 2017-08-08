function [s, str] = createStruct(h)
%CREATESTRUCT Return the response types

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

str = 'iircombFilterTypes';

s(1).construct = 'filtdes.iirnotch';
s(1).tag       = 'Notching';

s(2).construct = 'filtdes.iirpeak';
s(2).tag       = 'Peaking';

% [EOF]
