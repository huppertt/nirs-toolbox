function loadpublicinterface(this, s)
%LOADPUBLICINTERFACE   Load the public interface.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

abstract_loadpublicinterface(this, s);

set(this, 'BlockLength', s.BlockLength, ...
    'NonProcessedSamples', s.NonProcessedSamples);

% [EOF]
