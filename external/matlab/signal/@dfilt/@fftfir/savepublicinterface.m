function s = savepublicinterface(this)
%SAVEPUBLICINTERFACE   Save the public interface.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

s = abstract_savepublicinterface(this);

s.BlockLength         = get(this, 'BlockLength');
s.NonProcessedSamples = get(this, 'NonProcessedSamples');

% [EOF]
