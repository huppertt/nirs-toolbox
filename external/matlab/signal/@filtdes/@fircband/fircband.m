function h = fircband
%FIRCBAND Construct this object

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

h = filtdes.fircband;

abstractgremez_construct(h);

set(h, 'Tag', 'Constrained Band Equiripple FIR');

% [EOF]
