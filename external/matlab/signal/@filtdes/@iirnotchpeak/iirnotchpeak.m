function h = iirnotchpeak
%IIRNOTCHPEAK Construct an IIRNOTCHPEAK object

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

h = filtdes.iirnotchpeak;

set(h, 'Tag', 'IIR Notch/Peak');

designMethodwFs_construct(h);

% [EOF]
