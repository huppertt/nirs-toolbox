function cmd = maskinfo(h, d)
%MASKINFO

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

bw = getbandwidth(d);
fpeak = get(d, 'Fpeak');

cmd{1}.frequency  = fpeak + [-bw/2 bw/2];
cmd{1}.filtertype = 'bandpass';
cmd{1}.freqfcn    = 'wpass';

% [EOF]
