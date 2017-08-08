function cmd = maskinfo(h, d)
%MASKINFO

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

bw = getbandwidth(d);
fnotch = get(d, 'Fnotch');

cmd{1}.frequency  = [0 fnotch-bw/2];
cmd{1}.filtertype = 'lowpass';
cmd{1}.freqfcn    = 'wpass';

cmd{2}.frequency  = [fnotch+bw/2 getnyquist(d)];
cmd{2}.filtertype = 'highpass';
cmd{2}.freqfcn    = 'wpass';

% [EOF]
