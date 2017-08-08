function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

fc = get(d, 'Fc');

if strcmpi(get(d, 'TransitionMode'), 'Rolloff'),
    bw = 2*fc*get(d, 'Rolloff');
else
    bw = get(d, 'bandwidth');
end

cmd{1}.freqfcn    = 'wpass';
cmd{1}.frequency  = [0 fc-bw/2];
cmd{1}.filtertype = 'lowpass';

cmd{2}.freqfcn    = 'wstop';
cmd{2}.frequency  = [fc+bw/2 getnyquist(d)];
cmd{2}.filtertype = 'lowpass';

% [EOF]
