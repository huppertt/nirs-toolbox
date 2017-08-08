function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

fc = getnyquist(d)/get(d, 'Band');

if isprop(d, 'TransitionMode'),
    
    if strcmpi(get(d, 'TransitionMode'), 'Rolloff'),
        bw = 2*fc*get(d, 'Rolloff');
    else
        bw = get(d, 'bandwidth');
    end
else
    bw = 0;
end

cmd{1}.frequency = [0 fc-bw/2];
cmd{2}.frequency = [fc+bw/2 getnyquist(d)];

if bw == 0,
    cmd{2}.drawfreqbars = false;
end

cmd{1}.filtertype = 'lowpass';
cmd{2}.filtertype = 'lowpass';

cmd{1}.freqfcn    = 'wpass';
cmd{2}.freqfcn    = 'wstop';

% [EOF]
