function cmd = maskinfo(h, d)
%MASKINFO

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

nbands = get(d, 'Order');

bw = getbandwidth(d);

fn = getnyquist(d);

if nbands == 1,
    cmd{1}.frequency  = [fn-bw/2 fn];
    cmd{1}.filtertype = 'highpass';
    cmd{1}.freqfcn    = 'wpass';
else
    sep = fn/(nbands/2);
    
    for indx = 1:2:nbands-1
        cmd{(indx+1)/2}.frequency  = [sep/2-bw/2+sep*((indx+1)/2-1) sep/2+bw/2+sep*((indx+1)/2-1)];
        cmd{(indx+1)/2}.filtertype = 'bandpass';
        cmd{(indx+1)/2}.freqfcn    = 'wpass';
    end
    
    if rem(nbands, 2),
        cmd{end+1} = cmd{end};
        cmd{end}.frequency  = [fn-bw/2 fn];
        cmd{end}.filtertype = 'highpass';
    end
end

% [EOF]
