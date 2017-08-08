function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

cmd = base_maskinfo(hObj, d);

specobjs = get(hObj, 'SpecObjs');

fcmd = maskinfo(specobjs(1), d);
mcmd = maskinfo(specobjs(2), d);

fcmd{2} = setstructfields(fcmd{2}, mcmd{1});

fcmd{1} = fcmd{2}; % Copy defaults from passband

fcmd{1}.frequency = [1 1]*getnyquist(d) - fliplr(fcmd{2}.frequency);
fcmd{1}.magfcn    = 'stop';

if isdb(d),
    val = fcmd{2}.amplitude;
    val = (10^(val/20) - 1)/(10^(val/20) + 1); % perform the tolinear passband conversion
    val = -20*log10(val);                      % perform the todb stopband conversion
    
    fcmd{1}.amplitude = val;
else
    fcmd{2}.astop     = -fcmd{2}.amplitude;
end

cmd.bands    = fcmd;

% [EOF]
