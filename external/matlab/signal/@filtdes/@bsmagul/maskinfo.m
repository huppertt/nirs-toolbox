function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

astop  = get(d, 'DstopUpper');
astopp = -astop;

cmd{1}.magfcn     = 'pass';
cmd{1}.amplitude  = [d.Dpass1Upper -d.Dpass1Lower]/2;
cmd{1}.filtertype = 'lowpass';
cmd{1}.magunits   = 'linear';
cmd{1}.tag        = 'pass1';
cmd{1}.astop      = astopp;

cmd{2}.magfcn     = 'stop';
cmd{2}.amplitude  = [d.DstopUpper -d.DstopLower];
cmd{2}.filtertype = 'bandstop';
cmd{2}.magunits   = 'linear';

cmd{3}.magfcn     = 'pass';
cmd{3}.amplitude  = [d.Dpass2Upper -d.Dpass2Lower]/2;
cmd{3}.filtertype = 'highpass';
cmd{3}.magunits   = 'linear';
cmd{3}.tag        = 'pass2';
cmd{3}.astop      = astopp;

% [EOF]
