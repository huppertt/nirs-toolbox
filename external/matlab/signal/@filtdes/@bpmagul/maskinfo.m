function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

astop1 = get(d, 'Dstop1Upper');
apass  = get(d, 'DpassUpper');
astop2 = get(d, 'Dstop2Upper');
astopp = -max(astop1, astop2);

cmd{1}.magfcn     = 'stop';
cmd{1}.amplitude  = [astop1 -d.Dstop1Lower];
cmd{1}.filtertype = 'highpass';
cmd{1}.magunits   = 'linear';
cmd{1}.tag        = 'stop1';

cmd{2}.magfcn     = 'pass';
cmd{2}.amplitude  = [d.DpassUpper -d.DpassLower]/2;
cmd{2}.filtertype = 'bandpass';
cmd{2}.magunits   = 'linear';
cmd{2}.astop      = astopp;

cmd{3}.magfcn     = 'stop';
cmd{3}.amplitude  = [astop2 -d.Dstop2Lower];
cmd{3}.filtertype = 'lowpass';
cmd{3}.magunits   = 'linear';
cmd{3}.tag        = 'stop2';
% [EOF]
