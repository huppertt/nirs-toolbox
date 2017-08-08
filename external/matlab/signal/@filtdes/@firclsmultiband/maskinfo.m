function cmd = maskinfo(h, d)
%MASKINFO

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

cmd = base_maskinfo(h, d);

f = get(d, 'FrequencyVector');
a = get(d, 'MagnitudeVector');

band{1}.magfcn     = 'aline';
band{1}.frequency  = [f(1) reshape([f(2:end-1)' f(2:end-1)']', 1,(length(f)-2)*2), f(end)];
band{1}.amplitude  = [a(1) a(1) reshape([a(2:end-1)' a(2:end-1)']', 1,(length(a)-2)*2), a(end) a(end)];
band{1}.properties = {'Color', [0 0 0]};

cmd.bands = band;

% [EOF]
