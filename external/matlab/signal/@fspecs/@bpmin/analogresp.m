function ha = analogresp(h)
%ANALOGRESP   

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

% Compute analog frequencies
c = cparam(h);
wp = abs((c-cos(pi*h.Fpass2))/sin(pi*h.Fpass2));
ws1 = (c-cos(pi*h.Fstop1))/sin(pi*h.Fstop1);
ws2 = (c-cos(pi*h.Fstop2))/sin(pi*h.Fstop2);
ws = min(abs([ws1,ws2]));

% Construct analog specs object
ha = fspecs.alpmin(wp,ws,h.Apass,max(h.Astop1,h.Astop2));

% [EOF]
