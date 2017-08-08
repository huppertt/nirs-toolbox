function ha = analogresp(h)
%ANALOGRESP   

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

% Compute analog frequencies
c = cparam(h);
wp = abs(sin(pi*h.Fpass2)/(cos(pi*h.Fpass2)-c));
ws1 = sin(pi*h.Fstop1)/(cos(pi*h.Fstop1)-c);
ws2 = sin(pi*h.Fstop2)/(cos(pi*h.Fstop2)-c);
ws = min(abs([ws1,ws2]));

% Construct analog specs object
ha = fspecs.alpmin(wp,ws,min(h.Apass1,h.Apass2),h.Astop);

% [EOF]
