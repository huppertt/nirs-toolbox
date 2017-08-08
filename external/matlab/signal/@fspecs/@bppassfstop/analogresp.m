function ha = analogresp(h)
%ANALOGRESP   Compute analog response object.

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

% Compute analog frequency
c = cparam(h);
wp1 = (c-cos(pi*h.Fpass1))/sin(pi*h.Fpass1);
wp2 = (c-cos(pi*h.Fpass2))/sin(pi*h.Fpass2);
wp = min(abs([wp1,wp2]));

ws1 = (c-cos(pi*h.Fstop1))/sin(pi*h.Fstop1);
ws2 = (c-cos(pi*h.Fstop2))/sin(pi*h.Fstop2);
ws = min(abs([ws1,ws2]));

% Construct analog specs object
ha = fspecs.alppassfstop(h.FilterOrder,wp,ws,h.Apass);


% [EOF]
