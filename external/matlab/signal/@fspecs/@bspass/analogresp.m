function ha = analogresp(h)
%ANALOGRESP   Compute analog response object.

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

% Compute analog frequency
c = cparam(h);
wp1 = sin(pi*h.Fpass1)/(cos(pi*h.Fpass1)-c);
wp2 = sin(pi*h.Fpass2)/(cos(pi*h.Fpass2)-c);
wp = min(abs([wp1,wp2]));

% Construct analog specs object
ha = fspecs.alppass(h.FilterOrder,wp,h.Apass);


% [EOF]
