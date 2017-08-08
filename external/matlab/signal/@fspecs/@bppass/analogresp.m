function ha = analogresp(h)
%ANALOGRESP   Compute analog response object.

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

% Compute analog frequency
c = cparam(h);
wp1 = (c-cos(pi*h.Fpass1))/sin(pi*h.Fpass1);
wp2 = (c-cos(pi*h.Fpass2))/sin(pi*h.Fpass2);
wp = min(abs([wp1,wp2]));


% Construct analog specs object
ha = fspecs.alppass(h.FilterOrder,wp,h.Apass);


% [EOF]
