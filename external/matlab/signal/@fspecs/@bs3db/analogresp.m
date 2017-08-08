function ha = analogresp(h)
%ANALOGRESP   Compute analog response object.

%   Author(s): R. Losada
%   Copyright 1988-2005 The MathWorks, Inc.

% Compute analog frequency
c = cparam(h);
wc = abs(sin(pi*h.F3dB2)/(cos(pi*h.F3dB2) - c));

% Construct analog specs object
ha = fspecs.alpcutoff(h.FilterOrder,wc);


% [EOF]
