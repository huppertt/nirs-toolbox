function ha = analogresp(h)
%ANALOGRESP   Compute analog response object.

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

% Compute analog frequency
c = cparam(h);
wc = (c - cos(pi*h.Fcutoff2))/sin(pi*h.Fcutoff2);

% Construct analog specs object
ha = fspecs.alpcutoff(h.FilterOrder,wc);


% [EOF]
