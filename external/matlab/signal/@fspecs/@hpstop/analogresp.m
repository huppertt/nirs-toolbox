function ha = analogresp(h)
%ANALOGRESP   Compute analog response object.

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

% Compute analog frequency
ws = cot(pi*h.Fstop/2);

% Construct analog specs object
ha = fspecs.alpstop(h.FilterOrder,ws,h.Astop);


% [EOF]
