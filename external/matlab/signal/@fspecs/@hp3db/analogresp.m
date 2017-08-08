function ha = analogresp(h)
%ANALOGRESP   Compute analog response object.

%   Author(s): R. Losada
%   Copyright 1988-2005 The MathWorks, Inc.

% Compute analog frequency
wc = cot(pi*h.F3db/2);

% Construct analog specs object
ha = fspecs.alpcutoff(h.FilterOrder,wc);


% [EOF]
