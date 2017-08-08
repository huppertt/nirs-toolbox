function ha = analogresp(h)
%ANALOGRESP   Compute analog response object.

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

% Compute analog frequency
wp = cot(pi*h.Fpass/2);

% Construct analog specs object
ha = fspecs.alppass(h.FilterOrder,wp,h.Apass);


% [EOF]
