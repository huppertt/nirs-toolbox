function ha = analogresp(h)
%ANALOGRESP   Compute analog response object.

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

% Compute analog frequency
wp = tan(pi*h.Fpass/2);
ws = tan(pi*h.Fstop/2);

% Construct analog specs object
ha = fspecs.alppassfstop(h.FilterOrder,wp,ws,h.Apass);


% [EOF]
