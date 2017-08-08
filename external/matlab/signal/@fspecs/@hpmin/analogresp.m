function ha = analogresp(h)
%ANALOGRESP   

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

% Compute analog frequencies
wp = cot(pi*h.Fpass/2);
ws = cot(pi*h.Fstop/2);

% Construct analog specs object
ha = fspecs.alpmin(wp,ws,h.Apass,h.Astop);

% [EOF]
