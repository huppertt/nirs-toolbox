function [Hd, str] = singlesection(this)
%SINGLESECTION   Convert to a single section.

%   Author(s): J. Schickler
%   Copyright 2004 The MathWorks, Inc.

[b, a] = tf(this);

str = {'[b, a] = tf(Hd);'};

if isfir(this)
    
    % If the SOS filter is actually an FIR, we put its coefficients into a
    % DFFIR object because it is the "simplest" structure.
    Hd = dfilt.dffir(b/a);

    str = {str{:}, 'Hd     = dfilt.dffir(b/a);'};
else

    fs = class(this);
    fs(end-2:end) = [];

    Hd = feval(fs, b, a);
    str = {str{:}, sprintf('Hd     = %s(b, a);', fs)};
end

% [EOF]
