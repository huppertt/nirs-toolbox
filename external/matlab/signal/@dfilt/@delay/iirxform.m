function [Ht, anum, aden] = iirxform(Ho,fun,varargin)
%IIRXFORM   IIR Transformations.

%   Author(s): R. Losada
%   Copyright 2005-2006 The MathWorks, Inc.

% This should be private

[b, a] = tf(Ho);

[num, den, anum, aden] = feval(fun, b, a, varargin{:});

% Create the transformed and allpass filters
if isreal(num),
    if sign(num(end)) == sign(den(1)),        
        Ht  = dfilt.allpass(den(2:end));
    else
        % Coefficients are negative of each other, account for that
        Ht = cascade(dfilt.scalar(-1),dfilt.allpass(den(2:end)));
    end
else
    % Complex case; cannot use dfilt.allpass
    Ht  = dfilt.df2(num,den);
end

% [EOF]


