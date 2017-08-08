function [Ht,anum,aden] = iirxform(Hd,fun,varargin)
%IIRXFORM IIR Transformations
%
%   Inputs:
%       Hd - Handle to original filter
%       fun - function handle to transformation
%
%   Outputs:
%       Hout - Transformed filter
%       anum - Allpass numerator
%       aden - Allpass denominator

%   Author(s): R. Losada
%   Copyright 1988-2006 The MathWorks, Inc.

c = coefficients(Hd);
ns = nsections(Hd);
allpassden = cell(ns,1);
isre = true;
isneg = false;
for k = 1:ns,
    a = [1 c{k}];
    b = fliplr(a);
    [num,den,anum,aden] = feval(fun,b,a,varargin{:});
    if isreal(num),
        if sign(num(end)) ~= sign(den(1)),
            isneg = true;
        end
        allpassden{k} = den(2:end);
    else
        isre = false;
        allpassden{k} = dfilt.df2(num,den);
    end
end

if isre,
    if isneg,
        Ht = cascade(dfilt.scalar(-1),feval(class(Hd),allpassden{:}));
    else        
        Ht = feval(class(Hd),allpassden{:});
    end
else
    % Complex case; cannot use dfilt.allpass
    if ns > 1,
        Ht  = cascade(allpassden{:});
    else
        Ht = allpassden{1};
    end
end




