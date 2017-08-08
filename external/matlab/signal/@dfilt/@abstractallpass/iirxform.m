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


[b , a]= tf(Hd);
[num,den,anum,aden] = feval(fun,b,a,varargin{:});
if isreal(num),
    if sign(num(end)) == sign(den(1)),
        Ht = feval(class(Hd),den(2:end));
    else
        % Coefficients are negative of each other, account for that
        Ht = cascade(dfilt.scalar(-1),feval(class(Hd),den(2:end)));
    end
else
    % Complex case; cannot use dfilt.allpass
    Ht  = dfilt.df2(num,den);
end





