function [num, den] = thistf(Hd)
%THISTF  Convert to transfer function.
%   [NUM,DEN] = THISTF(Hd) converts discrete-time filter Hd to numerator and
%   denominator vectors.
%
%   See also DFILT.

%   Author: Thomas A. Bryan
%   Copyright 1988-2006 The MathWorks, Inc.

% This should be private

% Check if all stages have the same overall rate change factor
checkvalidparallel(Hd);

% The algorithm is equivalent to adding rational numbers:
%    b1/a1 + b2/a2 = (b1*a2 + a1*b2)/(a1*a2)
% where convolution takes the place of *.

if isempty(Hd.Stage)
    num = [];
    den = [];
    return
end

[num,den] = tf(Hd.Stage(1));

for k=2:length(Hd.Stage)
    [b,a] = tf(Hd.Stage(k));

    % (b1*a2 + a1*b2)
    num1 = conv(num,a);
    num2 = conv(den,b);
    if max(abs(num2)) == 0
        num = num1;
    else
        [num1,num2] = eqtflength(num1,num2);
        num = num1 + num2;
    end

    % (a1*a2)
    den = conv(den,a);
end

num = removetrailzeros(num);
den = removetrailzeros(den);
