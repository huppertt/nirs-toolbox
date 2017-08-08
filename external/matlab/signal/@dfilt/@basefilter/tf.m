function [num,den] = tf(Hb)
%TF  Convert to transfer function.
%   [NUM,DEN] = TF(Hb) converts discrete-time filter Hb to numerator and
%   denominator vectors.
%
%   See also DFILT.   

%   Author: V. Pellissier, J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

Hd         = dispatch(Hb(1));
[num, den] = thistf(Hd);

for indx = 2:length(Hb),
    Hd = dispatch(Hb(indx));
    [n,d] = thistf(Hd);
    
    num = lclpadwzeros(num, n);
    den = lclpadwzeros(den, d);
end

% ------------------------------------------------------
function num = lclpadwzeros(num, n)

if length(num) < length(n),
    num(end, length(n)) = 0;
elseif length(num) > length(n),
    n(length(num)) = 0;
end

num(end+1,:) = n;

% [EOF]
