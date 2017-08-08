function [Hout,anum,aden] = ciirxform(Hd,fun,varargin)
%CIIRXFORM Complex IIR Transformations.
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
%   Copyright 1988-2004 The MathWorks, Inc.

warnsv(Hd);

% Perform transformation of poles and zeros
[z,p,k] = zpk(Hd);
[zt,pt,kt,anum,aden] = feval(fun,z,p,k,varargin{:});
zot = leja(zt);
pot = leja(pt);

% Form numerator
numsos = formhalfsos(zot);

% Form denominator
densos = formhalfsos(pot);

% Form sos matrix
sosm = [numsos densos];

Hout = copy(Hd);
Hout.sosmatrix = sosm;

% Distribute kt evenly among the sections
nsec = size(numsos,1);
ktr = kt^(1/nsec);
Hout.ScaleValues = [ktr*ones(nsec,1);1];


%---------------------------------------------------------------------
function hsos = formhalfsos(r)
% Form half the sos matrix from the roots 
l = length(r);
if rem(l,2),
    hsos = zeros((l+1)/2,3);
    hsos(1,1:3) = [1 -r(end) 0];
    for n = 1:(l-1)/2,
        hsos(n+1,1:3) = poly(r(2*n-1:2*n));
    end
else
    hsos = zeros(l/2,3);
    for n = 1:size(hsos,1),
        hsos(n,1:3) = poly(r(2*n-1:2*n));
    end
end

