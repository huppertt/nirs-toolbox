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

warnsv(Hd);

% Form new SOS matrix
Href = reffilter(Hd);
sm   = Href.sosMatrix;

% Evaluate transformation on each section
sosm = [];
anum = [];
aden = [];
for k = 1:nsections(Href),
    [num,den,anum,aden] = feval(fun,sm(k,1:3),sm(k,4:6),varargin{:});
    s    = tf2sos(num,den);
    sosm = [sosm;s];
end

% Force leading numerator to one
sv = sosm(:,1);
sosm(:,1:3) = sosm(:,1:3)./repmat(sosm(:,1),1,3);

Ht = copy(Hd);
arith = Ht.Arithmetic; % Cache setting
Ht.Arithmetic = 'double';
Ht.sosmatrix = sosm;
Ht.scaleValues(1:size(sv,1)) = Ht.ScaleValues(1:size(sv,1)).*sv;
Ht.Arithmetic = arith; % Reset arithmetic



