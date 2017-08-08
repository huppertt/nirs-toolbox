function Ht = firxform(Hd,fun,varargin)
%FIRXFORM FIR Transformations

%   Author(s): R. Losada, J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

warnsv(Hd);

Hr = reffilter(Hd);

sosm = get(Hr, 'sosMatrix');

for indx = 1:size(sosm, 1)
    sosm(indx,1:3) = feval(fun, sosm(indx,1:3), varargin{:});
end

% Create the transformed filter
Ht = copy(Hd);
arith = Ht.Arithmetic; % Cache setting
Ht.Arithmetic = 'double';
Ht.sosMatrix = sosm;
Ht.Arithmetic = arith; % Reset arithmetic


% [EOF]
