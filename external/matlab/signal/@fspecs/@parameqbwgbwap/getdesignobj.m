function designobj = getdesignobj(this, str)
%GETDESIGNOBJ   Get the design object.

%   Author(s): J. Schickler
%   Copyright 2006 The MathWorks, Inc.


%#function fdfmethod.butterparameqminbwap
%#function fdfmethod.cheby1parameqminbwap
designobj.butter     = 'fdfmethod.butterparameqminbwap';
designobj.cheby1     = 'fdfmethod.cheby1parameqminbwap';


if nargin > 1
    designobj = designobj.(str);
end

% [EOF]
