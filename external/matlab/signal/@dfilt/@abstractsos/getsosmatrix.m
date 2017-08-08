function sosm = getsosmatrix(hObj, sosm)
%GETSOSMATRIX Get the sosmatrix from the object.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

sosm = getsosmatrix(hObj.filterquantizer, get(hObj, 'privNum'), get(hObj, 'privDen'));

% [EOF]
