function defaultmethod = getdefaultmethod(this)
%GETDEFAULTMETHOD Get the defaultmethod.

%   Copyright 2005-2011 The MathWorks, Inc.

switch this.Specification,
    case 'N,F,A'
        defaultmethod = 'freqsamp';
    case 'F,A,R'
        defaultmethod = 'equiripple'; 
    case 'Nb,Na,F,A'
        defaultmethod = 'iirlpnorm';
    case 'N,B,F,A'
        defaultmethod = 'equiripple';
    case 'N,B,F,A,C'
        defaultmethod = 'equiripple';        
    case 'B,F,A,R'
        defaultmethod = 'equiripple';        
    case 'Nb,Na,B,F,A'
        defaultmethod = 'iirlpnorm';
    otherwise,
        error(message('signal:fdesign:arbmag:getdefaultmethod:InternalError'));
end

% [EOF]
