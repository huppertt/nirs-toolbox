function defaultmethod = getdefaultmethod(this)
%GETDEFAULTMETHOD   Get the defaultmethod.

%   Copyright 2005-2010 The MathWorks, Inc.

switch lower(this.Specification)
    case 'n,f,gd',
        defaultmethod = 'iirlpnorm';
    case 'n,b,f,gd', 
        defaultmethod = 'iirlpnorm';
  otherwise
      error(message('signal:fdesign:arbgrpdelay:getdefaultmethod:InternalError'));      
end

% [EOF]
