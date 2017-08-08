function defaultmethod = getdefaultmethod(this)
%GETDEFAULTMETHOD   Get the defaultmethod.

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

switch this.Specification,
    case 'N,F0',
        defaultmethod = 'butter';
    otherwise,
        error(message('signal:fdesign:octave:getdefaultmethod:InternalError'));
end


% [EOF]
