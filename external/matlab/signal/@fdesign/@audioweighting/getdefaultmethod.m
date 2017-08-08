function defaultmethod = getdefaultmethod(this)
%GETDEFAULTMETHOD   Get the defaultmethod.

%   Copyright 2009 The MathWorks, Inc.

switch lower(this.Specification)
    case {'wt'}
        defaultmethod = 'freqsamp';
    case {'wt,class'}
        defaultmethod = 'ansis142';
    otherwise,
        error(message('signal:fdesign:audioweighting:getdefaultmethod:InternalError'));
end


% [EOF]
