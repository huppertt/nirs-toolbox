function logreport = getlogreport(this)
%GETLOGREPORT   Get the logreport.

%   Author(s): V. Pellissier
%   Copyright 2005-2006 The MathWorks, Inc.

logreport = this.filterquantizer.loggingreport;
if isempty(logreport),
    w = warning('on');
    warning(message('signal:dfilt:abstractfilter:getlogreport:LoggingOff', '''Arithmetic''', '''fixed'''));
    warning(w);
else
    logreport = copy(logreport);
end

% [EOF]
