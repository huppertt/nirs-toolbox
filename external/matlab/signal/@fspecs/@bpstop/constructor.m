function constructor(h,varargin)
%CONSTRUCTOR   

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.


respstr = h.ResponseType;
fstart = 2;
fstop = 3;
nargsnoFs = 4;
fsconstructor(h,respstr,fstart,fstop,nargsnoFs,varargin{:});


% [EOF]
