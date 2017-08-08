function constructor(h,varargin)
%CONSTRUCTOR   

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

respstr = h.ResponseType;
fstart = 2;
fstop = 2;
nargsnoFs = 3;
fsconstructor(h,respstr,fstart,fstop,nargsnoFs,varargin{:});


% [EOF]
