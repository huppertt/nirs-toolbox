function constructor(h,varargin)
%CONSTRUCTOR   

%   Copyright 2008-2011 The MathWorks, Inc.

respstr = h.ResponseType;
fstart = 3;
fstop = 3;
nargsnoFs = 3;
fsconstructor(h,respstr,fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
