function this = bsiir(varargin)
%BSIIR   Construct a BSIIR object.

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

this = fspecs.bsiir;
respstr = 'Bandstop';
fstart = 3;
fstop = 5;
nargsnoFs = 8;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
