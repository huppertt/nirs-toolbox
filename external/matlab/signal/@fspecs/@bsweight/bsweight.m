function this = bsweight(varargin)
%BSWEIGHT   Construct a BSWEIGHT object.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

this = fspecs.bsweight;

respstr = 'Bandstop';
fstart = 2;
fstop = 5;
nargsnoFs = 8;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
