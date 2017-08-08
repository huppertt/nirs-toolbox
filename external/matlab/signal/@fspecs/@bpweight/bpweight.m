function this = bpweight(varargin)
%BPWEIGHT   Construct a BPWEIGHT object.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

this = fspecs.bpweight;

respstr = 'Bandpass';
fstart = 2;
fstop = 5;
nargsnoFs = 8;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
