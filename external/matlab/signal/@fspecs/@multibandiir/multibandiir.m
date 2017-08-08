function this = multiband(varargin)
%MULTIBAND   Construct a MULTIBAND object.

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

this = fspecs.multibandiir;

respstr = 'Multi-Band Arbitrary Magnitude IIR';
fstart = 1;
fstop = 1;
nargsnoFs = 3;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
