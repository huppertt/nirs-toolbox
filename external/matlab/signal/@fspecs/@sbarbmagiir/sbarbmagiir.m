function this = sbarbmagiir(varargin)
%SBARBMAGIIR   Construct a SBARBMAGIIR object.

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

this = fspecs.sbarbmagiir;
respstr = 'Single-Band Arbitrary Magnitude IIR';
fstart = 1;
fstop = 1;
nargsnoFs = 4;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
