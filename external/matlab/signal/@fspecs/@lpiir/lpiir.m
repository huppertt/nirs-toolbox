function this = lpiir(varargin)
%LPIIR   Construct a LPIIR object.

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

this = fspecs.lpiir;

respstr = 'Lowpass';
fstart = 3;
fstop = 4;
nargsnoFs = 6;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
