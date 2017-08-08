function this = hpiir(varargin)
%HPIIR   Construct a HPIIR object.

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

this = fspecs.hpiir;

respstr = 'Highpass';
fstart = 3;
fstop = 4;
nargsnoFs = 6;
fsconstructor(this,respstr,fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
