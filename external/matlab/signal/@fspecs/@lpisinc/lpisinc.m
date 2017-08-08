function this = lpisinc(varargin)
%LPISINC   Construct a LPISINC object.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

this = fspecs.lpisinc;

fsconstructor(this,'Inverse-sinc Lowpass',2,2,5,varargin{:});

% [EOF]
