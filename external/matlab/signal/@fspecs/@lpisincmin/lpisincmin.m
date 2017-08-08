function this = lpisincmin(varargin)
%LPISINCMIN   Construct a LPISINCMIN object.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

this = fspecs.lpisincmin;

fsconstructor(this,'Inverse-sinc Lowpass',2,2,7,varargin{:});

% [EOF]
