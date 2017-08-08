function this = lpisincstopapass(varargin)
%LPISINCSTOPAPASS   Construct a LPISINCSTOPAPASS object.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

this = fspecs.lpisincstopapass;

fsconstructor(this,'Inverse-sinc lowpass',2,2,6,varargin{:});

% [EOF]
