function this = hpisincstopapass(varargin)
%HPISINCSTOPAPASS Construct a HPISINCSTOPAPASS object.

%   Copyright 2011 The MathWorks, Inc.

this = fspecs.hpisincstopapass;

fsconstructor(this,'Inverse-sinc highpass',2,2,6,varargin{:});

% [EOF]
