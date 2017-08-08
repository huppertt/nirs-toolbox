function this = hpisincpassastop(varargin)
%HPISINCPASSASTOP Construct a HPISINCPASSASTOP object.

%   Copyright 2011 The MathWorks, Inc.

this = fspecs.hpisincpassastop;

fsconstructor(this,'Inverse-sinc highpass',2,2,6,varargin{:});

% [EOF]
