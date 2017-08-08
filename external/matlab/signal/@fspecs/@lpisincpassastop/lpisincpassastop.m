function this = lpisincpassastop(varargin)
%LPISINCPASSASTOP   Construct a LPISINCPASSASTOP object.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

this = fspecs.lpisincpassastop;

fsconstructor(this,'Inverse-sinc lowpass',2,2,6,varargin{:});

% [EOF]
