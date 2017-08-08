function [hz, wz, phiz, opts] = computezerophase(Hd, varargin)
%COMPUTEZEROPHASE

%   Author: V. Pellissier, J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

% This should be private

inputs = freqzinputs(Hd, varargin{:});
[b,a]  = tf(Hd);
[hz,wz,phiz,opts] = zerophase(b,a,inputs{:});

% [EOF]
