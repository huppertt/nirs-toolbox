function [phi, w] = computephasez(Hd, varargin)
%COMPUTEPHASEZ Compute the phasez

%   Author: V. Pellissier, J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

% This should be private

inputs    = freqzinputs(Hd, varargin{:});
[b,a]     = tf(Hd);
[phi,w]   = phasez(b,a,inputs{:});

% [EOF]
