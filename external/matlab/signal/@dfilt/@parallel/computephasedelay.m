function [Phi, W] = computephasedelay(this, varargin)
%COMPUTEPHASEDELAY Phase Delay of a discrete-time filter

%   Copyright 2007 The MathWorks, Inc.

% This should be private

% Check if all stages have the same overall rate change factor
checkvalidparallel(this);

inputs = freqzinputs(this, varargin{:});
[b,a]  = tf(this);
[Phi,W]  = phasedelay(b,a,inputs{:});

% [EOF]

