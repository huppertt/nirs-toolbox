function [Phi, W] = computephasedelay(this, varargin)
%COMPUTEPHASEDELAY Phase Delay of a discrete-time filter

%   Copyright 2007 The MathWorks, Inc.

% This should be private

[Phi, W] = ms_freqresp(this, @phasedelay, @sum, varargin{:});

% [EOF]