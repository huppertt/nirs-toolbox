function [h, w] = computefreqz(Hd,varargin)
%COMPUTEFREQZ Compute the freqz

%   Author: Thomas A. Bryan, J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

% This should be private

% Cascades multiply the resulting frequency responses
[h, w] = ms_freqresp(Hd, @freqz, @prod, varargin{:});

% [EOF]
