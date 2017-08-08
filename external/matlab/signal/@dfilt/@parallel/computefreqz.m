function [h, w] = computefreqz(Hd,varargin)
%COMPUTEFREQZ  Discrete-time filter frequency response.

%   Author: Thomas A. Bryan, J. Schickler
%   Copyright 1988-2006 The MathWorks, Inc.

% This should be private

% Check if all stages have the same overall rate change factor
checkvalidparallel(Hd);

% Parallel structures add the resulting frequency responses
[h, w] = ms_freqresp(Hd, @freqz, @sum, varargin{:});

% [EOF]
