function [G,w] = computegrpdelay(Hd,varargin)
%COMPUTEGRPDELAY   

%   Author(s): R. Losada
%   Copyright 2005-2006 The MathWorks, Inc.


% This should be private

% Check if all stages have the same overall rate change factor
checkvalidparallel(Hd);

inputs = freqzinputs(Hd, varargin{:});
[b,a]  = tf(Hd);
[G,w]  = grpdelay(b,a,inputs{:});

% [EOF]
