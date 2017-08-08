function [h, w] = ms_freqresp(Hd, fcn, cfcn, varargin)
%MS_FREQRESP Compute the MultiSection Frequency Response
%   MS_FREQRESP(Hd, FCN, CFCN) Compute the multisection frequency response
%   specified by the function handle FCN on the filter Hd.  Use CFCN to
%   combine the responses.  CFCN should be a function handle to PROD or SUM
%   or a function with the same inputs.

%   Author: J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

% This should be private

[h, w] = feval(fcn, Hd.Stage, varargin{:});

[r, c] = size(w);
if r > c, dim = 2;
else,     dim = 1; end

h = feval(cfcn, h, dim);

% [EOF]
