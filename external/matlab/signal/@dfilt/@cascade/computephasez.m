function [phi, w] = computephasez(Hd,varargin)
%COMPUTEPHASEZ Phase response of a discrete-time filter.

%   Author: V. Pellissier
%   Copyright 1988-2005 The MathWorks, Inc.

% This should be private

if nargin < 2
    varargin{1} = 8192;
end

% Define the new N-point frequency vector where the frequency response is evaluated
[upn_or_w, upfactor, iswholerange, options, addpoint] = findfreqvector(varargin{:});

% Compute the frequency response (freqz) without forming the TF
[h, w] = freqz(Hd, upn_or_w, options{:});

% Extract phase from frequency response
[phi,w] = extract_phase(h,upn_or_w,iswholerange,upfactor,w,addpoint);

% [EOF]
