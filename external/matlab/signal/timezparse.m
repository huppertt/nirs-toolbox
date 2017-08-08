function [N, Fs] = timezparse(varargin)
%TIMEZPARSE Parser for the time responses

%   Author: J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

N  = [];
Fs = [];

if nargin > 0, N  = varargin{1}; end
if nargin > 1, Fs = varargin{2}; end

% [EOF]
