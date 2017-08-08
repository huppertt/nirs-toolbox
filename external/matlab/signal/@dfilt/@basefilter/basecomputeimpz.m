function [y,T] = basecomputeimpz(this,varargin)
%BASECOMPUTEIMPZ   

%   Author(s): R. Losada
%   Copyright 2003 The MathWorks, Inc.

% This should be private
Hd = double(this);
[N, Fs] = timezparse(varargin{:});
if isempty(N),  N  = impzlength(Hd); end
if isempty(Fs), Fs = 1;              end

T = (0:N-1)'/Fs;
x = [1;zeros(N-1,1)];

% Filter to compute impz
y = filter(Hd,x);



% [EOF]
