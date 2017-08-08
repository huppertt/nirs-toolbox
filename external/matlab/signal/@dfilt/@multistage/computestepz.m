function [y, t] = computetstepz(Hd, varargin)
%COMPUTESTEPZ Compute the stepz for the filter

%   Author(s): J. Schickler, R. Losada, T. Bryan
%   Copyright 1988-2004 The MathWorks, Inc.

% This should be private

[N, Fs] = timezparse(varargin{:});
if isempty(N),  N  = impzlength(Hd); end
if isempty(Fs), Fs = 1;              end

% Form time vector
t = (0:N-1)'./Fs;

% Form input vector
x = ones(size(t));

% Get current value of PersistentMemory
resetval = Hd.PersistentMemory;
Hd.PersistentMemory = false;
states = exportstates(Hd);
y = filter(Hd,x);
importstates(Hd,states);
% Set PersistentMemory back to what it was
Hd.PersistentMemory = resetval;

% [EOF]
