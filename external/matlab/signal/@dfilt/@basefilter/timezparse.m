function [Hd, opts] = timezparse(Hd, varargin)
%TIMEZPARSE Parse the time response inputs

%   Author: J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

% Call the nonmethod parser
[N, Fs] = timezparse(varargin{:});

% Convert the values to opts structure for FVTool
if isempty(N),
    opts.uselength = 'default';
else
    opts.impzlength = N;
    opts.uselength = 'specified';
end

% Set up FVTool's sampling Frequency
if isempty(Fs),
    opts.freqmode = 'on';
else
    opts.freqmode = 'off';
end

% Build the dfiltwfs for FVTool
for indx = 1:length(Hd),
    Hd(indx) = dfilt.dfiltwfs(Hd(indx), Fs);
end

% [EOF]
