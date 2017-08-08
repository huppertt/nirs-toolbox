function [Hd, opts] = freqzparse(Hd, varargin)
%FREQZPARSE Parser for freqz

%   Author: J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

% Call the nonmethod parser
[NorW, unitcircle, Fs] = freqzparse(varargin{:});

% Convet the values to opts structure for FVTool
if length(NorW) > 1,
    opts.freqvec    = NorW;
    opts.unitcircle = 4;
else
    opts.nfft = NorW;
    opts.unitcircle = find(strcmpi(unitcircle, {'half','whole'}));
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