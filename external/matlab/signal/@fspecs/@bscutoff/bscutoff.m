function h = bscutoff(varargin)
%BSCUTOFF   Construct a BSCUTOFF object.
%   H = BSCUTOFF(N,Fcutoff1,Fcutoff2,Fs) constructs a bandstop filter
%   specifications object with cutoff frequencies.
%
%   N is the filter order and must be an even positive integer.
%
%   Fcutoff1 is the lower cutoff frequency and must be a positive scalar
%   between 0 and 1 if no sampling frequency is specified or between 0 and
%   Fs/2 if a sampling frequency Fs is specified.
%
%   Fcutoff2 is the higher cutoff frequency and must be a positive scalar
%   between 0 and 1 if no sampling frequency is specified or between 0 and
%   Fs/2 if a sampling frequency Fs is specified.
%
%   Fs is the sampling frequency. If Fs is not specified, normalized
%   frequency is assumed. If Fs is specified, it must be a positive scalar.

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

h = fspecs.bscutoff;
constructor(h,varargin{:});
h.ResponseType = 'Bandstop with cutoff';
% [EOF]
