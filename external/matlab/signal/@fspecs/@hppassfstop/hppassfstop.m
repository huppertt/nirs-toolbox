function h = hppassfstop(varargin)
%HPPASSFSTOP   Construct a HPPASSFSTOP object.
%   H = HPPASSFSTOP(N,Fstop,Fpass,Apass,Fs) constructs a highpass filter
%   specifications object with passband-edge specs.
%
%   N is the filter order and must be a positive integer.
%
%   Fstop is the stopband-edge frequency and must be a positive scalar
%   between 0 and 1 if no sampling frequency is specified or between 0 and
%   Fs/2 if a sampling frequency Fs is specified.
%
%   Fpass is the passband-edge frequency and must be a positive scalar
%   greater than Fstop and between 0 and 1 if no sampling frequency is
%   specified or between 0 and Fs/2 if a sampling frequency Fs is
%   specified.
%
%   Apass is the maximum passband deviation and it must be a positive
%   scalar.
%
%   Fs is the sampling frequency. If Fs is not specified, normalized
%   frequency is assumed. If Fs is specified, it must be a positive scalar.

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

% Override defaults inherited from lowpass

if nargin < 1,
    varargin{1} = 10;
end

if nargin < 2,
    varargin{2} = .4;
end

if nargin < 3,
    varargin{3} = .6;
end

h = fspecs.hppassfstop;
respstr = 'Highpass with passband-edge specifications and stopband frequency';
fstart = 2;
fstop = 3;
nargsnoFs = 4;
fsconstructor(h,respstr,fstart,fstop,nargsnoFs,varargin{:});


% [EOF]
