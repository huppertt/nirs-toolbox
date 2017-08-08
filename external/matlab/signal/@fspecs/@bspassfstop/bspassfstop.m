function h = bspassfstop(varargin)
%BSPASSFSTOP   Construct a BSPASSFSTOP object.
%   H = BSPASSFSTOP(N,Fstop1,Fpass1,Fpass2,Fstop2,Apass,Fs) constructs a
%   bandstop filter specifications object with passband-edge specifications
%   and stopband frequencies.
%
%   N is the filter order and must be an even positive integer.
%
%   Fstop1 is the lower stopband-edge frequency and must be a positive
%   scalar between 0 and 1 if no sampling frequency is specified or between
%   0 and Fs/2 if a sampling frequency Fs is specified.
%
%   Fpass1 is the lower passband-edge frequency and must be a positive
%   scalar between 0 and 1 if no sampling frequency is specified or between
%   0 and Fs/2 if a sampling frequency Fs is specified.
%
%   Fpass2 is the higher passband-edge frequency and must be a positive
%   scalar larger than Fpass1 and between 0 and 1 if no sampling frequency
%   is specified or between 0 and Fs/2 if a sampling frequency Fs is
%   specified.
%
%   Fstop2 is the higher stopband-edge frequency and must be a positive
%   scalar between 0 and 1 if no sampling frequency is specified or between
%   0 and Fs/2 if a sampling frequency Fs is specified.
%
%   Apass is the maximum passband deviation and it must be a positive
%   scalar.
%
%   Fs is the sampling frequency. If Fs is not specified, normalized
%   frequency is assumed. If Fs is specified, it must be a positive scalar.

%   Author(s): R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.

h = fspecs.bspassfstop;
respstr = 'Bandstop with passband-edge specifications and stopband frequencies.';
fstart = 2;
fstop = 5;
nargsnoFs = 6;

% Fix the frequency specifications
if nargin < 5
    varargin{5} = .65;
    if nargin < 4
        varargin{4} = .55;
        if nargin < 3
            varargin{3} = .45;
            if nargin < 2
                varargin{2} = .35;
                if nargin < 1
                    varargin{1} = 10;
                end
            end
        end
    end
end

fsconstructor(h,respstr,fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
