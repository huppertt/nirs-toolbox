function [H,W] = getdata(this,varargin)
%GETDATA   Get the complex data and frequency from the object.
%  [H,W] = GETDATA(H,ISPSD,ISDB,ISNORMALIZED,FREQRANGE,CENTERDC) returns
%  the data and frequencies from the object H, modified according to the
%  inputs ISPSD, ISDB, ISNORMALIZED, FREQRANGE, and CENTERDC.
%
%  Inputs:
%    this         - handle to data object.
%    ispsd        - boolean indicating if data should be returned as a
%                   power spectral density.
%    isdb         - boolean indicating if data should be in dB scale
%                   (default = false).
%    isnormalized - boolean indicating if frequency should be normalized
%                   (default = state of NormalizedFrequency).
%    freqrange    - string indicating if spectrum should be calculated over
%                   half or the whole Nyquist interval.
%                   The possible options are:
%                        'half'  - convert spectrum to half or onesided
%                        'whole' - convert spectrum to whole or twosided
%   centerdc      - shift the spectrum so that 0 (the DC value) is in the
%                   center of the frequency grid.

%   Author(s): P. Pacheco
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(1,6,nargin,'struct'));

% Set default values and parse inputs.
[ispsd,isdb,isnormalized,freqrange,centerdc] = parseinputs(this,varargin{:});

% Cache object property values.
Fs = this.getfs;
normFlag = this.NormalizedFrequency;

% Define a boolean flag representing the frequency range requested.
ishalfrange = false;
if strcmpi(freqrange,'half'),
    ishalfrange = true;
end

% Calculate the response for the frequency range selected by the user.
[H,W] = computeresp4freqrange(this,ishalfrange,ispsd,isnormalized,centerdc);

if isdb,
    H = convert2db(this,H);
end

%--------------------------------------------------------------------------
function [ispsd,isdb,isnormalized,freqrange,centerdc] = parseinputs(this,varargin)
%PARSEINPUTS   Set default values and parse the input argument list.

% Defaults
ispsd        = isdensity(this);
isdb         = false;
isnormalized = this.NormalizedFrequency;
freqrange    = 'whole';
centerdc     = false;

if ishalfnyqinterval(this),
    freqrange = 'half';
end

% Parse inputs
if nargin >= 2,
    ispsd = varargin{1};
    if nargin >= 3,
        isdb = varargin{2};
        if nargin >= 4,
            isnormalized = varargin{3};
            if nargin >= 5,
                freqrange = varargin{4};
                if nargin >= 6,
                    centerdc = varargin{5};
                end
            end
        end
    end
end

validStrs = {'half','whole'};
if ~any(strcmpi(freqrange,validStrs)),
    error(message('signal:dspdata:abstractfreqresp:getdata:invalidFrequencyRangeStr', validStrs{ 1 }, validStrs{ 2 })); 
end

% [EOF]
