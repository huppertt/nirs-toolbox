function normalizefreq(this,varargin)
%NORMALIZEFREQ   Normalize/un-normalize the frequency of the data object.

%   Author(s): P. Pacheco
%   Copyright 1988-2003 The MathWorks, Inc.

[normFlag,Fs] = parseinputs(this,varargin{:});

freq = this.Frequencies;
oldFs = getfs(this);  % Cache Fs stored in the object before it gets updated.
newFsFlag = false;

% If already in the units requested, and Fs hasn't changed return early.
if ~xor(this.NormalizedFrequency, normFlag), 
    % Only proceed if user specified a different Fs.
    if isequal(oldFs,Fs),
        return;
    else
        % Convert to normalized frequency in order to scale by new Fs.
        newFsFlag = true;
        freq = freq/oldFs*(2*pi);      
    end
end

if normFlag,    freq = freq/Fs*(2*pi);   % Convert to normalized frequency.
else            freq = freq/(2*pi)*Fs;   % Convert to linear frequency.
end

if normFlag,
    this.Fs = Fs; % Set Fs first since you can't do it after it's in normalized mode.
    this.privNormalizedFrequency = normFlag; 
else
    this.privNormalizedFrequency = normFlag;  % Change to linear to allow us to set Fs.
    this.Fs = Fs;
end
this.Frequencies = freq;

% Allow concrete classes to do further manipulation of the data if necessary.
thisnormalizefreq(this,oldFs,newFsFlag);

%--------------------------------------------------------------------------
function [normFlag,Fs] = parseinputs(this,varargin)
% Parse and validate inputs.

% Setup defaults
normFlag = true;
Fs = getfs(this);

if nargin >= 2,
    normFlag = varargin{1};
    if nargin == 3,
        Fs = varargin{2};
    end
end
    
if nargin == 3 && normFlag
    error(message('signal:dspdata:abstractfreqresp:normalizefreq:invalidInputArgumentFs', 'Fs'));
end

if ~islogical(normFlag),
    error(message('signal:dspdata:abstractfreqresp:normalizefreq:invalidLogicalFlag'));
end

% [EOF]
