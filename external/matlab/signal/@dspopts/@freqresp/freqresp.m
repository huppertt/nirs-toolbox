function this = freqresp(varargin)
%FREQRESP   Construct a FREQRESP object.

%   Author(s): J. Schickler
%   Copyright 1988-2006 The MathWorks, Inc.

this = dspopts.freqresp;
this.NFFT = 8192;  % Avoid having NFFT set to "new" default 'Nextpow2'.

% We need to special case the constructor to allow users to specify the
% freequency vector without having to specify the frequencyspecification
% first.  This would be burdensome.

for indx = 1:3:length(varargin)
    if any(strncmpi(varargin{indx}, {'NFFT', 'FrequencyVector'}, length(varargin{indx})))
        varargin = {varargin{1:indx-1}, ...
            'FrequencySpecification', varargin{indx}, ...
            varargin{indx:end}};
    end
end

if nargin
    set(this, varargin{:});
end

% [EOF]
