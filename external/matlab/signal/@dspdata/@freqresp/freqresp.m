function this = freqresp(varargin)
%FREQRESP   Construct a FREQRESP object.

%   Author(s): J. Schickler
%   Copyright 2004 The MathWorks, Inc.

error(nargchk(0,8,nargin,'struct'));

this = dspdata.freqresp;

set(this, 'Name', 'Frequency Response');

% Construct a metadata object.
set(this,'Metadata',dspdata.powermetadata);
set(this.Metadata,...
    'FrequencyUnits','Hz',...
    'DataUnits','volts^2/Hz');

% Initialize Data and Frequencies with defaults or user specified values.
initialize(this,varargin{:});

% [EOF]
