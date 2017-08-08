function this = cpsd(varargin)
%CPSD   Construct a CPSD object.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(0,8,nargin,'struct'));

% Create object and set the properties specific to this object.
this = dspdata.cpsd;
set(this,'Name','Cross Power Spectral Density');

% Construct a metadata object.
set(this,'Metadata',dspdata.powermetadata);
set(this.Metadata,...
    'FrequencyUnits','Hz',...
    'DataUnits','volts^2/Hz');

% Initialize Data and Frequencies with defaults or user specified values.
initialize(this,varargin{:});

% [EOF]
