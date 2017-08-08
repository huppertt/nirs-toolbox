function this = freqz(varargin)
%FREQZ   Discrete-time frequency response object.

%   Author(s): P. Pacheco
%   Copyright 1988-2003 The MathWorks, Inc.

error(nargchk(0,8,nargin,'struct'));

% Create object and set the properties specific to this object.
this = dspdata.freqz;
set(this,'Name','Transfer Function Estimate');

% Construct a metadata object.
set(this,'Metadata',dspdata.powermetadata);
set(this.Metadata,'FrequencyUnits','Hz');
set(this.Metadata,'DataUnits','');

% Initialize Data and Frequencies with defaults or user specified values.
initialize(this,varargin{:});

%--------------------------------------------------------------------------
function validate(spectrumRange)
% This error checking should be done in the object's set method, but for
% enum datatypes UDD first checks the list before calling the set method.

validStrs = {'half','whole'};
if ~ischar(spectrumRange) | ~any(strcmpi(spectrumRange,validStrs)),
    error(message('signal:dspdata:freqz:freqz:invalidSpectrumRange', 'SpectrumRange', validStrs{ 1 }, validStrs{ 2 }));
end

% [EOF]
