function varargout = multistage(this, varargin)
%MULTISTAGE   Design a multistage equiripple decimator.

%   Author(s): J. Schickler
%   Copyright 2005-2014 The MathWorks, Inc.

% Make sure that MULTISTAGE is valid for this set of specifications.
if ~isdesignmethod(this, 'multistage')
    error(message('signal:fdesign:decimator:multistage:invalidMethod', 'MULTISTAGE', this.Specification));
end

% Parse the inputs for the multirate filter structure.
[filtstruct, varargin] = parsestruct(this, 'firdecim', 'multistage', varargin{:});

dfactor = get(this, 'DecimationFactor');

Hm = multistage(this.CurrentFDesign, varargin{:}, 'rcf', -dfactor,...
    'FilterStructure',filtstruct);

% Use THIS object instead of the 'CurrentFDesign'.
if isa(Hm, 'dsp.private.FilterAnalysis')
    setMetaData(Hm,this);
else
    Hm.setfdesign(this);
end

if nargout
    varargout = {Hm};
else
    if this.NormalizedFrequency,
        inputs = {'NormalizedFrequency', 'On'};
    else
        inputs = {'Fs', this.Fs};
    end

    fvtool(Hm, inputs{:});
end

% [EOF]
