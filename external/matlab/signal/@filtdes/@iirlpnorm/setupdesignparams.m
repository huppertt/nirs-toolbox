function args = setupdesignparams(d)
%Design  Method to design the filter given the specs.

%   Author(s): R. Losada, J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

% Normalize frequencies temporarily
freqUnits = get(d,'freqUnits');
freqUnitsOpts = set(d,'freqUnits');
set(d,'freqUnits',freqUnitsOpts{1}); % Normalized (0 to 1)

% Get handle to the filter type object
hft = get(d,'responseTypeSpecs');
[F,E,A,W] = getNumericSpecs(hft,d);

% Set frequencies back to what they were
set(d,'freqUnits',freqUnits);

% Set up design params
N = get(d,'numOrder');
M = get(d,'denOrder');
IN = get(d,'initNum');
ID = get(d,'initDen');
[P,DENS] = getNumericSpecs(d);

args = {N,M,F,E,A,W,P,{DENS}};
if ~isempty(IN),
	args = {args{:},IN,ID};
end

% [EOF]
