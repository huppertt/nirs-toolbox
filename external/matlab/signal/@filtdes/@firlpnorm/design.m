function varargout = design(d)
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
N = get(d,'order');
I = get(d,'initNum');
[P,DENS] = getNumericSpecs(d);

args = {N,F,E,A,W,P,{DENS}};
if ~isempty(I),
    args = {args{:},I};
end

% Determine optional design args
minPhase = get(d,'minphase');

switch minPhase,
    case 'on',
        % Design minimum-phase filter
        b = firlpnorm(args{:},'minphase');
    case 'off'
        b = firlpnorm(args{:});
end

% Construct object
h = dfilt.dffir(b);

% Add the dynamic property to store the mask commands
p = schema.prop(h, 'MaskInfo', 'MATLAB array');
set(p, 'Visible', 'Off');
set(h, 'MaskInfo', maskinfo(d));

if nargout,
    
    varargout = {h};
else,
    drawmasknresp(d, h);
end

% [EOF]
