function varargout = design(d)
%Design  Method to design the filter given the specs.

%   Author(s): R. Losada, J. Schickler
%   Copyright 1988-2008 The MathWorks, Inc.

validate(d);

% Normalize frequencies temporarily
freqUnits = get(d,'freqUnits');
freqUnitsOpts = set(d,'freqUnits');
set(d,'freqUnits',freqUnitsOpts{1}); % Normalized (0 to 1)

try
    
    % Call class specific design
    Hd = thisdesign(d);

catch ME
    
    % Set frequencies back to what they were
    set(d,'freqUnits',freqUnits);
    throw(ME);
end

% Set frequencies back to what they were
set(d,'freqUnits',freqUnits);


% Add the dynamic property to store the mask commands
p = schema.prop(Hd, 'MaskInfo', 'MATLAB array');
set(p, 'Visible', 'Off');
set(Hd, 'MaskInfo', maskinfo(d));

if nargout,
    
    varargout = {Hd};
else
    
    % Eventually we will just call fvtool on Hd with masks on.
    drawmasknresp(d, Hd);
end

% [EOF]
