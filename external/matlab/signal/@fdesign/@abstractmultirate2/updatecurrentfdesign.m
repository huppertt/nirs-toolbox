function updatecurrentfdesign(this)
%UPDATECURRENTFDESIGN   Update the current FDesign object.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

% Get the constructor for the current specification type.
cFDesignCon = getconstructor(this);

% If the CurrentSpecs is already correct, just return.
if ~strcmpi(class(this.CurrentFDesign), cFDesignCon),
    % If there are any stored SPEC objects see if our constructor matches.
    allFDesign = get(this, 'AllFDesign');
    if isempty(allFDesign), cFDesign = [];
    else,                   cFDesign = find(allFDesign, '-class', cFDesignCon); end

    % If we could not find the needed spec object, create it and store it.
    if isempty(cFDesign),
        cFDesign = feval(cFDesignCon);
        
        % Add smart defaults.
        multiratedefaults(cFDesign, max(getratechangefactors(this)));
        
        set(this, 'AllFDesign', [allFDesign; cFDesign]);
    end

    % Set the current specs, this will fire the pre-set to update the props.
    set(this, 'CurrentFDesign', cFDesign);
    l = [handle.listener(cFDesign, 'FaceChanged', {@facechanged_listener, cFDesign}); ...
        handle.listener(cFDesign, 'FaceChanging', {@facechanging_listener, cFDesign})];

    set(l, 'CallbackTarget', this);

    set(this, 'SpecificationTypeListeners', l);
end

updatefdesignfactors(this);

% -------------------------------------------------------------------------
function facechanging_listener(this, eventData, cFDesign)

rmprops(this, cFDesign);

% -------------------------------------------------------------------------
function facechanged_listener(this, eventData, cFDesign)

addprops(this, cFDesign);

% [EOF]
