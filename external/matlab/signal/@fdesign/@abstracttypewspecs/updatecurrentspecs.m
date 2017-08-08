function updatecurrentspecs(this)
%UPDATECURRENTSPECS   Update the currentSpecs object.

%   Author(s): J. Schickler
%   Copyright 1999-2005 The MathWorks, Inc.

% Get the constructor for the current specification type.
cSpecCon = getconstructor(this);

% % If the CurrentSpecs is already correct, just return.
% if strcmpi(class(this.CurrentSpecs), cSpecCon), return; end

% If there are any stored SPEC objects see if our constructor matches.
allSpecs = get(this, 'AllSpecs');
if isempty(allSpecs), cSpec = [];
else,                 cSpec = find(allSpecs, '-class', cSpecCon); end

% If we could not find the needed spec object, create it and store it.
if isempty(cSpec),
    cSpec = feval(cSpecCon);
    set(this, 'AllSpecs', [allSpecs; cSpec]);
end

% Set the current specs, this will fire the pre-set to update the props.
set(this, 'CurrentSpecs', cSpec);

% [EOF]
