function syncspecs(this, newspecs)
%SYNCSPECS Sync specs from the current specs to a new specs object.

%   Copyright 1999-2011 The MathWorks, Inc.

% Grab the old Fs information.
if ~isempty(this.CurrentSpecs),

    syncfs(this, newspecs);
    syncotherprops(this, newspecs);
end

% -------------------------------------------------------------------------
function syncfs(this, newspecs)

oldspecs = get(this, 'CurrentSpecs');
normalized = oldspecs.NormalizedFrequency;
if oldspecs.NormalizedFrequency,
    
    % If we are coming from a normalized setting, unnormalize, grab the fs
    % and renormalize.
    normalizefreq(oldspecs,false);
    fs = get(oldspecs, 'Fs');
    normalizefreq(oldspecs,true);
else
    fs = get(oldspecs, 'Fs');
end

% Use the fs from the oldspecs so that if the user unnormalizes they will
% get what they had set in Fs in the old specs.
normalizefreq(newspecs,false,fs);
normalizefreq(newspecs,normalized);

% -------------------------------------------------------------------------
function syncotherprops(this, newspecs)

oldspecs = get(this, 'CurrentSpecs');

% Sync all other props.
p = propstosync(newspecs);
prop_modified = false(length(p), 1);
for indx = 1:length(p),
    if isprop(oldspecs, p{indx}),
        set(newspecs, p{indx}, get(oldspecs, p{indx}));
        prop_modified(indx) = true;
    end
end

p(prop_modified) = [];
allspecs         = get(this, 'AllSpecs');

for indx = 1:length(p)
    for jndx = 1:length(allspecs)
        if isprop(allspecs(jndx), p{indx})
            % Objects in allspecs might have the properties that we want to
            % sync in the newspecs object. However, their values might not
            % be appropriate for the newspecs object settings. We do not
            % want to error out by setting an invalid value so we need to
            % catch errors and move on if we find one.
            try
              set(newspecs, p{indx}, get(allspecs(jndx), p{indx}));
              break
            catch ME %#ok<NASGU>
            end
        end
    end
end

% [EOF]
