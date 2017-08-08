function dopts = designoptions(this, method, sigonlyflag)
%DESIGNOPTIONS Return the design options.

%   Copyright 2004-2011 The MathWorks, Inc.

if nargin == 3
  hd = feval(getdesignobj(this, method, sigonlyflag));
  dopts = designopts(this, method, sigonlyflag);
else
  hd = feval(getdesignobj(this, method));
  dopts = designopts(this, method);
end

% Add the SystemObject property if structures are supported for this object
addsysobjdesignopt(hd);

if isprop(hd,'SystemObject')
  dopts.SystemObject = hd.SystemObject;
end

if isempty(fieldnames(dopts))
    return;
end

fn = setdiff(fieldnames(dopts), 'FilterStructure');

% The 'FilterStructure' property is not an enumerated value, so we need to
% get the options from the GETVALIDSTRUCTS method.
if isfield(dopts,'FilterStructure')
  dopts.DefaultFilterStructure = dopts.FilterStructure;
  dopts.FilterStructure        = getvalidstructs(hd);
end

% Loop over each of the fields to fix any problems.
for indx = 1:length(fn)
    
    dopts.(sprintf('Default%s', fn{indx})) = dopts.(fn{indx});
    dopts.(fn{indx}) = set(hd, fn{indx});
    
    if isempty(dopts.(fn{indx}))
        
        % If the field is empty we do not have an enumerated type, so we
        % need to get the valid values from the DataType of the property.
        p = findprop(hd, fn{indx});
        dopts.(fn{indx}) = get(p, 'DataType');
    elseif size(dopts.(fn{indx}), 2) == 1
        
        % Make sure that the strings show up by making the cell array a row.
        dopts.(fn{indx}) = transpose(dopts.(fn{indx}));
    end
end

% [EOF]
