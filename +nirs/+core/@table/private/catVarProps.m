function prop = catVarProps(aprop,bprop,avars,bvars)
% Concatenate per-variable cellstr properties

%   Copyright 2012 The MathWorks, Inc. 

% If both property values are not empty, concatenate them together.  If
% only one is not empty, concatenate it with empty values.  If both are
% empty, return empty.
if ~isempty(aprop) || ~isempty(bprop)
    nvars = length(avars) + length(bvars);
    prop = repmat({''},1,nvars);
    if ~isempty(aprop)
        prop(1:length(avars)) = aprop(avars);
    end
    if ~isempty(bprop)
        prop((length(avars)+1):nvars) = bprop(bvars);
    end
else
    prop = {};
end
