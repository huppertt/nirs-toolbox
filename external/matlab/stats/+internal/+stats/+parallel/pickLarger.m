function val = pickLarger(val, update)
%PICKLARGER is an argmax reduction operator.
%
%   PICKLARGER is an internal utility and is not meant for
%   general purpose use. Its functionality may change and should not be
%   relied upon.

%   Copyright 2010-2014 The MathWorks, Inc.

if isempty(val) || update{1} > val{1}
    val = update;
    return;
end
if update{1} == val{1}
    if update{2} > val{2}
        val = update;
    end
end

