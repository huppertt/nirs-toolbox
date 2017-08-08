function n = thisnstates(this)
%THISNSTATES   

%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.

if isempty(this.refallpasscoeffs),
    n = 0;
else
    % Initialize with uppermost states
    n = length(this.refallpasscoeffs{1});

    % Add the states of each section
    for k = 2:length(this.refallpasscoeffs),
        % States may depend on previous section
        n = n + max(length(this.refallpasscoeffs{k}),length(this.refallpasscoeffs{k-1}));
    end

    % Add the states of the last section again since this is not shared with
    % any other section
    n = n + length(this.refallpasscoeffs{end});
end

% [EOF]
