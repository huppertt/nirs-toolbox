function vectordisp(this)
%VECTORDISP   

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

for indx = 1:length(this),
    disp(class(this(indx)));
end

if strcmpi(get(0, 'formatspacing'), 'loose')
    fprintf(1, '\n');
end

% [EOF]
