function checkvalidparallel(this)
%CHECKVALIDPARALLEL   Check if parallel is valid and error if not.

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

if ~isvalidparallel(this),
    error(message('signal:dfilt:parallel:checkvalidparallel:invalidParallel'));
end

% [EOF]
