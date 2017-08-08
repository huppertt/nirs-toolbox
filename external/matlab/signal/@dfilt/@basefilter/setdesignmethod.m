function setdesignmethod(this, designmethod)
%SETDESIGNMETHOD   Set the designmethod.

%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.

if ~ischar(designmethod),
    error(message('signal:dfilt:basefilter:setdesignmethod:MustBeAString'));
end

set(this, 'privdesignmethod', designmethod);

% [EOF]
