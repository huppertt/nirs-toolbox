function setmaskspecs(this)
%SETMASKSPECS   Set mask specs.

%   Copyright 2009 The MathWorks, Inc.

switch lower(this.WeightingType)
    case 'cmessage'
        setcmessageweightingmask(this);
    case 'itut041'
        setitut041weightingmask(this);
    case 'itur4684'
        setitur4684weightingmask(this);
end

% [EOF]
