function setmaskspecs(this)
%SETMASKSPECS   Set mask specs.

%   Copyright 2009 The MathWorks, Inc.

switch lower(this.WeightingType)
    case 'a'
        setaweightingmask(this);                
    case 'c'        
        setcweightingmask(this);
end

% [EOF]
