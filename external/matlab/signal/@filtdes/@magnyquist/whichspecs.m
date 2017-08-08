function specs = whichspecs(h)
%WHICHSPECS Determine which specs are required for this class.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.


specs(1) = cell2struct({'DesignType','magnyquistDesignType','Normal',...
        [],'magspec'},specfields(h),2);

% specs(end+1) = cell2struct({'Decay','udouble',0,...
%         [],'magspec'},specfields(h),2);

