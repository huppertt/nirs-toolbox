function specs = whichspecs(h)
%WHICHSPECS Determine which specs are required for this class.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

% Prop name, data type, default value, listener callback
specs(1) = cell2struct({'MagnitudeVector','double_vector',...
        [0 1 0],[],'magspec'},specfields(h),2);

specs(2) = cell2struct({'UpperVector','double_vector',...
        [.01 1.02 .01],[],'magspec'},specfields(h),2);

specs(3) = cell2struct({'LowerVector','double_vector',...
        [-.01 .98 -.01],[],'magspec'},specfields(h),2);

% [EOF]
