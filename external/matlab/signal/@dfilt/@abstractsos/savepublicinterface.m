function s = savepublicinterface(this)
%SAVEPUBLICINTERFACE   Save the public interface.

%   Copyright 2008 The MathWorks, Inc.

s = abstract_savepublicinterface(this);
s.OptimizeScaleValues = get(this, 'OptimizeScaleValues');

% [EOF]
