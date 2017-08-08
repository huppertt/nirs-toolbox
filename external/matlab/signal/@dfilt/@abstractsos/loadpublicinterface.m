function loadpublicinterface(this, s)
%LOADPUBLICINTERFACE   Load the public interface.

%   Copyright 2008 The MathWorks, Inc.

abstract_loadpublicinterface(this, s);

if s.version.number > 3
   set(this, 'OptimizeScaleValues', s.OptimizeScaleValues); 
end

% [EOF]
