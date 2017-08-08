function R = validatespecs(this)
%VALIDATESPECS Validate the specs

%   Copyright 2011 The MathWorks, Inc.

super_validatespecs(this);

for idx = 1:this.Nbands
  R(1,idx) = this.(sprintf('%s%d%s','B',idx,'Ripple')); %#ok<AGROW>
end

% [EOF]
