function filterorder = set_filterorder(this, filterorder)
%SET_FILTERORDER PreSet function for the 'filterorder' property.

%   Copyright 2005-2013 The MathWorks, Inc.

if rem(filterorder,2),
  if isfromdesignfilt(this)
    error(message('signal:fspecs:abstractdiffordmb:set_filterorder:InvalidOrderDesignfilt'));
  else
    error(message('signal:fspecs:abstractdiffordmb:set_filterorder:InvalidOrder'));
  end      
end


% [EOF]
