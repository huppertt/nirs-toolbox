function gremez_filterType_listener(h, varargin)
%GREMEZ_FILTERTYPE_LISTENER

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

super_filterType_listener(h, varargin{:});

% Turn all the properties on and turn them off according to the output of
% DISABLEDPROPS of the ResponseType.
enabdynprop(h, 'FIRType', 'On');
enabdynprop(h, 'SinglePointBands', 'On');
enabdynprop(h, 'ForcedFreqPoints', 'On');
enabdynprop(h, 'IndeterminateFreqPoints', 'On');

p = disabledprops(h.ResponseTypeSpecs);

indx = find(strcmpi('initorder', p));
if ~isempty(indx), p(indx) = []; end

for indx = 1:length(p)
    enabdynprop(h, p{indx}, 'Off');
end

% [EOF]
