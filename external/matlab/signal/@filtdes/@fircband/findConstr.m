function consStr = findConstr(h,ft,orderMode)
%FINDCONSTR Find the appropriate constructor.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

% This method is overloaded so that it can intercept 'minimum even' and
% 'minimum odd' OrderModes.  We can probably change this to:

% if nargin > 2 && any(strcmpi(orderMode, 'minimum even', 'minimum odd')),
%     orderMOde = 'minimum';
% end
% consStr = super::findConstr(h, ft, orderMode);

if nargin < 3, orderMode = 'specify'; end

indx = findConstrIndx(h,ft);

s = get(h,'availableTypes');

cb = sort(h.ConstrainedBands); % make sure the bands are in ascending order.
cbfield = sprintf('c%s', sprintf('%d', cb));
consStr = [];

if isfield(s(indx).construct, cbfield),
    consStr = s(indx).construct.(cbfield);
end

if strcmpi(orderMode, 'minimum')
    consStr = [];
elseif isempty(consStr),
    consStr = s(indx).construct.(orderMode);
end

% [EOF]
