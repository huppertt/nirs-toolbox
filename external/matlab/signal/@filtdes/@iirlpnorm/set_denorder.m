function dummy = set_denorder(h,M)
%SET_DENORDER Set the denominator order property.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

error(nargchk(2,2,nargin,'struct'));

% Get handle to numDenFilterOrder object
g = get(h,'numDenFilterOrderObj');

% Set value
set(g,'denOrder',M);

% Return a dummy value for now
dummy = 20;
