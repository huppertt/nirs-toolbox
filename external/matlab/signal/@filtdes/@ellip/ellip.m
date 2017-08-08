function d = ellip
%ELLIP  Constructor for this design method object.
%
%   Outputs:
%       d - Handle to the design method object

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

error(nargchk(0,3,nargin,'struct'));

d = filtdes.ellip;

% Call super's constructor
classiciir_construct(d);

% Set the tag
set(d,'Tag','Elliptic', 'MatchExactly', 'both');

