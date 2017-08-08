function designMethod_construct(d,varargin)
%DESIGNMETHOD_CONSTRUCT  'Real' constructor for the design method object.


%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Add filter types supported by this design method
addTypes(d,varargin{:});

% Call the listener here so it fires the first time
filterType_listener(d);







