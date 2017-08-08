function numden_construct(d,varargin)
%NUMDEN_CONSTRUCT  Base constructor for all subclasses.


%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.


% Construct a numDenFilterOrder object and store it
h = filtdes.numDenFilterOrder;
set(d,'numDenFilterOrderObj',h);


% Call super's constructor
designMethodwFs_construct(d,varargin{:});

