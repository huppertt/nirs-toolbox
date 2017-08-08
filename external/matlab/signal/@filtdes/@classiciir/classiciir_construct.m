function classiciir_construct(this)
%CLASSICIIR_CONSTRUCT   Constructor for the classiciir

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

schema.prop(this, 'MatchExactly', getmatchexactlytype(this));

set(this, 'MatchExactly', 'passband');

dynMinOrder_construct(this);

% [EOF]
