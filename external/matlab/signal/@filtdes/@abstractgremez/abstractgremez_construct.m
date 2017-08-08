function abstractgremez_construct(h, varargin)
%ABSTRACTGREMEZ_CONSTRUCT

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

schema.prop(h, 'FIRType', 'gremezFIRType');
schema.prop(h, 'SinglePointBands', 'posint_vector');
schema.prop(h, 'ForcedFreqPoints', 'posint_vector');
schema.prop(h, 'IndeterminateFreqPoints', 'posint_vector');

dynMinOrder_construct(h, varargin{:});

% [EOF]
