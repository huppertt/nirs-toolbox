function isincprops(c)
%ISINCPROPS Add the Inverse Sinc properties.

%   Copyright 2005-2011 The MathWorks, Inc.

p = schema.prop(c, 'FrequencyFactor', 'double');
set(p, 'FactoryValue', .5);

p = schema.prop(c, 'Power', 'double');
set(p, 'FactoryValue', 1);

p = schema.prop(c, 'CICRateChangeFactor', 'double');
set(p, 'FactoryValue', 1);

% [EOF]
