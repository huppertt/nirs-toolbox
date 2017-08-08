function p = props2normalize(this)
%PROPS2NORMALIZE   Return the property name to normalize.

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

p = [];
for i=1:this.NBands,
    p = [p {sprintf('%s%d%s','B',i,'Frequencies')}];
end

% [EOF]
