function bool = isminord(d)
%ISMINORD Returns true if the object is minimum order.

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

try
    bool = strcmpi(d.ordermode, 'minimum');
catch
    
    % The ordermode property is private or not there when the object is not
    % minimum order, so the above get will fail.
    bool = false;
end    

% [EOF]
