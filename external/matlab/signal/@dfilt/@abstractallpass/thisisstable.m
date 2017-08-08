function b = thisisstable(this)
%THISISSTABLE   True if the object is stable.

%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.

b = false; % Assume unstable by default

if isfir(this),
    
    % Section is FIR, always stable
    b = true;
    
elseif signalpolyutils('isstable',[1,this.AllpassCoefficients]),

    b = true;
end

% [EOF]
