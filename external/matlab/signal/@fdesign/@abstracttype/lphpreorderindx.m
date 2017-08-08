function reorderindx = lphpreorderindx(this,Hd,nsections)
%LPHPREORDERINDX  Determine indices for lowpass/highpass reorder. 

%   Author(s): R. Losada
%   Copyright 1999-2004 The MathWorks, Inc.

% Form indices. Use nsections passed in, in case used for bandpass/bandstop
indx = 1:nsections;

% Take odd indices
oddindx = indx(1:2:end);
reorderodd = partialreorder(oddindx);

% Take even indices
evenindx = indx(2:2:end);
reordereven = partialreorder(evenindx);

reorderindx = [reorderodd,reordereven];

%---------------------------------------------------------------------
function reorderindx = partialreorder(indx)

% Initialize
reorderindx = zeros(size(indx));

for n = 1:ceil(length(indx)/2),
    reorderindx(end-n+1) = indx(2*n-1);    
end
for n = 1:floor(length(indx)/2),
    reorderindx(n) = indx(2*n);    
end


% [EOF]
