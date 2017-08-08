function disp(this)
%DISP Object display.
  
%   Author: V. Pellissier
%   Copyright 1988-2005 The MathWorks, Inc.

if length(this) > 1
    vectordisp(this);
    return;
end

fn = fieldnames(this);
N = length(fn);
% Reorder the fields. NumSamplesProcessed, ResetStates and States in
% the end.
nidx = [3, 5, 1];
if N>5, 
    nidx = [3, 6, 5, 1];
end

if this.PersistentMemory,
    % display states
    nidx = [nidx, 4];
end
fn = fn(nidx);

siguddutils('dispstr', this, fn, 20);

disp(this.filterquantizer, 20);

% [EOF]
