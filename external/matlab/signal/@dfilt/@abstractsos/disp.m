function disp(this)
%DISP Object display.
  
%   Author: V. Pellissier
%   Copyright 1988-2008 The MathWorks, Inc.

if length(this) > 1
    vectordisp(this);
    return;
end

fn = fieldnames(this);
N = length(fn);
% Reorder the fields. NumSamplesProcessed, ResetStates and States in
% the end.
if N>8, N=8; end
nidx = [3, 5:N, 1];
if this.PersistentMemory,
    % display states
    nidx = [nidx, 4];
end
fn = fn(nidx);

siguddutils('dispstr', this, fn, 24);

disp(this.filterquantizer, 24)

% [EOF]
