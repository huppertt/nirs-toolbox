function disp(this)
%DISP Object display.
  
%   Author: V. Pellissier
%   Copyright 1988-2009 The MathWorks, Inc.

if length(this) > 1
    vectordisp(this);
    return;
end

fn = fieldnames(this);

nidx = [3, 7, 5, 6, 1];


if this.PersistentMemory,
    % display states
    nidx = [nidx, 4];
end
fn = fn(nidx);

siguddutils('dispstr', this, fn, 20);

disp(this.filterquantizer, 20);

% [EOF]
